#!/usr/bin/env python3
import argparse
import csv
import math
import re
from pathlib import Path

from readlif.reader import LifFile


def safe_image_name(name):
    return re.sub(r'[\\/*?:"<>|]', '_', str(name)).strip()


def safe_output_name(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name)).strip('_')


def nm_per_pixel(scale_px_per_um):
    if scale_px_per_um is None:
        return None
    scale_px_per_um = float(scale_px_per_um)
    if scale_px_per_um <= 0:
        return None
    return 1000.0 / scale_px_per_um


def req_mem_gb(width, height, channels, bytes_per_pixel=2):
    raw_bytes = max(1, int(width)) * max(1, int(height)) * max(1, int(channels)) * bytes_per_pixel
    # EuBI/BioFormats conversions need overhead beyond raw pixel bytes.
    return int(min(max(math.ceil(raw_bytes * 6 / (1024 ** 3)), 16), 128))


def lif_context(lif_path, raw_lif_dir, run_dir):
    try:
        rel_parts = lif_path.relative_to(run_dir).parts[:-1]
    except ValueError:
        rel_parts = lif_path.relative_to(raw_lif_dir).parts[:-1]
    run_name = run_dir.name
    run_date = ""
    condition = lif_path.stem

    for part in rel_parts:
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", part):
            run_date = part
    for part in reversed(rel_parts):
        if part.lower() in {"control", "treatment"}:
            condition = part
            break
    else:
        if lif_path.stem.lower() in {"control", "treatment"}:
            condition = lif_path.stem
        elif rel_parts:
            condition = rel_parts[-1]

    return run_name, run_date, condition


def tif_candidates(lif_path, raw_lif_dir, run_dir, condition, safe_name):
    candidates = [
        lif_path.parent / f"{safe_name}.tif",
        run_dir / condition / f"{safe_name}.tif",
    ]
    try:
        rel_parent = lif_path.parent.relative_to(raw_lif_dir)
        candidates.append(run_dir / rel_parent / f"{safe_name}.tif")
    except ValueError:
        pass
    seen = set()
    unique = []
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def first_existing_tif(candidates, run_dir, safe_name):
    for candidate in candidates:
        if candidate.exists():
            return candidate

    matches = sorted(run_dir.rglob(f"{safe_name}.tif")) + sorted(run_dir.rglob(f"{safe_name}.tiff"))
    if matches:
        return matches[0]

    return candidates[0]


def row_for_image(lif_path, raw_lif_dir, run_dir, extracted_tif_dir, default_z_scale_nm, image):
    run_name, run_date, condition = lif_context(lif_path, raw_lif_dir, run_dir)
    safe_name = safe_image_name(image.name)
    if extracted_tif_dir is not None:
        tif_path = extracted_tif_dir / run_name / run_date / condition / f"{safe_name}.tif"
    else:
        tif_path = first_existing_tif(tif_candidates(lif_path, raw_lif_dir, run_dir, condition, safe_name), run_dir, safe_name)
    output_name = safe_output_name(f"{run_name}_{condition}_{safe_name}")

    scale_x = image.scale[0] if len(image.scale) > 0 else None
    scale_y = image.scale[1] if len(image.scale) > 1 else None
    scale_z = image.scale[2] if len(image.scale) > 2 else None

    x_scale_nm = nm_per_pixel(scale_x)
    y_scale_nm = nm_per_pixel(scale_y)
    z_scale_nm = nm_per_pixel(scale_z) if scale_z is not None else default_z_scale_nm

    width = int(image.dims.x)
    height = int(image.dims.y)
    z_size = int(image.dims.z)
    channels = int(image.channels)

    return {
        "run_name": run_name,
        "run_date": run_date,
        "condition": condition,
        "lif_path": str(lif_path),
        "image_name": image.name,
        "safe_name": safe_name,
        "tif_path": str(tif_path),
        "tif_exists": str(tif_path.exists()).upper(),
        "output_name": output_name,
        "x_scale_nm": f"{x_scale_nm:.10g}" if x_scale_nm is not None else "",
        "y_scale_nm": f"{y_scale_nm:.10g}" if y_scale_nm is not None else "",
        "z_scale_nm": f"{z_scale_nm:.10g}" if z_scale_nm is not None else "",
        "scale_x_px_per_um": f"{float(scale_x):.10g}" if scale_x is not None else "",
        "scale_y_px_per_um": f"{float(scale_y):.10g}" if scale_y is not None else "",
        "scale_z_px_per_um": f"{float(scale_z):.10g}" if scale_z is not None else "",
        "width_px": str(width),
        "height_px": str(height),
        "z_size": str(z_size),
        "channels": str(channels),
        "req_mem": str(req_mem_gb(width, height, channels)),
    }


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Build osFISH TIFF conversion table from LIF metadata.")
    parser.add_argument("--raw-lif-dir", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--all-output", default="all_images.tsv", type=Path)
    parser.add_argument("--process-output", default="images_to_process.tsv", type=Path)
    parser.add_argument("--default-z-scale-nm", default=1000.0, type=float)
    parser.add_argument("--extracted-tif-dir", default="", type=Path)
    args = parser.parse_args()

    extracted_tif_dir = args.extracted_tif_dir if str(args.extracted_tif_dir) else None

    lif_paths = sorted(args.raw_lif_dir.rglob("*.lif"))
    if not lif_paths:
        raise SystemExit(f"No .lif files found in {args.raw_lif_dir}")

    rows = []
    for lif_path in lif_paths:
        lif = LifFile(str(lif_path))
        for image in lif.get_iter_image():
            rows.append(row_for_image(lif_path, args.raw_lif_dir, args.run_dir, extracted_tif_dir, args.default_z_scale_nm, image))

    fieldnames = [
        "run_name", "run_date", "condition", "lif_path", "image_name", "safe_name", "tif_path", "tif_exists",
        "output_name", "x_scale_nm", "y_scale_nm", "z_scale_nm",
        "scale_x_px_per_um", "scale_y_px_per_um", "scale_z_px_per_um",
        "width_px", "height_px", "z_size", "channels", "req_mem",
    ]
    write_tsv(args.all_output, rows, fieldnames)

    process_rows = [
        row for row in rows
        if row["tif_exists"] == "TRUE" and row["x_scale_nm"] and row["y_scale_nm"] and row["z_scale_nm"]
    ]
    write_tsv(args.process_output, process_rows, fieldnames)

    if not process_rows:
        raise SystemExit("No TIFF files with complete scale metadata were found for conversion.")


if __name__ == "__main__":
    main()
