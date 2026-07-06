import argparse
import csv
import json
import re
from pathlib import Path


DEFAULT_COLORS = ["red", "green", "yellow", "blue", "magenta", "cyan", "white"]


def public_s3_prefix(bucket):
    bucket_path = re.sub(r"^s3embl/", "", bucket).rstrip("/")
    return f"https://s3.embl.de/{bucket_path}"


def compact(value):
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def color_for_display(display, fallback):
    value = compact(display)
    if value in {"gray", "grey"}:
        return "white"
    if value in {"magenta", "red"}:
        return "magenta"
    if value == "green" or "gfp" in value or "fitc" in value:
        return "cyan"
    if value == "tl":
        return "white"
    if "chloa" in value or "chlorophyll" in value or "cy5" in value:
        return "magenta"
    if "cy3" in value or "pe" in value:
        return "yellow"
    return fallback


def sanitize_display(value, index):
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "")).strip("_")
    return cleaned or f"channel_{index}"


def metadata_channels(metadata):
    channels = metadata.get("channels") or []
    if not channels:
        return [
            {
                "index": 0,
                "label": "channel_0",
                "display": "channel_0",
                "color": DEFAULT_COLORS[0],
            }
        ]

    rows = []
    seen = set()
    for fallback_index, channel in enumerate(channels):
        index = int(channel.get("index", fallback_index))
        label = str(channel.get("label") or f"channel_{index}")
        display = sanitize_display(channel.get("display") or label, index)
        if display.lower() in seen:
            continue
        seen.add(display.lower())
        color = color_for_display(
            display,
            str(channel.get("color") or DEFAULT_COLORS[index % len(DEFAULT_COLORS)]),
        )
        rows.append(
            {
                "index": len(rows),
                "label": label,
                "display": display,
                "color": color,
                "contrast_limits": str(channel.get("contrast_limits") or ""),
                "excitation_wavelength_nm": str(channel.get("excitation_wavelength_nm") or ""),
                "emission_wavelength_nm": str(channel.get("emission_wavelength_nm") or ""),
                "emission_begin_nm": str(channel.get("emission_begin_nm") or ""),
                "emission_end_nm": str(channel.get("emission_end_nm") or ""),
            }
        )
    return rows


def near_square_columns(n):
    if n <= 1:
        return 1
    columns = 1
    while columns * columns < n:
        columns += 1
    return columns


def load_metadata_rows(path):
    metadata_path = Path(path)
    if metadata_path.is_dir():
        files = sorted(metadata_path.glob("*_metadata.json"))
    else:
        files = [metadata_path]
    rows = []
    for fallback_index, file_path in enumerate(files):
        metadata = json.loads(file_path.read_text(encoding="utf-8"))
        metadata.setdefault("grid_index", fallback_index)
        rows.append(metadata)
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Create a one-dataset MoBIE collection table from extracted LIF metadata."
    )
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--dataset-name", required=True)
    parser.add_argument("--s3-bucket", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    metadata_rows = load_metadata_rows(args.metadata_json)
    grid_columns = near_square_columns(len(metadata_rows))

    fieldnames = [
        "uri",
        "name",
        "view",
        "grid",
        "grid_position",
        "channel",
        "display",
        "color",
        "excitation_wavelength_nm",
        "emission_wavelength_nm",
        "emission_begin_nm",
        "emission_end_nm",
        "contrast_limits",
        "blend",
        "format",
        "exclusive",
    ]

    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for metadata in metadata_rows:
            dataset_name = metadata.get("name") or args.dataset_name
            grid_index = int(metadata.get("grid_index", 0))
            grid_position = f"({grid_index % grid_columns},{grid_index // grid_columns})"
            uri = f"{public_s3_prefix(args.s3_bucket)}/{dataset_name}.ome.zarr/"
            for channel in metadata_channels(metadata):
                writer.writerow(
                    {
                        "uri": uri,
                        "name": f"{dataset_name}_c{channel['index']}_{channel['display']}",
                        "view": "osFISH",
                        "grid": "osFISH",
                        "grid_position": grid_position,
                        "channel": channel["index"],
                        "display": channel["display"],
                        "color": channel["color"],
                        "excitation_wavelength_nm": channel["excitation_wavelength_nm"],
                        "emission_wavelength_nm": channel["emission_wavelength_nm"],
                        "emission_begin_nm": channel["emission_begin_nm"],
                        "emission_end_nm": channel["emission_end_nm"],
                        "contrast_limits": channel["contrast_limits"],
                        "blend": "sum",
                        "format": "OmeZarr",
                        "exclusive": "TRUE",
                    }
                )


if __name__ == "__main__":
    main()
