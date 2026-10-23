import argparse
import csv
import json
import re
from pathlib import Path


DEFAULT_COLORS = ["red", "green", "yellow", "blue", "magenta", "cyan", "white"]


def image_uri(data_dir, name, uri_from, uri_to):
    """Path of the OME-Zarr as seen by the MoBIE machine, e.g. /g/schwab/... -> A:\\..."""
    path = f"{str(data_dir).rstrip('/')}/{name}.ome.zarr"
    if not uri_from:
        return path
    prefix = uri_from.rstrip("/")
    if path != prefix and not path.startswith(f"{prefix}/"):
        raise ValueError(f"{path} is not below --uri-from {uri_from}")
    rest = path[len(prefix):].lstrip("/")
    separator = "\\" if "\\" in uri_to else "/"
    return uri_to.rstrip("/\\") + separator + rest.replace("/", separator)


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
                "contrast_limits": "",
                "excitation_wavelength_nm": "",
                "emission_wavelength_nm": "",
                "emission_begin_nm": "",
                "emission_end_nm": "",
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


def load_metadata_rows(paths):
    rows = []
    for path in paths:
        metadata_path = Path(path)
        if metadata_path.is_dir():
            files = sorted(metadata_path.glob("*_metadata.json"))
        else:
            files = [metadata_path]
        for fallback_index, file_path in enumerate(files):
            metadata = json.loads(file_path.read_text(encoding="utf-8"))
            metadata.setdefault("grid_index", fallback_index)
            rows.append(metadata)
    return rows


def load_inventory(path):
    if not path:
        return {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return {row["dataset_name"]: row for row in csv.DictReader(handle, delimiter="\t")}


def run_name(metadata, inventory_row, fallback):
    return (inventory_row or {}).get("run") or metadata.get("parent_lif_name") or fallback


def assign_grid_positions(metadata_rows, inventory, fallback_view):
    """Lay out all scenes of one run on one near-square grid, ordered by LIF and scene."""
    by_run = {}
    for metadata in metadata_rows:
        inventory_row = inventory.get(metadata.get("parent_lif_name"))
        by_run.setdefault(run_name(metadata, inventory_row, fallback_view), []).append(metadata)
    positions = {}
    for run, rows in by_run.items():
        rows.sort(key=lambda row: (str(row.get("parent_lif_name")), int(row.get("grid_index", 0))))
        columns = near_square_columns(len(rows))
        for index, metadata in enumerate(rows):
            positions[metadata["name"]] = (run, f"({index % columns},{index // columns})")
    return positions


def main():
    parser = argparse.ArgumentParser(
        description="Create one MoBIE collection table, one grid view per run, from extracted LIF scene metadata."
    )
    parser.add_argument("--metadata-json", required=True, nargs="+")
    parser.add_argument("--inventory", default="", help="Inventory TSV from build_run_inventory.py.")
    parser.add_argument("--dataset-name", default="osFISH", help="View name when no inventory is given.")
    parser.add_argument("--data-dir", required=True, help="Directory holding the converted OME-Zarr files.")
    parser.add_argument("--uri-from", default="", help="Path prefix to replace, e.g. /g/schwab.")
    parser.add_argument("--uri-to", default="", help="Replacement prefix, e.g. 'A:\\'.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    metadata_rows = load_metadata_rows(args.metadata_json)
    inventory = load_inventory(args.inventory)
    positions = assign_grid_positions(metadata_rows, inventory, args.dataset_name)

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
        "run",
        "run_date",
        "source_file",
        "format",
        "scene",
        "sample",
        "condition",
    ]

    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for metadata in metadata_rows:
            dataset_name = metadata.get("name") or args.dataset_name
            inventory_row = inventory.get(metadata.get("parent_lif_name")) or {}
            run, grid_position = positions[dataset_name]
            uri = image_uri(args.data_dir, dataset_name, args.uri_from, args.uri_to)
            for channel in metadata_channels(metadata):
                writer.writerow(
                    {
                        "uri": uri,
                        "name": f"{dataset_name}_c{channel['index']}_{channel['display']}",
                        "view": run,
                        "grid": run,
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
                        "run": run,
                        "run_date": inventory_row.get("run_date", ""),
                        "source_file": inventory_row.get("source_file") or Path(str(metadata.get("raw_path", ""))).name,
                        "format": inventory_row.get("format") or str(metadata.get("source_suffix", "")).lstrip("."),
                        "scene": metadata.get("selected_scene", ""),
                        "sample": inventory_row.get("sample", ""),
                        "condition": inventory_row.get("condition", ""),
                    }
                )


if __name__ == "__main__":
    main()
