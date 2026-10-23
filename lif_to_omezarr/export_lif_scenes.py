import argparse
import csv
import json
import re
from pathlib import Path

try:
    import tifffile
except ImportError:
    tifffile = None

try:
    from bioio import BioImage
except ImportError:
    BioImage = None

from extract_lif_metadata import bioio_physical_size, bioio_scale_to_nm


def sanitize_name(value):
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "")).strip("_")
    return cleaned or "scene"


def scene_dataset_name(base_name, scene, index, total):
    if total <= 1:
        return sanitize_name(base_name)
    label = sanitize_name(scene.get("name") or f"s{index:02d}")
    return sanitize_name(f"{base_name}_{index:02d}_{label}")


def physical_size_um(metadata, axis):
    value = metadata.get(f"{axis.lower()}_scale_nm")
    if value is None:
        return None
    return float(value) / 1000.0


def is_ome_tiff(path):
    name = Path(path).name.lower()
    return name.endswith(".ome.tif") or name.endswith(".ome.tiff")


def merged_scene_keys(scenes):
    keys = set()
    for scene in scenes:
        match = re.fullmatch(r"(.*?)[ _]?merged", str(scene), flags=re.IGNORECASE)
        if match:
            keys.add(match.group(1).lower())
    return keys


def drop_unmerged_tiles(scenes):
    """Keep "R 1_Merged" and drop the raw "R 1" tiles it was stitched from."""
    merged = merged_scene_keys(scenes)
    return [scene for scene in scenes if str(scene).lower() not in merged]


def scene_scales_nm(image, metadata):
    """Pixel sizes of the current scene, falling back to the LIF-level values."""
    physical_pixel_sizes = getattr(image, "physical_pixel_sizes", None)
    scales = {}
    for axis in ("x", "y", "z"):
        value = bioio_scale_to_nm(bioio_physical_size(physical_pixel_sizes, axis.upper()))
        scales[f"{axis}_scale_nm"] = value if value is not None else metadata.get(f"{axis}_scale_nm")
    return scales


def scene_channels(image, metadata):
    """Channels of the current scene; scenes in one LIF can have different channel sets."""
    channels = list(metadata.get("channels") or [])
    names = [str(name) for name in (getattr(image, "channel_names", []) or [])]
    if not names or names == [str(channel.get("label")) for channel in channels]:
        return channels
    by_label = {str(channel.get("label")).lower(): channel for channel in channels}
    scene_rows = []
    for index, name in enumerate(names):
        channel = dict(by_label.get(name.lower(), {"label": name, "display": sanitize_name(name)}))
        channel["index"] = index
        scene_rows.append(channel)
    return scene_rows


def write_scene_ome_tiff(image, output_path, metadata):
    data = image.get_image_data("TCZYX")
    channel_names = [
        channel.get("display") or channel.get("label") or f"channel_{index}"
        for index, channel in enumerate(metadata.get("channels") or [])
    ]
    ome_metadata = {
        "axes": "TCZYX",
        "PhysicalSizeX": physical_size_um(metadata, "x"),
        "PhysicalSizeXUnit": "um",
        "PhysicalSizeY": physical_size_um(metadata, "y"),
        "PhysicalSizeYUnit": "um",
        "PhysicalSizeZ": physical_size_um(metadata, "z"),
        "PhysicalSizeZUnit": "um",
    }
    if channel_names:
        ome_metadata["Channel"] = {"Name": channel_names}

    tifffile.imwrite(
        output_path,
        data,
        ome=True,
        metadata=ome_metadata,
        bigtiff=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Export each scene of a LIF/CZI file (or pass through an OME-TIFF) to a scene-specific OME-TIFF."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--base-name", required=True)
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--scene-dir", required=True)
    parser.add_argument("--scene-metadata-dir", required=True)
    parser.add_argument("--scene-table", required=True)
    parser.add_argument(
        "--keep-unmerged-tiles",
        action="store_true",
        help="Also export raw tile scenes that have a *_Merged counterpart.",
    )
    args = parser.parse_args()

    metadata = json.loads(Path(args.metadata_json).read_text(encoding="utf-8"))
    passthrough = is_ome_tiff(args.input)
    if passthrough:
        # An OME-TIFF is already a single scene that eubi can read directly.
        image = None
        all_scenes = [Path(args.input).name]
    else:
        if BioImage is None:
            raise RuntimeError("bioio and its reader plugin (bioio-lif / bioio-czi) are required to export scenes")
        if tifffile is None:
            raise RuntimeError("tifffile is required to export scenes")
        image = BioImage(args.input)
        all_scenes = list(getattr(image, "scenes", []) or []) or ["scene_0"]
    scenes = all_scenes if args.keep_unmerged_tiles else drop_unmerged_tiles(all_scenes)
    skipped = [str(scene) for scene in all_scenes if scene not in scenes]
    if skipped:
        print(f"Skipping raw tiles with a merged version: {', '.join(skipped)}")

    scene_dir = Path(args.scene_dir)
    scene_metadata_dir = Path(args.scene_metadata_dir)
    scene_dir.mkdir(parents=True, exist_ok=True)
    scene_metadata_dir.mkdir(parents=True, exist_ok=True)

    scene_rows = []
    scene_lookup = {
        str(scene.get("name")): scene
        for scene in metadata.get("scenes") or []
        if scene.get("name") is not None
    }
    total = len(all_scenes)

    for index, scene_name in enumerate(all_scenes):
        if scene_name not in scenes:
            continue
        source_scene = scene_lookup.get(str(scene_name), {"index": index, "name": str(scene_name)})
        scene_dataset = scene_dataset_name(args.base_name, source_scene, index, total)
        scene_tiff = scene_dir / f"{scene_dataset}.ome.tif"
        scene_metadata = scene_metadata_dir / f"{scene_dataset}_metadata.json"

        scene_json = dict(metadata)
        if passthrough:
            scene_tiff.symlink_to(Path(args.input).resolve())
        else:
            image.set_scene(scene_name)
            scene_json.update(scene_scales_nm(image, metadata))
            scene_json["channels"] = scene_channels(image, metadata)
            scene_json["size_c"] = len(scene_json["channels"])
            write_scene_ome_tiff(image, scene_tiff, scene_json)

        scene_json["name"] = scene_dataset
        scene_json["parent_lif_name"] = args.base_name
        scene_json["raw_path"] = str(Path(args.input))
        scene_json["source_suffix"] = "".join(Path(args.input).suffixes[-2:]).lower() if passthrough else Path(args.input).suffix.lower()
        scene_json["selected_scene"] = str(scene_name)
        scene_json["selected_scene_index"] = index
        scene_json["scene_count"] = total
        scene_json["grid_index"] = index
        scene_json["scene_ome_tiff"] = str(scene_tiff)
        scene_metadata.write_text(
            json.dumps(scene_json, indent=2, sort_keys=True),
            encoding="utf-8",
        )

        scene_rows.append(
            {
                "dataset_name": scene_dataset,
                "scene_index": index,
                "scene_name": str(scene_name),
                "scene_ome_tiff": str(scene_tiff),
                "scene_metadata": str(scene_metadata),
                "x_scale_nm": scene_json.get("x_scale_nm") or "",
                "y_scale_nm": scene_json.get("y_scale_nm") or "",
                "z_scale_nm": scene_json.get("z_scale_nm") or "",
            }
        )

    with Path(args.scene_table).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset_name",
                "scene_index",
                "scene_name",
                "scene_ome_tiff",
                "scene_metadata",
                "x_scale_nm",
                "y_scale_nm",
                "z_scale_nm",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(scene_rows)


if __name__ == "__main__":
    main()
