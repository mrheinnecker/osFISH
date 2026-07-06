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


def write_scene_ome_tiff(image, scene_name, output_path, metadata):
    image.set_scene(scene_name)
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
        description="Export each scene in a Leica LIF file to a scene-specific OME-TIFF."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--base-name", required=True)
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--scene-dir", required=True)
    parser.add_argument("--scene-metadata-dir", required=True)
    parser.add_argument("--scene-table", required=True)
    args = parser.parse_args()

    if BioImage is None:
        raise RuntimeError("bioio and bioio-lif are required to export LIF scenes")
    if tifffile is None:
        raise RuntimeError("tifffile is required to export LIF scenes")

    metadata = json.loads(Path(args.metadata_json).read_text(encoding="utf-8"))
    image = BioImage(args.input)
    scenes = list(getattr(image, "scenes", []) or [])
    if not scenes:
        scenes = ["scene_0"]

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
    total = len(scenes)

    for index, scene_name in enumerate(scenes):
        source_scene = scene_lookup.get(str(scene_name), {"index": index, "name": str(scene_name)})
        scene_dataset = scene_dataset_name(args.base_name, source_scene, index, total)
        scene_tiff = scene_dir / f"{scene_dataset}.ome.tif"
        scene_metadata = scene_metadata_dir / f"{scene_dataset}_metadata.json"

        write_scene_ome_tiff(image, scene_name, scene_tiff, metadata)

        scene_json = dict(metadata)
        scene_json["name"] = scene_dataset
        scene_json["parent_lif_name"] = args.base_name
        scene_json["raw_path"] = str(Path(args.input))
        scene_json["source_suffix"] = ".lif"
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
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(scene_rows)


if __name__ == "__main__":
    main()
