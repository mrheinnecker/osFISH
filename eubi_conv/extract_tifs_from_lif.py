#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

import numpy as np
import tifffile
from readlif.reader import LifFile


def safe_image_name(name):
    return re.sub(r'[\\/*?:"<>|]', '_', str(name)).strip()


def write_tiff_for_image(image, out_path):
    channels = []
    for channel in range(int(image.channels)):
        arr = np.asarray(image.get_frame(c=channel))
        channels.append(arr)
    out_arr = np.stack(channels, axis=0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        out_path,
        out_arr,
        imagej=True,
        metadata={"axes": "CYX"},
    )


def main():
    parser = argparse.ArgumentParser(description="Extract CYX TIFF files from one Leica LIF file.")
    parser.add_argument("--lif", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--manifest", default="extracted_tifs.tsv", type=Path)
    parser.add_argument("--overwrite", default="FALSE")
    args = parser.parse_args()

    overwrite = str(args.overwrite).strip().lower() in {"true", "1", "yes", "y"}
    lif = LifFile(str(args.lif))
    rows = []

    for index, image in enumerate(lif.get_iter_image()):
        safe_name = safe_image_name(image.name)
        out_path = args.outdir / f"{safe_name}.tif"
        if overwrite or not out_path.exists():
            write_tiff_for_image(image, out_path)
        rows.append({
            "lif_path": str(args.lif),
            "image_index": index,
            "image_name": image.name,
            "safe_name": safe_name,
            "tif_path": str(out_path),
            "channels": int(image.channels),
            "width_px": int(image.dims.x),
            "height_px": int(image.dims.y),
        })

    with open(args.manifest, "w", newline="", encoding="utf-8") as handle:
        fieldnames = ["lif_path", "image_index", "image_name", "safe_name", "tif_path", "channels", "width_px", "height_px"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
