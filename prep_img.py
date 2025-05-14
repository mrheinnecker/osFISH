import os
import numpy as np
from aicspylibczi import CziFile
from czifile import CziFile as CzifileLegacy
import xml.etree.ElementTree as ET
from skimage.io import imsave
from glob import glob
import random


# Set your image directory
image_dir = "/scratch/rheinnec/osFISH/image_analysis"
output_dir = "/g/schwab/Marco/projects/osFISH/SF/13052025/prep"
#os.makedirs(output_dir, exist_ok=True)


# === HELPER: Extract channel names from metadata ===
def get_channel_names(czi_path):
    with CzifileLegacy(czi_path) as czi:
        metadata_xml = czi.metadata()

    root = ET.fromstring(metadata_xml)
    channels = root.findall(".//Channels/Channel")

    # Extract unique channel names, preserving order
    channel_names = [ch.attrib.get("Name", f"Channel_{i}") for i, ch in enumerate(channels)]
    channel_names = list(dict.fromkeys(channel_names))
    return channel_names

# === LOAD FILES ===
czi_files = sorted(glob(os.path.join(image_dir, "*.czi")))

# === EXTRACT CHANNEL NAMES ===
channel_names = get_channel_names(czi_files[0])
num_channels = len(channel_names)
print(f"Detected channels: {channel_names}")

# === PASS 1: Compute global max per channel ===
channel_max = np.zeros(num_channels)
print("Computing global max values...")

for fpath in czi_files:
    czi = CziFile(fpath)
    img, shp = czi.read_image()
    img = np.squeeze(img)  # Expected shape: (C, Z, Y, X)

    for c in range(num_channels):
        channel_max[c] = max(channel_max[c], np.max(img[c]))

print(f"Max values per channel: {channel_max}")

# === PASS 2: Process each image and save random Z slice per channel ===
print("Processing images...")
for fpath in czi_files:
    fname = os.path.basename(fpath).replace(".czi", "")
    czi = CziFile(fpath)
    img, shp = czi.read_image()
    img = np.squeeze(img)  # shape: (C, Z, Y, X)

    z_dim = img.shape[1]
    random_z = random.randint(0, z_dim - 1)

    for c in range(num_channels):
        slice_img = img[c, random_z]
        norm_img = (slice_img / channel_max[c]) * 255
        norm_img = np.clip(norm_img, 0, 255).astype(np.uint8)

        out_path = os.path.join(output_dir, f"{fname}_z{random_z}_{channel_names[c]}.png")
        imsave(out_path, norm_img)

print("Done!")





