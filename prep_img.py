import os
import numpy as np
from aicspylibczi import CziFile
from czifile import CziFile as CzifileLegacy
import xml.etree.ElementTree as ET
from skimage.io import imsave
from glob import glob
import random
from PIL import Image, ImageDraw, ImageFont


# Set your image directory
image_dir = "/g/schwab/Marco/projects/osFISH/SF/13052025/raw"
output_dir = "/scratch/rheinnec/osFISH/image_analysis/test"
#os.makedirs(output_dir, exist_ok=True)


def add_scale_bar(image_rgb, pixel_size_um, bar_length_um=100, margin=50, bar_height=20):
    """
    Draws a white scale bar of `bar_length_um` in bottom-right of image_rgb.
    """
    h, w, _ = image_rgb.shape
    bar_length_px = int(bar_length_um / pixel_size_um)

    pil_img = Image.fromarray(image_rgb)
    draw = ImageDraw.Draw(pil_img)

    # Bar position and size
    #margin = 10
    #bar_height = 20
    x1 = w - bar_length_px - margin
    y1 = h - margin - bar_height
    x2 = w - margin
    y2 = h - margin

    # Draw bar
    draw.rectangle([x1, y1, x2, y2], fill=(255, 255, 255))

    # Draw label above bar
    font_size = 50
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
    except:
        font = ImageFont.load_default()

    label = f"{bar_length_um} μm"
    bbox = draw.textbbox((0, 0), label, font=font)
    text_width, text_height = bbox[2] - bbox[0], bbox[3] - bbox[1]

    text_x = x1 + (bar_length_px - text_width) // 2
    text_y = y1 - text_height - 20

    draw.text((text_x, text_y), label, fill=(255, 255, 255), font=font)
    return np.array(pil_img)


def get_pixel_size_um(czi_path):
    with CzifileLegacy(czi_path) as czi:
        metadata_xml = czi.metadata()

    root = ET.fromstring(metadata_xml)
    scaling = root.find(".//Scaling")
    x_scaling = scaling.find(".//Distance[@Id='X']")
    pixel_size_um = float(x_scaling.find("Value").text) * 1e6  # meters to microns

    return pixel_size_um



def apply_color_map(gray_img, hex_color):
    """
    Convert a grayscale image to an RGB image using a hex color code.
    """
    # Convert hex to RGB tuple (0–255)
    hex_color = hex_color.lstrip('#')
    rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))  # (R, G, B)

    # Normalize grayscale to 0–1 and scale each channel
    colored = np.stack([
        (gray_img / 255.0) * rgb[0],
        (gray_img / 255.0) * rgb[1],
        (gray_img / 255.0) * rgb[2],
    ], axis=-1)

    return np.clip(colored, 0, 255).astype(np.uint8)


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


def stitch_channels_horizontally(image_list):
    """
    Takes a list of RGB images (as numpy arrays) and stitches them side by side.
    Assumes all images are the same height.
    """
    return np.hstack(image_list)


# === LOAD FILES ===
czi_files = sorted(glob(os.path.join(image_dir, "*.czi")))

czi_files = czi_files[0:2]


# === EXTRACT CHANNEL NAMES ===
channel_names = get_channel_names(czi_files[0])
num_channels = len(channel_names)
print(f"Detected channels: {channel_names}")

# === PASS 1: Compute global max per channel ===
channel_max = np.zeros(num_channels)
print("Computing global max values...")

for fpath in czi_files:
    print("extract channel max:", fpath)
    czi = CziFile(fpath)
    img, shp = czi.read_image()
    img = np.squeeze(img)  # Expected shape: (C, Z, Y, X)

    for c in range(num_channels):
        channel_max[c] = max(channel_max[c], np.max(img[c]))

print(f"Max values per channel: {channel_max}")

# === PASS 2: Process each image and save random Z slice per channel ===

# Map channel names to desired hex colors
channel_colors = {
    'DAPI': '#0000FF',        # Blue
    'EGFP': '#00FF00',        # Green
    'Cy3': '#FF9900',         # Orange
    'Cy5': '#FF00FF',         # Magenta
    'At590': '#FF0000',       # Red
    'Bright': '#CCCCCC',      # Gray
}


rel_channels = ["DAPI", "Cy3", "Cy5"]

for fpath in czi_files:

    print("Processing images:", fpath)

    fname = os.path.basename(fpath).replace(".czi", "")
    img_output_dir = os.path.join(output_dir, fname)
    os.makedirs(img_output_dir, exist_ok=True)
    czi = CziFile(fpath)
    img, shp = czi.read_image()
    img = np.squeeze(img)  # shape: (C, Z, Y, X)

    z_dim = img.shape[1]
    #random_z = random.randint(0, z_dim - 1)

    for random_z in range(z_dim):
        print("z:", random_z)
        stitched_images = []  # collect for stitching

        #composite = np.zeros_like(img, dtype=np.float32)  

        for c in range(num_channels):
            slice_img = img[c, random_z]
            norm_img = (slice_img / channel_max[c]) * 255
            norm_img = np.clip(norm_img, 0, 255).astype(np.uint8)

            channel_name = channel_names[c]

            if channel_name in rel_channels:

                hex_color = channel_colors.get(channel_name, '#FFFFFF')  # fallback: white

                # Color + scale bar
                colored_img = apply_color_map(norm_img, hex_color)
                pixel_size_um = get_pixel_size_um(fpath)
                colored_img_with_bar = add_scale_bar(colored_img, pixel_size_um, bar_length_um=50)

                # Save individual image
                out_path = os.path.join(output_dir, f"{fname}_z{random_z}_{channel_name}.png")
                #imsave(out_path, colored_img_with_bar)

                # Add to stitched panel
                stitched_images.append(colored_img_with_bar)
            else:
                print("channel not selected")

        # === Save stitched panel ===
        panel = stitch_channels_horizontally(stitched_images)
        panel_path = os.path.join(img_output_dir, f"{fname}_z{random_z}_panel.png")
        imsave(panel_path, panel)



print("Done!")





