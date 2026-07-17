#!/usr/bin/env python3
"""
Colorize a tau-contrast image with a colormap and modulate it by an intensity image
(same FOV) to suppress background.

Output is an RGB(A) PNG.

Examples
--------
# tau in channel 2 (axis 0), intensity in channel 0 (axis 0)
python tau_with_intensity_mask.py \
  --tau tau.tif --tau-channel-axis 0 --tau-channel 2 \
  --int intensity.tif --int-channel-axis 0 --int-channel 0 \
  --out tau_masked.png --mode alpha --colorbar

# zstack: (Z,C,Y,X) pick Z=5 by choosing --z 5 (optional), else uses z=0
python tau_with_intensity_mask.py \
  --tau flim.tif --tau-channel-axis 1 --tau-channel 1 --tau-z 5 \
  --int flim.tif --int-channel-axis 1 --int-channel 0 --int-z 5 \
  --out out.png
"""

import argparse
from pathlib import Path

import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt
import matplotlib.cm as cm
def extract_2d(arr: np.ndarray, channel_axis: int, channel: int, z_axis: int | None, z: int | None) -> np.ndarray:
    """
    Extract a 2D plane from an N-D array.
    If arr is already 2D, return it unchanged.
    Otherwise: select given channel (and optional z), and fix all other axes to 0.
    """
    # ✅ If it's already a 2D image, nothing to select
    if arr.ndim == 2:
        return arr

    if channel_axis < 0 or channel_axis >= arr.ndim:
        raise ValueError(f"channel_axis {channel_axis} out of bounds for ndim={arr.ndim}")
    if channel < 0 or channel >= arr.shape[channel_axis]:
        raise ValueError(f"channel {channel} out of bounds for axis size {arr.shape[channel_axis]}")

    if (z_axis is None) ^ (z is None):
        raise ValueError("Provide both --*-z-axis and --*-z, or neither.")

    if z_axis is not None:
        if z_axis < 0 or z_axis >= arr.ndim:
            raise ValueError(f"z_axis {z_axis} out of bounds for ndim={arr.ndim}")
        if z < 0 or z >= arr.shape[z_axis]:
            raise ValueError(f"z {z} out of bounds for axis size {arr.shape[z_axis]}")

    slc = []
    for ax in range(arr.ndim):
        if ax == channel_axis:
            slc.append(channel)
        elif z_axis is not None and ax == z_axis:
            slc.append(z)
        else:
            slc.append(0)

    plane = arr[tuple(slc)]

    # Sometimes slicing can still leave singleton dimensions, squeeze them out
    plane = np.squeeze(plane)

    if plane.ndim != 2:
        raise ValueError(f"Extracted plane has shape {plane.shape}, expected 2D (Y,X).")

    return plane


def robust_normalize(img: np.ndarray, pmin: float, pmax: float) -> tuple[np.ndarray, float, float]:
    img = img.astype(np.float32)
    finite = np.isfinite(img)
    if not np.any(finite):
        raise ValueError("Image has no finite values.")
    vmin = np.percentile(img[finite], pmin)
    vmax = np.percentile(img[finite], pmax)
    if vmax <= vmin:
        vmin = float(np.min(img[finite]))
        vmax = float(np.max(img[finite]))
        if vmax <= vmin:
            return np.zeros_like(img, dtype=np.float32), float(vmin), float(vmax)
    norm = (np.clip(img, vmin, vmax) - vmin) / (vmax - vmin)
    return norm, float(vmin), float(vmax)


def gamma_correct(x: np.ndarray, gamma: float) -> np.ndarray:
    if gamma == 1.0:
        return x
    x = np.clip(x, 0, 1)
    return np.power(x, gamma)


def main():
    ap = argparse.ArgumentParser(description="Colorize tau with intensity masking (same FOV).")

    # tau inputs
    ap.add_argument("--tau", required=True, help="Tau TIFF path")
    ap.add_argument("--tau-channel-axis", type=int, required=True, help="Tau channel axis (0-based)")
    ap.add_argument("--tau-channel", type=int, required=True, help="Tau channel index (0-based)")
    ap.add_argument("--tau-z-axis", type=int, default=None, help="Optional tau Z axis (0-based)")
    ap.add_argument("--tau-z", type=int, default=None, help="Optional tau Z index (0-based)")

    # intensity inputs
    ap.add_argument("--int", dest="intp", required=True, help="Intensity TIFF path")
    ap.add_argument("--int-channel-axis", type=int, required=True, help="Intensity channel axis (0-based)")
    ap.add_argument("--int-channel", type=int, required=True, help="Intensity channel index (0-based)")
    ap.add_argument("--int-z-axis", type=int, default=None, help="Optional intensity Z axis (0-based)")
    ap.add_argument("--int-z", type=int, default=None, help="Optional intensity Z index (0-based)")

    # output / viz
    ap.add_argument("--out", required=True, help="Output PNG path")
    ap.add_argument("--cmap", default="turbo", help="Colormap for tau (default: turbo)")
    ap.add_argument("--tau-pmin", type=float, default=1.0, help="Tau lower percentile clip")
    ap.add_argument("--tau-pmax", type=float, default=99.0, help="Tau upper percentile clip")
    ap.add_argument("--int-pmin", type=float, default=1.0, help="Intensity lower percentile clip")
    ap.add_argument("--int-pmax", type=float, default=99.5, help="Intensity upper percentile clip")
    ap.add_argument("--int-gamma", type=float, default=0.7, help="Gamma for intensity mask (default: 0.7)")
    ap.add_argument("--mode", choices=["alpha", "multiply"], default="alpha",
                    help="alpha: intensity as transparency; multiply: intensity scales RGB (default: alpha)")
    ap.add_argument("--alpha-min", type=float, default=0.0,
                    help="Minimum alpha/brightness for background (0..1). (default: 0)")
    ap.add_argument("--alpha-max", type=float, default=1.0,
                    help="Maximum alpha/brightness (0..1). (default: 1)")
    ap.add_argument("--colorbar", action="store_true", help="Add tau colorbar (separate render, still saved into same PNG)")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    # Load
    tau_arr = tiff.imread(args.tau)
    int_arr = tiff.imread(args.intp)
    print(f"Tau TIFF shape: {tau_arr.shape}")
    print(f"Intensity TIFF shape: {int_arr.shape}")

    tau_plane = extract_2d(tau_arr, args.tau_channel_axis, args.tau_channel, args.tau_z_axis, args.tau_z)
    int_plane = extract_2d(int_arr, args.int_channel_axis, args.int_channel, args.int_z_axis, args.int_z)

    if tau_plane.shape != int_plane.shape:
        raise ValueError(f"Tau plane shape {tau_plane.shape} != intensity plane shape {int_plane.shape} (must match).")

    # Normalize tau for colormap
    tau_norm, tau_vmin, tau_vmax = robust_normalize(tau_plane, args.tau_pmin, args.tau_pmax)

    # Normalize intensity for masking
    int_norm, int_vmin, int_vmax = robust_normalize(int_plane, args.int_pmin, args.int_pmax)
    int_mask = gamma_correct(int_norm, args.int_gamma)
    int_mask = np.clip(args.alpha_min + int_mask * (args.alpha_max - args.alpha_min), 0, 1)

    # Map tau to RGB via colormap
    cmap = cm.get_cmap(args.cmap)
    rgba = cmap(tau_norm)  # (H,W,4) floats 0..1

    if args.mode == "alpha":
        rgba[..., 3] = rgba[..., 3] * int_mask  # modulate alpha
    else:  # multiply
        rgba[..., :3] = rgba[..., :3] * int_mask[..., None]
        rgba[..., 3] = 1.0

    # Save
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not args.colorbar:
        plt.imsave(str(out_path), rgba)
    else:
        # Render with a colorbar in a single figure
        fig = plt.figure()
        ax = plt.gca()
        ax.imshow(rgba)
        ax.axis("off")

        # Build a scalar mappable so colorbar corresponds to tau in original units
        sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=tau_vmin, vmax=tau_vmax))
        cbar = plt.colorbar(sm, fraction=0.046, pad=0.04)
        cbar.set_label("τ (original units)")

        plt.tight_layout()
        plt.savefig(str(out_path), dpi=args.dpi, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)

    print(f"Saved: {out_path}")
    print(
        f"Tau clip: {args.tau_pmin}-{args.tau_pmax} (vmin={tau_vmin:.4g}, vmax={tau_vmax:.4g}); "
        f"Intensity clip: {args.int_pmin}-{args.int_pmax} (vmin={int_vmin:.4g}, vmax={int_vmax:.4g}); "
        f"mode={args.mode}, gamma={args.int_gamma}"
    )


if __name__ == "__main__":
    main()
