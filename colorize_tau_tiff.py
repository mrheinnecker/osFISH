#!/usr/bin/env python3
"""
Colorize a grayscale tau-contrast image stored in a TIFF (e.g., a channel in a multi-channel TIFF).

Examples
--------
# simplest (single-plane TIFF):
python colorize_tau_tiff.py --in tau.tif --out tau_color.png

# multi-page / multi-channel TIFF: choose which page to colorize:
python colorize_tau_tiff.py --in multi.tif --page 3 --out tau_color.png

# percentile clipping (robust contrast) + add colorbar:
python colorize_tau_tiff.py --in multi.tif --page 3 --pmin 1 --pmax 99 --colorbar
"""

import argparse
from pathlib import Path

import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt


def pick_2d_plane(arr: np.ndarray, page: int | None) -> np.ndarray:
    """
    Try to get a 2D plane from a TIFF array.
    - If arr is already 2D -> return it
    - If arr is 3D/4D -> pick 'page' from the first axis by default
    """
    if arr.ndim == 2:
        return arr

    if page is None:
        page = 0

    if arr.ndim >= 3:
        # Common cases:
        # - (pages, y, x)
        # - (channels, y, x)
        # - (z, y, x)
        # We'll treat axis 0 as the stack/channel axis.
        if page < 0 or page >= arr.shape[0]:
            raise ValueError(f"--page {page} out of range for axis0 with size {arr.shape[0]}")
        plane = arr[page]
        if plane.ndim != 2:
            # Sometimes you might get (y, x, c) etc.
            raise ValueError(
                f"Selected plane has shape {plane.shape} (ndim={plane.ndim}), not 2D. "
                "If your TIFF is (y,x,c), you’ll need to split channels differently."
            )
        return plane

    raise ValueError(f"Unsupported TIFF array ndim={arr.ndim} with shape={arr.shape}")


def robust_normalize(img: np.ndarray, pmin: float, pmax: float) -> tuple[np.ndarray, float, float]:
    """
    Normalize image to [0,1] using percentile clipping (robust).
    Returns (norm_img, vmin, vmax) where vmin/vmax are in original units.
    """
    img = img.astype(np.float32)

    finite = np.isfinite(img)
    if not np.any(finite):
        raise ValueError("Image has no finite values.")

    vmin = np.percentile(img[finite], pmin)
    vmax = np.percentile(img[finite], pmax)
    if vmax <= vmin:
        # Fallback to min/max
        vmin = float(np.min(img[finite]))
        vmax = float(np.max(img[finite]))
        if vmax <= vmin:
            # Completely flat image
            norm = np.zeros_like(img, dtype=np.float32)
            return norm, vmin, vmax

    clipped = np.clip(img, vmin, vmax)
    norm = (clipped - vmin) / (vmax - vmin)
    return norm, float(vmin), float(vmax)


def main():
    ap = argparse.ArgumentParser(description="Colorize tau-contrast grayscale TIFF into a multicolor PNG.")
    ap.add_argument("--in", dest="inp", required=True, help="Input TIFF path")
    ap.add_argument("--out", dest="out", required=True, help="Output image path (e.g., .png)")
    ap.add_argument("--page", type=int, default=None,
                    help="Which plane/page/channel to use if TIFF is a stack. (default: 0)")
    ap.add_argument("--cmap", default="turbo",
                    help="Matplotlib colormap name (default: turbo). Good options: turbo, viridis, plasma, magma, inferno")
    ap.add_argument("--pmin", type=float, default=1.0, help="Lower percentile for clipping (default: 1)")
    ap.add_argument("--pmax", type=float, default=99.0, help="Upper percentile for clipping (default: 99)")
    ap.add_argument("--colorbar", action="store_true", help="Include a colorbar in the output")
    ap.add_argument("--title", default=None, help="Optional title for the figure")
    ap.add_argument("--dpi", type=int, default=200, help="Output DPI (default: 200)")
    args = ap.parse_args()

    in_path = Path(args.inp)
    out_path = Path(args.out)

    # Read TIFF (this handles many TIFF layouts)
    arr = tiff.imread(str(in_path))

    plane = pick_2d_plane(arr, args.page)

    norm, vmin, vmax = robust_normalize(plane, args.pmin, args.pmax)

    # Render
    plt.figure()
    im = plt.imshow(norm, cmap=args.cmap, vmin=0.0, vmax=1.0)
    plt.axis("off")

    if args.title:
        plt.title(args.title)

    if args.colorbar:
        # show colorbar with original-value scaling
        cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
        cbar.set_label("τ (a.u. / original units)")
        # relabel ticks from [0..1] to [vmin..vmax]
        ticks = cbar.get_ticks()
        ticklabels = [f"{vmin + t*(vmax-vmin):.3g}" for t in ticks]
        cbar.set_ticklabels(ticklabels)

    plt.tight_layout()
    plt.savefig(str(out_path), dpi=args.dpi, bbox_inches="tight", pad_inches=0.02)
    plt.close()

    print(f"Saved colorized image to: {out_path}")
    print(f"Used plane/page: {args.page if args.page is not None else 0} | "
          f"clip percentiles: {args.pmin}-{args.pmax} | vmin={vmin:.4g}, vmax={vmax:.4g} | cmap={args.cmap}")


if __name__ == "__main__":
    main()
