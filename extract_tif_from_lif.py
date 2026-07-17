## execute singularity shell -B /scratch -B /g /g/schwab/marco/container_legacy/py_liffile.sif

from pathlib import Path
import numpy as np
import tifffile
from readlif.reader import LifFile
import re

lif_path = Path("/g/schwab/marco/projects/tax")
out_dir = Path("/g/schwab/marco/projects/osFISH/runs/regFISH02/control")
out_dir.mkdir(exist_ok=True)

lif = LifFile(str(lif_path))




for i, img in enumerate(lif.get_iter_image()):
    print(i, img.name, img.dims)
    safe_name = re.sub(r'[\\/*?:"<>|]', "_", img.name).strip()
    dims = img.dims
    n_channels = img.channels
    channels = []
    for c in range(n_channels):
        print(c)
        frame = img.get_frame(c=c)
        arr = np.array(frame)
        channels.append(arr)
    out_arr = np.stack(channels, axis=0)
    out_path = out_dir / f"{safe_name}.tif"
    tifffile.imwrite(
        out_path,
        out_arr,
        imagej=True,
        metadata={"axes": "CYX"}
    )


