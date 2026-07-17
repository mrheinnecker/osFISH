



out="/g/schwab/marco/projects/osFISH/runs/regFISH02/treatment/R1_omezarr"
image="/g/schwab/marco/projects/osFISH/runs/regFISH02/treatment/R1_Merged.tif"


apptainer exec -B /home -B /scratch -B /g /g/schwab/marco/container_devel/eubi_v5.sif \
  

eubi to_zarr \
      $image \
      $out \
      --x_unit nm \
      --y_unit nm \
      --z_unit nm \
      --x_scale "10" \
      --y_scale "10" \
      --z_scale "50" \
      --concatenation_axes z \
      --z_tag "s" \
      --save_omexml True \
      --autochunk True \
      --zarr_format "2" 


