# osFISH single-LIF to OME-Zarr workflow

This is a downsized copy of the PLASTIC LIF conversion path for one Leica `.lif`
file at a time. It does not read an input table.

Pipeline:

```text
single .lif
  -> extract LIF metadata and channel names
  -> export each LIF scene to an intermediate OME-TIFF
  -> eubi to_zarr for each scene
  -> patch OME-Zarr omero/channel metadata for each scene
  -> optional S3 upload
  -> MoBIE collection table with LIF scenes arranged on a grid
```

## Run

Convert locally on an interactive node without upload:

```bash
bash lif_to_omezarr/osfish_lif_main.sh interactive \
  --input_lif /path/to/image.lif \
  --main_dir /scratch/rheinnec/osFISH/lif_to_omezarr \
  --upload FALSE
```

Convert and upload:

```bash
bash lif_to_omezarr/osfish_lif_main.sh cluster \
  --input_lif /path/to/image.lif \
  --dataset_name SF02_example \
  --main_dir /scratch/rheinnec/osFISH/lif_to_omezarr \
  --s3_bucket s3embl/imatrec/central_data_processing/osfish \
  --collection_table_url "https://docs.google.com/spreadsheets/d/.../edit?gid=0#gid=0" \
  --collection_table_sheet osfish_collection_table \
  --upload TRUE
```

If LIF metadata is missing physical scales, pass fallbacks:

```bash
--default_x_scale 250 --default_y_scale 250 --default_z_scale 1000 --scale_unit nm
```

## Outputs

The converted data are published to:

```text
<main_dir>/processed/<dataset_name>_<scene_index>_<scene_name>.ome.zarr
```

The MoBIE collection table is written to the workflow log directory as:

```text
mobie_collection_table.tsv
```

Channels are detected from LIF metadata through `BioImage`/`bioio-lif`, with
optional `.lifext` sidecar metadata used when present.
