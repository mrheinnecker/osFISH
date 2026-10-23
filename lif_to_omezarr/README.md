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
  -> MoBIE collection table with local file paths, scenes arranged on a grid
```

## Run

Convert one file on an interactive node:

```bash
bash lif_to_omezarr/osfish_lif_main.sh interactive \
  --input_lif /path/to/image.lif \
  --main_dir /scratch/rheinnec/osFISH/lif_to_omezarr \
  --outdir /g/schwab/marco/projects/osFISH/omezarr \
  --uri_from /g/schwab --uri_to 'A:\'
```

`--uri_from`/`--uri_to` rewrite the collection-table paths from the cluster
path to the path the MoBIE machine sees (here the Windows drive `A:\` mapped
to `/g/schwab`). Without them the table holds the cluster paths. `--outdir`
must be below `--uri_from`, so it should be on `/g`, not `/scratch`.

## Batch: all runs into one collection table

Build an inventory of every image file below the runs directory, linked to the
run's sheet in `samples.xlsx` (sheet name = run folder name):

```bash
python3 lif_to_omezarr/build_run_inventory.py \
  --runs-dir /g/schwab/marco/projects/osFISH/runs \
  --samples-xlsx /g/schwab/marco/projects/osFISH/samples.xlsx \
  --output /g/schwab/marco/projects/osFISH/lif_inventory.tsv
```

Which files are included:

- `.lif` (Leica) and `.czi` (Zeiss) everywhere.
- `.ome.tif` (LAS X project exports) only in runs without LIFs; plain TIFFs and
  OME-TIFFs next to LIFs are exports of LIF scenes and are skipped.
- Raw tiles are skipped when a `*_Merged` version exists, both for OME-TIFF
  files and for scenes inside a LIF (`--keep-unmerged-tiles` in
  `export_lif_scenes.py` turns this off).
- A file copied into another run (same folder, name and size) is kept once.
- A folder without a sheet whose files are all named `<sheet>_*` uses that
  sheet (e.g. `13052025` → `SF02`).

Samples are matched from the file name to the sheet's sample columns (also
Zeiss numbering like `SF02_03_01` → `s3_...`). Check
the TSV and fill in `sample`/`condition` by hand where the match is empty or
wrong; delete rows for LIFs that should not be converted. `--runs SF09,SF10`
limits the inventory to some runs.

Then convert everything in the inventory:

```bash
bash lif_to_omezarr/osfish_lif_main.sh cluster \
  --input_table /g/schwab/marco/projects/osFISH/lif_inventory.tsv \
  --main_dir /scratch/rheinnec/osFISH/lif_to_omezarr \
  --outdir /g/schwab/marco/projects/osFISH/omezarr \
  --uri_from /g/schwab --uri_to 'A:\' \
  --default_z_scale 1000
```

Add `--collection_table_url ...` to also write the table to a Google Sheet.

The collection table has one MoBIE view and grid per run (pick the run from
the view dropdown) and annotation columns `run`, `run_date`, `source_file`,
`format`, `scene`, `sample` and `condition`. A file whose metadata cannot be
read, or that still fails conversion after two retries with more memory/time,
is skipped and listed in the Nextflow report.

If metadata is missing physical scales, pass fallbacks. They only fill values
missing from the file (and per-scene pixel sizes always win):

```bash
--default_x_scale 250 --default_y_scale 250 --default_z_scale 1000 --scale_unit nm
```

## Outputs

The converted data are published to:

```text
<outdir>/<dataset_name>_<scene_index>_<scene_name>.ome.zarr
```

`<outdir>` defaults to `<main_dir>/processed`. The MoBIE collection table is
written next to the data and to the workflow log directory as:

```text
<outdir>/mobie_collection_table.tsv
```

Channels are detected from LIF metadata through `BioImage`/`bioio-lif`, with
optional `.lifext` sidecar metadata used when present.
