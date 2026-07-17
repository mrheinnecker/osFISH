#!/usr/bin/env bash
#SBATCH --job-name=regFISH02_omezarr
#SBATCH --output=/scratch/rheinnec/regFISH02_omezarr_%j.out
#SBATCH --error=/scratch/rheinnec/regFISH02_omezarr_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
set -euo pipefail

repo_dir="/g/schwab/marco/repos/osFISH"
run_dir="/g/schwab/marco/projects/osFISH/runs/regFISH02"
scratch_dir="/scratch/rheinnec/osFISH/regFISH02"
timestamp="$(date +%Y-%m-%d_%H-%M)"

mkdir -p "${run_dir}/logs" "${scratch_dir}"

bash "${repo_dir}/eubi_conv/osfish_omezarr_main.sh" cluster \
  --run_dir "${run_dir}" \
  --raw_lif_dir "${run_dir}/2026-06-22" \
  --outdir "${run_dir}/omezarr" \
  --extracted_tif_dir "${run_dir}/extracted_tifs" \
  --logdir "${run_dir}/logs/wfOMEZARR_${timestamp}" \
  --work_dir "${scratch_dir}/work_omezarr" \
  --default_z_scale_nm 1000 \
  --zarr_format 2 \
  --workflow_stage all \
  --s3_bucket "s3embl/temscreen/osFISH" \
  --sheet_mode google \
  --google_key "/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/1vFMQKq8MDs3nURapyu6odc58IlkWsWXj0sKeb06NLnE/edit?gid=0#gid=0" \
  --collection_table_sheet "ct" \
  --resume TRUE




