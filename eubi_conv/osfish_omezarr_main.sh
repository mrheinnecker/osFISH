#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp="$(date +%Y-%m-%d_%H-%M)"

mode="${MODE:-interactive}"
resume="${RESUME:-TRUE}"
run_dir="${RUN_DIR:-/g/schwab/marco/projects/osFISH/runs/regFISH02}"
raw_lif_dir="${RAW_LIF_DIR:-${run_dir}/2026-06-22}"
outdir="${OUTDIR:-${run_dir}/omezarr}"
extracted_tif_dir="${EXTRACTED_TIF_DIR:-${run_dir}/extracted_tifs}"
logdir="${LOGDIR:-${run_dir}/logs/wfOMEZARR_${timestamp}}"
work_dir="${WORK_DIR:-${run_dir}/work_omezarr}"
default_z_scale_nm="${DEFAULT_Z_SCALE_NM:-1000}"
zarr_format="${ZARR_FORMAT:-2}"
workflow_stage="${WORKFLOW_STAGE:-all}"
s3_bucket="${S3_BUCKET:-s3embl/temscreen/osFISH}"
sheet_mode="${SHEET_MODE:-google}"
google_key="${GOOGLE_KEY:-/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json}"
collection_table_url="${COLLECTION_TABLE_URL:-https://docs.google.com/spreadsheets/d/1vFMQKq8MDs3nURapyu6odc58IlkWsWXj0sKeb06NLnE/edit?gid=0#gid=0}"
collection_table_sheet="${COLLECTION_TABLE_SHEET:-ct}"
metadata_root="${METADATA_ROOT:-$(dirname "$run_dir")}"
eubi_extra_args="${EUBI_EXTRA_ARGS:-}"

usage() {
  cat <<EOF
Usage:
  bash eubi_conv/osfish_omezarr_main.sh [local|interactive|cluster] [options]

Options:
  --run_dir PATH              Run directory. Default: /g/schwab/marco/projects/osFISH/runs/regFISH02
  --raw_lif_dir PATH          Directory with source .lif files. Default: <run_dir>/2026-06-22
  --outdir PATH               Output directory for OME-Zarrs. Default: <run_dir>/omezarr
  --extracted_tif_dir PATH    Output directory for extracted TIFFs. Default: <run_dir>/extracted_tifs
  --logdir PATH               Log/table directory. Default: <run_dir>/logs/wfOMEZARR_<timestamp>
  --work_dir PATH             Nextflow work directory. Default: <run_dir>/work_omezarr
  --default_z_scale_nm VALUE  Fallback Z scale for single-plane images. Default: 1000
  --zarr_format VALUE         OME-Zarr format. Default: 2
  --workflow_stage VALUE      process, all, or collection. Default: all
  --s3_bucket VALUE           Upload bucket/prefix. Default: s3embl/temscreen/osFISH
  --sheet_mode VALUE          local or google. Default: google
  --google_key PATH           Google service-account JSON key.
  --collection_table_url URL  Google Sheet URL for MoBIE collection table.
  --collection_table_sheet    Google Sheet tab. Default: ct
  --metadata_root PATH        Root containing osFISH run dirs/logs to combine. Default: dirname(<run_dir>)
  --eubi_extra_args VALUE     Extra arguments appended to eubi to_zarr.
  --resume TRUE|FALSE         Add Nextflow -resume. Default: TRUE
  --help                      Show this message.
EOF
}

to_upper_bool() {
  case "$1" in
    TRUE|true|1|yes|YES) echo "TRUE" ;;
    FALSE|false|0|no|NO) echo "FALSE" ;;
    *) echo "Boolean value must be TRUE or FALSE, got: $1" >&2; exit 1 ;;
  esac
}

if [[ $# -gt 0 && "${1:-}" != --* ]]; then
  mode="$1"
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h) usage; exit 0 ;;
    --run_dir|--run-dir) run_dir="${2:?--run_dir requires a path}"; shift 2 ;;
    --run_dir=*|--run-dir=*) run_dir="${1#*=}"; shift ;;
    --raw_lif_dir|--raw-lif-dir) raw_lif_dir="${2:?--raw_lif_dir requires a path}"; shift 2 ;;
    --raw_lif_dir=*|--raw-lif-dir=*) raw_lif_dir="${1#*=}"; shift ;;
    --outdir) outdir="${2:?--outdir requires a path}"; shift 2 ;;
    --outdir=*) outdir="${1#*=}"; shift ;;
    --extracted_tif_dir|--extracted-tif-dir) extracted_tif_dir="${2:?--extracted_tif_dir requires a path}"; shift 2 ;;
    --extracted_tif_dir=*|--extracted-tif-dir=*) extracted_tif_dir="${1#*=}"; shift ;;
    --logdir) logdir="${2:?--logdir requires a path}"; shift 2 ;;
    --logdir=*) logdir="${1#*=}"; shift ;;
    --work_dir|--work-dir) work_dir="${2:?--work_dir requires a path}"; shift 2 ;;
    --work_dir=*|--work-dir=*) work_dir="${1#*=}"; shift ;;
    --default_z_scale_nm|--default-z-scale-nm) default_z_scale_nm="${2:?--default_z_scale_nm requires a value}"; shift 2 ;;
    --default_z_scale_nm=*|--default-z-scale-nm=*) default_z_scale_nm="${1#*=}"; shift ;;
    --zarr_format|--zarr-format) zarr_format="${2:?--zarr_format requires a value}"; shift 2 ;;
    --zarr_format=*|--zarr-format=*) zarr_format="${1#*=}"; shift ;;
    --workflow_stage|--workflow-stage) workflow_stage="${2:?--workflow_stage requires a value}"; shift 2 ;;
    --workflow_stage=*|--workflow-stage=*) workflow_stage="${1#*=}"; shift ;;
    --s3_bucket|--s3-bucket) s3_bucket="${2:?--s3_bucket requires a value}"; shift 2 ;;
    --s3_bucket=*|--s3-bucket=*) s3_bucket="${1#*=}"; shift ;;
    --sheet_mode|--sheet-mode) sheet_mode="${2:?--sheet_mode requires a value}"; shift 2 ;;
    --sheet_mode=*|--sheet-mode=*) sheet_mode="${1#*=}"; shift ;;
    --google_key|--google-key) google_key="${2:?--google_key requires a path}"; shift 2 ;;
    --google_key=*|--google-key=*) google_key="${1#*=}"; shift ;;
    --collection_table_url|--collection-table-url) collection_table_url="${2:?--collection_table_url requires a URL}"; shift 2 ;;
    --collection_table_url=*|--collection-table-url=*) collection_table_url="${1#*=}"; shift ;;
    --collection_table_sheet|--collection-table-sheet) collection_table_sheet="${2:?--collection_table_sheet requires a value}"; shift 2 ;;
    --collection_table_sheet=*|--collection-table-sheet=*) collection_table_sheet="${1#*=}"; shift ;;
    --metadata_root|--metadata-root) metadata_root="${2:?--metadata_root requires a path}"; shift 2 ;;
    --metadata_root=*|--metadata-root=*) metadata_root="${1#*=}"; shift ;;
    --eubi_extra_args|--eubi-extra-args) eubi_extra_args="${2:?--eubi_extra_args requires a value}"; shift 2 ;;
    --eubi_extra_args=*|--eubi-extra-args=*) eubi_extra_args="${1#*=}"; shift ;;
    --resume) resume="$(to_upper_bool "${2:?--resume requires TRUE or FALSE}")"; shift 2 ;;
    --resume=*) resume="$(to_upper_bool "${1#*=}")"; shift ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

case "$mode" in
  local|interactive|cluster) ;;
  *) echo "Unknown profile: $mode" >&2; usage >&2; exit 1 ;;
esac

if command -v module >/dev/null 2>&1; then
  module load Nextflow/24.10.4 || true
fi

mkdir -p "$outdir" "$extracted_tif_dir" "$logdir" "$work_dir"

nextflow_args=(
  run "${script_dir}/wfOSFISH_OMEZARR.nf"
  -c "${script_dir}/nextflow.config"
  -profile "$mode"
  -work-dir "$work_dir"
  --script_dir "$script_dir"
  --run_dir "$run_dir"
  --raw_lif_dir "$raw_lif_dir"
  --outdir "$outdir"
  --extracted_tif_dir "$extracted_tif_dir"
  --logdir "$logdir"
  --default_z_scale_nm "$default_z_scale_nm"
  --zarr_format "$zarr_format"
  --workflow_stage "$workflow_stage"
  --s3_bucket "$s3_bucket"
  --sheet_mode "$sheet_mode"
  --google_key "$google_key"
  --collection_table_url "$collection_table_url"
  --collection_table_sheet "$collection_table_sheet"
  --metadata_root "$metadata_root"
)

if [[ -n "$eubi_extra_args" ]]; then
  nextflow_args+=(--eubi_extra_args "$eubi_extra_args")
fi

if [[ "$resume" == "TRUE" ]]; then
  nextflow_args+=(-resume)
fi

nextflow "${nextflow_args[@]}"
