#!/usr/bin/env bash
set -euo pipefail

PATH="/usr/local/bin:/usr/bin:/bin:${PATH:-}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
timestamp="$(date +%Y-%m-%d_%H-%M)"

mode="${MODE:-interactive}"
input_lif="${INPUT_LIF:-}"
dataset_name="${DATASET_NAME:-}"
main_dir="${OSFISH_OMEZARR_DIR:-/scratch/rheinnec/osFISH/lif_to_omezarr}"
outdir="${OUTDIR:-}"
logdir="${LOGDIR:-}"
work_dir="${WORK_DIR:-}"
s3_bucket="${S3_BUCKET:-s3embl/imatrec/central_data_processing/osfish}"
collection_table_url="${COLLECTION_TABLE_URL:-}"
collection_table_sheet="${COLLECTION_TABLE_SHEET:-osfish_collection_table}"
google_key="${GOOGLE_KEY:-${script_dir}/trec-tem-screen-e98a2e03f58b.json}"
upload="${UPLOAD:-FALSE}"
zarr_format="${ZARR_FORMAT:-2}"
eubi_extra_args="${EUBI_EXTRA_ARGS:-}"
default_x_scale="${DEFAULT_X_SCALE:-}"
default_y_scale="${DEFAULT_Y_SCALE:-}"
default_z_scale="${DEFAULT_Z_SCALE:-}"
scale_unit="${SCALE_UNIT:-nm}"
resume="${RESUME:-TRUE}"

usage() {
  cat <<EOF
Usage:
  bash osfish_lif_main.sh [local|interactive|cluster] --input_lif /path/to/image.lif [options]

Options:
  --input_lif PATH              Single Leica LIF file to convert.
  --dataset_name VALUE          Optional output dataset name. Defaults to sanitized LIF basename.
  --main_dir PATH               Base run directory. Default: /scratch/rheinnec/osFISH/lif_to_omezarr.
  --outdir PATH                 Converted-output directory.
  --logdir PATH                 Workflow log directory.
  --work_dir PATH               Nextflow work directory.
  --s3_bucket VALUE             S3 bucket/prefix for upload and collection-table URI.
  --collection_table_url URL    Optional Google Sheet URL to write the MoBIE table.
  --collection_table_sheet NAME Google Sheet tab name. Default: osfish_collection_table.
  --google_key PATH             Google service-account JSON key.
  --upload TRUE|FALSE           Upload OME-Zarr to S3. Default: FALSE.
  --zarr_format VALUE           OME-Zarr format version. Default: 2.
  --eubi_extra_args VALUE       Extra raw arguments appended to eubi to_zarr.
  --default_x_scale VALUE       Fallback X pixel size if missing from LIF metadata.
  --default_y_scale VALUE       Fallback Y pixel size if missing from LIF metadata.
  --default_z_scale VALUE       Fallback Z spacing if missing from LIF metadata.
  --scale_unit VALUE            Unit for fallback scales: nm, um, or mm. Default: nm.
  --resume TRUE|FALSE           Add Nextflow -resume when TRUE. Default: TRUE.
  --help                        Show this message.
EOF
}

to_upper_bool() {
  case "$1" in
    TRUE|true|1|yes|YES) echo "TRUE" ;;
    FALSE|false|0|no|NO) echo "FALSE" ;;
    *)
      echo "Boolean value must be TRUE or FALSE, got: $1" >&2
      exit 1
      ;;
  esac
}

if [[ $# -gt 0 && "${1:-}" != --* && "${1:-}" != "-resume" ]]; then
  mode="$1"
  shift
fi

case "$mode" in
  local)
    profile="local"
    main_dir="${OSFISH_OMEZARR_DIR:-${main_dir}}"
    ;;
  interactive)
    profile="interactive"
    if command -v module >/dev/null 2>&1; then
      module load Nextflow/24.10.4
    fi
    ;;
  cluster)
    profile="cluster"
    if command -v module >/dev/null 2>&1; then
      module load Nextflow/24.10.4
    fi
    ;;
  *)
    echo "Unknown profile: $mode" >&2
    usage
    exit 1
    ;;
esac

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --input_lif|--input-lif)
      input_lif="${2:?--input_lif requires a path}"
      shift 2
      ;;
    --input_lif=*|--input-lif=*)
      input_lif="${1#*=}"
      shift
      ;;
    --dataset_name|--dataset-name)
      dataset_name="${2:?--dataset_name requires a value}"
      shift 2
      ;;
    --dataset_name=*|--dataset-name=*)
      dataset_name="${1#*=}"
      shift
      ;;
    --main_dir|--main-dir)
      main_dir="${2:?--main_dir requires a path}"
      shift 2
      ;;
    --main_dir=*|--main-dir=*)
      main_dir="${1#*=}"
      shift
      ;;
    --outdir)
      outdir="${2:?--outdir requires a path}"
      shift 2
      ;;
    --outdir=*)
      outdir="${1#*=}"
      shift
      ;;
    --logdir)
      logdir="${2:?--logdir requires a path}"
      shift 2
      ;;
    --logdir=*)
      logdir="${1#*=}"
      shift
      ;;
    --work_dir|--work-dir)
      work_dir="${2:?--work_dir requires a path}"
      shift 2
      ;;
    --work_dir=*|--work-dir=*)
      work_dir="${1#*=}"
      shift
      ;;
    --s3_bucket|--s3-bucket)
      s3_bucket="${2:?--s3_bucket requires a value}"
      shift 2
      ;;
    --s3_bucket=*|--s3-bucket=*)
      s3_bucket="${1#*=}"
      shift
      ;;
    --collection_table_url|--collection-table-url)
      collection_table_url="${2:?--collection_table_url requires a URL}"
      shift 2
      ;;
    --collection_table_url=*|--collection-table-url=*)
      collection_table_url="${1#*=}"
      shift
      ;;
    --collection_table_sheet|--collection-table-sheet)
      collection_table_sheet="${2:?--collection_table_sheet requires a value}"
      shift 2
      ;;
    --collection_table_sheet=*|--collection-table-sheet=*)
      collection_table_sheet="${1#*=}"
      shift
      ;;
    --google_key|--google-key)
      google_key="${2:?--google_key requires a path}"
      shift 2
      ;;
    --google_key=*|--google-key=*)
      google_key="${1#*=}"
      shift
      ;;
    --upload)
      upload="$(to_upper_bool "${2:?--upload requires TRUE or FALSE}")"
      shift 2
      ;;
    --upload=*)
      upload="$(to_upper_bool "${1#*=}")"
      shift
      ;;
    --zarr_format|--zarr-format)
      zarr_format="${2:?--zarr_format requires a value}"
      shift 2
      ;;
    --zarr_format=*|--zarr-format=*)
      zarr_format="${1#*=}"
      shift
      ;;
    --eubi_extra_args|--eubi-extra-args)
      eubi_extra_args="${2:?--eubi_extra_args requires a value}"
      shift 2
      ;;
    --eubi_extra_args=*|--eubi-extra-args=*)
      eubi_extra_args="${1#*=}"
      shift
      ;;
    --default_x_scale|--default-x-scale)
      default_x_scale="${2:?--default_x_scale requires a value}"
      shift 2
      ;;
    --default_x_scale=*|--default-x-scale=*)
      default_x_scale="${1#*=}"
      shift
      ;;
    --default_y_scale|--default-y-scale)
      default_y_scale="${2:?--default_y_scale requires a value}"
      shift 2
      ;;
    --default_y_scale=*|--default-y-scale=*)
      default_y_scale="${1#*=}"
      shift
      ;;
    --default_z_scale|--default-z-scale)
      default_z_scale="${2:?--default_z_scale requires a value}"
      shift 2
      ;;
    --default_z_scale=*|--default-z-scale=*)
      default_z_scale="${1#*=}"
      shift
      ;;
    --scale_unit|--scale-unit)
      scale_unit="${2:?--scale_unit requires a value}"
      shift 2
      ;;
    --scale_unit=*|--scale-unit=*)
      scale_unit="${1#*=}"
      shift
      ;;
    --resume)
      resume="$(to_upper_bool "${2:?--resume requires TRUE or FALSE}")"
      shift 2
      ;;
    --resume=*)
      resume="$(to_upper_bool "${1#*=}")"
      shift
      ;;
    -resume)
      resume="TRUE"
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [ -z "$input_lif" ]; then
  echo "--input_lif is required" >&2
  usage
  exit 1
fi

upload="$(to_upper_bool "$upload")"
resume="$(to_upper_bool "$resume")"
outdir="${outdir:-${main_dir}/processed}"
logdir="${logdir:-${main_dir}/logs/wfOSFISH_LIF_${timestamp}}"
work_dir="${work_dir:-${main_dir}/work}"

mkdir -p "$main_dir" "$outdir" "$logdir" "$work_dir"
cd "$main_dir"

nextflow_args=(
  run "${script_dir}/wfOSFISH_LIF.nf"
  -c "${script_dir}/nextflow.config"
  -work-dir "$work_dir"
  --script_dir "$script_dir"
  --input_lif "$input_lif"
  --dataset_name "$dataset_name"
  --outdir "$outdir"
  --logdir "$logdir"
  --s3_bucket "$s3_bucket"
  --collection_table_url "$collection_table_url"
  --collection_table_sheet "$collection_table_sheet"
  --google_key "$google_key"
  --upload "$upload"
  --zarr_format "$zarr_format"
  --eubi_extra_args "$eubi_extra_args"
  --default_x_scale "$default_x_scale"
  --default_y_scale "$default_y_scale"
  --default_z_scale "$default_z_scale"
  --scale_unit "$scale_unit"
  -profile "$profile"
)

if [ "$resume" = "TRUE" ]; then
  nextflow_args+=("-resume")
fi

nextflow "${nextflow_args[@]}"
