params.input_lif = params.input_lif ?: ""
params.dataset_name = params.dataset_name ?: ""
params.script_dir = params.script_dir ?: baseDir.toString()
params.outdir = params.outdir ?: "osfish_omezarr"
params.logdir = params.logdir ?: "osfish_logs"
params.s3_bucket = params.s3_bucket ?: "s3embl/osfish"
params.collection_table_url = params.collection_table_url ?: ""
params.collection_table_sheet = params.collection_table_sheet ?: "osfish_collection_table"
params.google_key = params.google_key ?: "${params.script_dir}/trec-tem-screen-e98a2e03f58b.json"
params.upload = params.upload ?: "FALSE"
params.zarr_format = params.zarr_format ?: 2
params.default_x_scale = params.default_x_scale ?: ""
params.default_y_scale = params.default_y_scale ?: ""
params.default_z_scale = params.default_z_scale ?: ""
params.scale_unit = params.scale_unit ?: "nm"
params.eubi_extra_args = params.eubi_extra_args ?: ""


def sanitizeDatasetName(value) {
    def base = new File(value.toString()).getName()
    base = base.replaceFirst(/(?i)\.lif$/, "")
    base = base.replaceAll(/[^A-Za-z0-9]+/, "_")
    base = base.replaceAll(/^_+|_+$/, "")
    return base ?: "osfish_lif"
}


process EXTRACTOSFISHLIFMETADATA {

    cpus 1
    memory "2GB"
    time "20m"

    publishDir "${params.logdir}/metadata", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    tuple val(dataset_name), val(input_lif)

    output:
    tuple val(dataset_name), val(input_lif), path("${dataset_name}_metadata.json"), path("${dataset_name}_pixel_size.tsv"), emit: metadata
    path "${dataset_name}_metadata.json", emit: metadata_json

    script:
    """
    set -euo pipefail

    python3 "${params.script_dir}/extract_lif_metadata.py" \
      --input "${input_lif}" \
      --name "${dataset_name}" \
      --metadata-json "${dataset_name}_metadata.json" \
      --pixel-size-tsv "${dataset_name}_pixel_size.tsv" \
      --x-scale "${params.default_x_scale}" \
      --y-scale "${params.default_y_scale}" \
      --z-scale "${params.default_z_scale}" \
      --scale-unit "${params.scale_unit}"
    """
}


process CONVERTOSFISHLIFTOOMEZARR {

    cpus 1
    memory "32GB"
    time "2h"

    publishDir "${params.outdir}", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    tuple val(dataset_name), val(input_lif), path(metadata_json), path(pixel_size_tsv)

    output:
    tuple val(dataset_name), path("${dataset_name}.ome.zarr"), path(metadata_json), emit: omezarr
    path "${dataset_name}_conversion_done.txt"

    script:
    """
    set -euo pipefail

    pixel_scale_x=\$(awk 'NR==2 {print \$1}' "${pixel_size_tsv}")
    pixel_scale_y=\$(awk 'NR==2 {print \$2}' "${pixel_size_tsv}")
    pixel_scale_z=\$(awk 'NR==2 {print \$3}' "${pixel_size_tsv}")

    if [ -z "\$pixel_scale_x" ] || [ -z "\$pixel_scale_y" ] || [ -z "\$pixel_scale_z" ]; then
      echo "Missing x/y/z scale values for ${dataset_name}; add metadata to the file or pass fallback scales." >&2
      exit 1
    fi

    rm -rf "${dataset_name}.ome.zarr"

    extra_args=()
    eubi_extra_args="${params.eubi_extra_args}"
    case "\${eubi_extra_args}" in
      ""|TRUE|true|FALSE|false|0|1|yes|YES|no|NO)
        eubi_extra_args=""
        ;;
    esac
    if [ -n "\${eubi_extra_args}" ]; then
      if [[ "\${eubi_extra_args}" != --* ]]; then
        echo "Ignoring eubi_extra_args because it does not look like CLI flags: \${eubi_extra_args}" >&2
      else
        read -r -a extra_args <<< "\${eubi_extra_args}"
      fi
    fi

    echo "Input LIF: ${input_lif}"
    echo "Output OME-Zarr: ${dataset_name}.ome.zarr"
    echo "Extra EuBI args: \${extra_args[*]:-<none>}"

    eubi to_zarr \
      "${input_lif}" \
      "${dataset_name}.ome.zarr" \
      --x_unit nm \
      --y_unit nm \
      --z_unit nm \
      --x_scale "\${pixel_scale_x}" \
      --y_scale "\${pixel_scale_y}" \
      --z_scale "\${pixel_scale_z}" \
      --dimension_order xyzct \
      --squeeze True \
      --save_omexml True \
      --zar_format "${params.zarr_format}" \
      --auto_chunk True \
      --jvm_memory 8GB \
      --max_workers 1 \
      "\${extra_args[@]}"

    touch "${dataset_name}_conversion_done.txt"
    """
}


process PATCHOSFISHOMEZARRMETADATA {

    cpus 1
    memory "1GB"
    time "10m"

    publishDir "${params.logdir}/omezarr_metadata", mode:"copy", pattern:"*_omezarr_metadata.tsv"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    tuple val(dataset_name), path(omezarr), path(metadata_json)

    output:
    tuple val(dataset_name), path(omezarr), path(metadata_json), emit: patched_omezarr
    path "${dataset_name}_omezarr_metadata.tsv"

    script:
    """
    set -euo pipefail

    python3 "${params.script_dir}/patch_omezarr_metadata.py" \
      --omezarr "${omezarr}" \
      --metadata-json "${metadata_json}" \
      --log "${dataset_name}_omezarr_metadata.tsv"
    """
}


process UPLOADOSFISHOMEZARR {

    cpus 1
    memory "1GB"
    time "30m"

    publishDir "${params.logdir}/upload", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    tuple val(dataset_name), path(omezarr), path(metadata_json)

    output:
    tuple val(dataset_name), path(metadata_json), path("${dataset_name}_s3_upload_done.txt"), emit: uploaded

    script:
    """
    set -euo pipefail

    image_zarr="${omezarr}"
    if [ ! -e "\${image_zarr}/.zattrs" ] && [ ! -e "\${image_zarr}/.zgroup" ]; then
      for candidate in "\${image_zarr}"/*.zarr "\${image_zarr}"/*.ome.zarr "\${image_zarr}"/*/*.zarr "\${image_zarr}"/*/*.ome.zarr; do
        if [ -d "\$candidate" ]; then
          image_zarr="\$candidate"
          break
        fi
      done
    fi

    if [ ! -e "\${image_zarr}/.zattrs" ] && [ ! -e "\${image_zarr}/.zgroup" ]; then
      echo "Could not find an OME-Zarr root marker under ${omezarr}" >&2
      exit 1
    fi

    mc cp "\${image_zarr}/" "${params.s3_bucket}/${dataset_name}.ome.zarr/" --recursive
    touch "${dataset_name}_s3_upload_done.txt"
    """
}


process MAKEOSFISHMOBIETABLE {

    cpus 1
    memory "1GB"
    time "10m"

    publishDir "${params.logdir}", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    tuple val(dataset_name), path(metadata_json)

    output:
    path "mobie_collection_table.tsv"

    script:
    """
    set -euo pipefail

    python3 "${params.script_dir}/make_mobie_collection_table.py" \
      --metadata-json "${metadata_json}" \
      --dataset-name "${dataset_name}" \
      --s3-bucket "${params.s3_bucket}" \
      --output "mobie_collection_table.tsv"
    """
}


process UPLOADOSFISHMOBIETABLE {

    cpus 1
    memory "1GB"
    time "10m"

    publishDir "${params.logdir}", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    path collection_table

    output:
    path "collection_table_upload_done.txt"

    script:
    """
    set -euo pipefail

    Rscript "${params.script_dir}/upload_collection_table.R" \
      --collection_table "${collection_table}" \
      --google_key "${params.google_key}" \
      --collection_table_url "${params.collection_table_url}" \
      --collection_table_sheet "${params.collection_table_sheet}"
    """
}


workflow {
    if (!params.input_lif) {
        error "Missing required parameter: --input_lif /path/to/file.lif"
    }

    dataset_name = params.dataset_name ?: sanitizeDatasetName(params.input_lif)
    input_ch = Channel.value(tuple(dataset_name, params.input_lif))

    EXTRACTOSFISHLIFMETADATA(input_ch)
    CONVERTOSFISHLIFTOOMEZARR(EXTRACTOSFISHLIFMETADATA.out.metadata)
    PATCHOSFISHOMEZARRMETADATA(CONVERTOSFISHLIFTOOMEZARR.out.omezarr)

    upload_enabled = params.upload.toString().toLowerCase() in ["true", "1", "yes"]
    if (upload_enabled) {
        UPLOADOSFISHOMEZARR(PATCHOSFISHOMEZARRMETADATA.out.patched_omezarr)
        table_input = UPLOADOSFISHOMEZARR.out.uploaded.map { name, metadata_json, done -> tuple(name, metadata_json) }
        MAKEOSFISHMOBIETABLE(table_input)
    } else {
        table_input = PATCHOSFISHOMEZARRMETADATA.out.patched_omezarr.map { name, omezarr, metadata_json -> tuple(name, metadata_json) }
        MAKEOSFISHMOBIETABLE(table_input)
    }

    if (params.collection_table_url) {
        UPLOADOSFISHMOBIETABLE(MAKEOSFISHMOBIETABLE.out)
    }
}
