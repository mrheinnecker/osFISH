params.input_lif = params.input_lif ?: ""
params.input_table = params.input_table ?: ""
params.dataset_name = params.dataset_name ?: ""
params.script_dir = params.script_dir ?: baseDir.toString()
params.outdir = params.outdir ?: "osfish_omezarr"
params.logdir = params.logdir ?: "osfish_logs"
params.uri_from = params.uri_from ?: ""
params.uri_to = params.uri_to ?: ""
params.collection_table_url = params.collection_table_url ?: ""
params.collection_table_sheet = params.collection_table_sheet ?: "osfish_collection_table"
params.google_key = params.google_key ?: "${params.script_dir}/trec-tem-screen-e98a2e03f58b.json"
params.zarr_format = params.zarr_format ?: 2
params.default_x_scale = params.default_x_scale ?: ""
params.default_y_scale = params.default_y_scale ?: ""
params.default_z_scale = params.default_z_scale ?: ""
params.scale_unit = params.scale_unit ?: "nm"
params.eubi_extra_args = params.eubi_extra_args ?: ""


def sanitizeDatasetName(value) {
    def base = new File(value.toString()).getName()
    base = base.replaceFirst(/(?i)(\.ome)?\.(lif|czi|tiff?)$/, "")
    base = base.replaceAll(/[^A-Za-z0-9]+/, "_")
    base = base.replaceAll(/^_+|_+$/, "")
    return base ?: "osfish_lif"
}


process EXTRACTOSFISHLIFMETADATA {

    cpus 1
    memory "2GB"
    time "20m"
    errorStrategy "ignore"
    tag "${dataset_name}"

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


process CONVERTOSFISHLIFSCENESTOOMEZARR {

    cpus 1
    memory { 32.GB * task.attempt }
    time { 4.h * task.attempt }
    errorStrategy { task.attempt < 3 ? "retry" : "ignore" }
    maxRetries 2
    tag "${dataset_name}"

    publishDir "${params.outdir}", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    tuple val(dataset_name), val(input_lif), path(metadata_json), path(pixel_size_tsv)

    output:
    tuple val(dataset_name), path("*.ome.zarr"), path("scene_metadata"), emit: omezarr
    path "${dataset_name}_conversion_done.txt"

    script:
    """
    set -euo pipefail

    pixel_scale_x=\$(awk 'NR==2 {print \$1}' "${pixel_size_tsv}")
    pixel_scale_y=\$(awk 'NR==2 {print \$2}' "${pixel_size_tsv}")
    pixel_scale_z=\$(awk 'NR==2 {print \$3}' "${pixel_size_tsv}")

    rm -rf *.ome.zarr scene_inputs scene_metadata scenes_to_convert.tsv

    python3 "${params.script_dir}/export_lif_scenes.py" \
      --input "${input_lif}" \
      --base-name "${dataset_name}" \
      --metadata-json "${metadata_json}" \
      --scene-dir "scene_inputs" \
      --scene-metadata-dir "scene_metadata" \
      --scene-table "scenes_to_convert.tsv"

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
    echo "Scene table:"
    cat scenes_to_convert.tsv
    echo "Extra EuBI args: \${extra_args[*]:-<none>}"

    tail -n +2 scenes_to_convert.tsv | while IFS=\$'\t' read -r scene_dataset scene_index scene_name scene_ome_tiff scene_metadata scene_x scene_y scene_z; do
      echo "Converting scene \${scene_index}: \${scene_name} -> \${scene_dataset}.ome.zarr"

      # Scenes in one LIF can have different pixel sizes; LIF-level values are the fallback.
      scene_x="\${scene_x:-\$pixel_scale_x}"
      scene_y="\${scene_y:-\$pixel_scale_y}"
      scene_z="\${scene_z:-\$pixel_scale_z}"
      if [ -z "\$scene_x" ] || [ -z "\$scene_y" ] || [ -z "\$scene_z" ]; then
        echo "Missing x/y/z scale values for \${scene_dataset}; add metadata to the file or pass fallback scales." >&2
        exit 1
      fi

      eubi to_zarr \
        "\${scene_ome_tiff}" \
        "\${scene_dataset}.ome.zarr" \
        --x_unit nm \
        --y_unit nm \
        --z_unit nm \
        --x_scale "\${scene_x}" \
        --y_scale "\${scene_y}" \
        --z_scale "\${scene_z}" \
        --dimension_order xyzct \
        --squeeze True \
        --save_omexml True \
        --zar_format "${params.zarr_format}" \
        --auto_chunk True \
        --jvm_memory 8GB \
        --max_workers 1 \
        "\${extra_args[@]}"

      python3 "${params.script_dir}/patch_omezarr_metadata.py" \
        --omezarr "\${scene_dataset}.ome.zarr" \
        --metadata-json "\${scene_metadata}" \
        --log "\${scene_dataset}_omezarr_metadata.tsv"
    done

    touch "${dataset_name}_conversion_done.txt"
    """
}


process MAKEOSFISHMOBIETABLE {

    cpus 1
    memory "1GB"
    time "10m"

    publishDir "${params.logdir}", mode:"copy"
    publishDir "${params.outdir}", mode:"copy"
    containerOptions "--bind /g --bind /scratch --bind /home"

    input:
    path(metadata_dirs, stageAs: "scene_metadata_*")

    output:
    path "mobie_collection_table.tsv"

    script:
    inventory_arg = params.input_table ? "--inventory \"${params.input_table}\"" : ""
    """
    set -euo pipefail

    python3 "${params.script_dir}/make_mobie_collection_table.py" \
      --metadata-json scene_metadata_* \
      ${inventory_arg} \
      --data-dir "${file(params.outdir).toAbsolutePath()}" \
      --uri-from "${params.uri_from}" \
      --uri-to "${params.uri_to}" \
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
    if (params.input_table) {
        // One row per LIF/CZI/OME-TIFF, as written by build_run_inventory.py.
        input_ch = Channel
            .fromPath(params.input_table)
            .splitCsv(header: true, sep: "\t")
            .map { row -> tuple(row.dataset_name, row.input_file) }
    } else if (params.input_lif) {
        dataset_name = params.dataset_name ?: sanitizeDatasetName(params.input_lif)
        input_ch = Channel.value(tuple(dataset_name, params.input_lif))
    } else {
        error "Missing required parameter: --input_lif /path/to/file.lif or --input_table inventory.tsv"
    }

    EXTRACTOSFISHLIFMETADATA(input_ch)
    CONVERTOSFISHLIFSCENESTOOMEZARR(EXTRACTOSFISHLIFMETADATA.out.metadata)

    metadata_dirs = CONVERTOSFISHLIFSCENESTOOMEZARR.out.omezarr.map { name, omezarr, metadata_dir -> metadata_dir }
    // One collection table for all LIFs, written after every conversion has finished.
    MAKEOSFISHMOBIETABLE(metadata_dirs.collect())

    if (params.collection_table_url) {
        UPLOADOSFISHMOBIETABLE(MAKEOSFISHMOBIETABLE.out)
    }
}
