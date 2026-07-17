nextflow.enable.dsl=2

params.run_dir = params.run_dir ?: "/g/schwab/marco/projects/osFISH/runs/regFISH02"
params.raw_lif_dir = params.raw_lif_dir ?: "${params.run_dir}/2026-06-22"
params.outdir = params.outdir ?: "${params.run_dir}/omezarr"
params.extracted_tif_dir = params.extracted_tif_dir ?: "${params.run_dir}/extracted_tifs"
params.logdir = params.logdir ?: "${params.run_dir}/logs/wfOMEZARR"
params.script_dir = params.script_dir ?: baseDir.toString()
params.default_z_scale_nm = params.default_z_scale_nm ?: 1000
params.zarr_format = params.zarr_format ?: 2
params.eubi_extra_args = params.eubi_extra_args ?: ""
params.workflow_stage = params.workflow_stage ?: "all"
params.s3_bucket = params.s3_bucket ?: "s3embl/temscreen/osFISH"
params.sheet_mode = params.sheet_mode ?: "google"
params.google_key = params.google_key ?: "/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json"
params.collection_table_url = params.collection_table_url ?: "https://docs.google.com/spreadsheets/d/1vFMQKq8MDs3nURapyu6odc58IlkWsWXj0sKeb06NLnE/edit?gid=0#gid=0"
params.collection_table_sheet = params.collection_table_sheet ?: "ct"
params.metadata_root = params.metadata_root ?: new File(params.run_dir.toString()).getParent()


process EXTRACTTIFFROMLIF {
    cpus 1
    memory "128GB"
    time "8h"

    publishDir "${params.logdir}/extraction", mode:"copy", pattern:"*_extracted_tifs.tsv"
    containerOptions "--bind /g --bind /home --bind /scratch"
    errorStrategy "retry"
    maxRetries 1

    input:
    tuple val(run_name), val(run_date), val(condition), path(lif_path)

    output:
    path "*_extracted_tifs.tsv", emit: manifest

    script:
    """
    set -euo pipefail

    outdir="${params.extracted_tif_dir}/${run_name}/${run_date}/${condition}"
    mkdir -p "\$outdir"

    python3 "${params.script_dir}/extract_tifs_from_lif.py" \
      --lif "${lif_path}" \
      --outdir "\$outdir" \
      --manifest "${run_name}_${run_date}_${condition}_extracted_tifs.tsv" \
      --overwrite FALSE
    """
}

process BUILDCONVERSIONTABLE {
    cpus 1
    memory "2GB"
    time "20m"

    publishDir "${params.logdir}", mode:"copy"
    containerOptions "--bind /g --bind /home --bind /scratch"

    input:
    path extraction_manifests

    output:
    path "images_to_process.tsv", emit: to_process
    path "all_images.tsv", emit: all_images

    script:
    """
    set -euo pipefail

    python3 "${params.script_dir}/extract_lif_tif_metadata.py" \
      --raw-lif-dir "${params.raw_lif_dir}" \
      --run-dir "${params.run_dir}" \
      --all-output "all_images.tsv" \
      --process-output "images_to_process.tsv" \
      --default-z-scale-nm "${params.default_z_scale_nm}" \
      --extracted-tif-dir "${params.extracted_tif_dir}"
    """
}

process EUBITIFFTOOMEZARR {
    cpus 1
    memory { "${Math.min(Math.max((req_mem as Integer), 16), 128)}GB" }
    time "4h"

    publishDir "${params.outdir}/${condition}", mode:"copy"
    containerOptions "--bind /g --bind /home --bind /scratch"
    errorStrategy "ignore"
    maxRetries 1

    input:
    tuple val(condition), val(output_name), val(raw_tif_path), path(tif_path), val(x_scale_nm), val(y_scale_nm), val(z_scale_nm), val(req_mem)

    output:
    tuple val(condition), val(output_name), path("${output_name}.ome.zarr"), emit: omezarr
    path "${output_name}_conversion_done.txt"

    script:
    """
    set -euo pipefail

    rm -rf "${output_name}.ome.zarr"

    input_path="${raw_tif_path}"
    if [ ! -f "\$input_path" ]; then
      echo "Expected source TIFF does not exist: \$input_path" >&2
      ls -la >&2
      exit 1
    fi

    eubi_extra_args="${params.eubi_extra_args}"

    if [ -n "\$eubi_extra_args" ]; then
      eubi to_zarr \
        "\$input_path" \
        --output_path "${output_name}.ome.zarr" \
        --x_unit nm \
        --y_unit nm \
        --z_unit nm \
        --x_scale "${x_scale_nm}" \
        --y_scale "${y_scale_nm}" \
        --z_scale "${z_scale_nm}" \
        --save_omexml True \
        --autochunk True \
        --zarr_format "${params.zarr_format}" \
        --max_workers 1 \
        \$eubi_extra_args
    else
      eubi to_zarr \
        "\$input_path" \
        --output_path "${output_name}.ome.zarr" \
        --x_unit nm \
        --y_unit nm \
        --z_unit nm \
        --x_scale "${x_scale_nm}" \
        --y_scale "${y_scale_nm}" \
        --z_scale "${z_scale_nm}" \
        --save_omexml True \
        --autochunk True \
        --zarr_format "${params.zarr_format}" \
        --max_workers 1
    fi

    touch "${output_name}_conversion_done.txt"
    """
}


process S3UPLOADOSFISH {
    cpus 1
    memory "1GB"
    time "30m"

    publishDir "${params.logdir}/upload", mode:"copy"
    containerOptions "--bind /g --bind /home --bind /scratch"
    errorStrategy "retry"
    maxRetries 1

    input:
    tuple val(condition), val(output_name), path(omezarr)

    output:
    path "${output_name}_s3_upload_done.txt"

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
      find "${omezarr}" -maxdepth 3 -type f -o -type d >&2 || true
      exit 1
    fi

    mc cp "\${image_zarr}/" "${params.s3_bucket}/${output_name}.ome.zarr/" --recursive
    touch "${output_name}_s3_upload_done.txt"
    """
}

process COLLECTOSFISHS3FILES {
    cpus 1
    memory "1GB"
    time "10m"

    publishDir "${params.logdir}", mode:"copy"
    containerOptions "--bind /g --bind /home --bind /scratch"

    input:
    val trigger

    output:
    path "all_s3_entries.txt", emit: all_s3

    script:
    """
    set -euo pipefail
    mc ls "${params.s3_bucket}" > all_s3_entries.txt
    """
}

process MAKEOSFISHCOLLECTIONTABLE {
    cpus 1
    memory "1GB"
    time "10m"

    publishDir "${params.logdir}", mode:"copy"
    containerOptions "--bind /g --bind /home --bind /scratch"

    input:
    path all_s3
    path all_images

    output:
    path "done.tsv"
    path "osfish_collection_table.tsv", emit: collection_table

    script:
    """
    set -euo pipefail

    Rscript "${params.script_dir}/make_collection_table.R" \
      --all_s3 "${all_s3}" \
      --all_images "${all_images}" \
      --s3_bucket "${params.s3_bucket}" \
      --local_collection_table "osfish_collection_table.tsv" \
      --sheet_mode "${params.sheet_mode}" \
      --google_key "${params.google_key}" \
      --collection_table_url "${params.collection_table_url}" \
      --collection_table_sheet "${params.collection_table_sheet}" \
      --metadata_root "${params.metadata_root}"
    """
}

workflow {
    lif_ch = Channel
        .fromPath("${params.raw_lif_dir}/**/*.lif")
        .map { lif_file ->
            def relative_parts = params.run_dir ? lif_file.getParent().toString().replaceFirst("^" + java.util.regex.Pattern.quote(params.run_dir.toString()), "").tokenize("/") : []
            def run_date = relative_parts.find { it ==~ /\d{4}-\d{2}-\d{2}/ } ?: ""
            def condition = relative_parts.reverse().find { it.toLowerCase() in ["control", "treatment"] } ?: lif_file.baseName
            tuple(params.run_dir.tokenize("/")[-1], run_date, condition, lif_file)
        }

    EXTRACTTIFFROMLIF(lif_ch)
    extraction_done_ch = EXTRACTTIFFROMLIF.out.manifest.collect()
    BUILDCONVERSIONTABLE(extraction_done_ch)

    if (params.workflow_stage != "collection") {
        BUILDCONVERSIONTABLE.out.to_process
            .splitCsv(header:true, sep:'\t')
            .map { row ->
                tuple(
                    row.condition,
                    row.output_name,
                    row.tif_path,
                    file(row.tif_path),
                    row.x_scale_nm,
                    row.y_scale_nm,
                    row.z_scale_nm,
                    row.req_mem
                )
            }
            .set { conversion_ch }

        EUBITIFFTOOMEZARR(conversion_ch)

        if (params.workflow_stage == "all") {
            upload_done_ch = S3UPLOADOSFISH(EUBITIFFTOOMEZARR.out.omezarr).collect()
            COLLECTOSFISHS3FILES(upload_done_ch)
            MAKEOSFISHCOLLECTIONTABLE(
                COLLECTOSFISHS3FILES.out.all_s3,
                BUILDCONVERSIONTABLE.out.all_images
            )
        }
    }

    if (params.workflow_stage == "collection") {
        COLLECTOSFISHS3FILES(Channel.value("collection_table_only"))
        MAKEOSFISHCOLLECTIONTABLE(
            COLLECTOSFISHS3FILES.out.all_s3,
            BUILDCONVERSIONTABLE.out.all_images
        )
    }
}
