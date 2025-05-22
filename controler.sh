

wrkdir="/g/schwab/marco/repos/osFISH"
container="/g/schwab/marco/container_legacy/python_latest.sif"
logdir="/scratch/rheinnec/logs"


sf_run="SF02"
anal_run="SF02_2205_2"

main_out_dir="/scratch/rheinnec/osFISH/$anal_run"
raw_img_dir="/g/schwab/marco/projects/osFISH/SF/13052025/raw"




run_template="/g/schwab/marco/projects/osFISH/image_run_templates/SF02.tsv"


mkdir $main_out_dir

tail -n +2  "$run_template" | while IFS=$'\t' read -r prefix numbers channels species bf_scaling dapi_scaling
do
# echo $channels
# done
    echo $prefix
    image_dir="${raw_img_dir}/$prefix"
    IFS=',' read -ra species_array <<< "$numbers"

    # Create the target directory
    mkdir -p "${raw_img_dir}/$prefix"

    # Loop over each species number
    for num_raw in "${species_array[@]}"; do
        #echo $num_raw
        num=$(printf "%02d" "$num_raw")

        # Find matching files
        pattern="${sf_run}_${num}*"
        
        # Use globbing to match files
        for filepath in "${raw_img_dir}"/$pattern; do
            echo $filepath
            # Check if file actually exists (to avoid issues if no match)
            [ -e "$filepath" ] || continue
            
            # Get just the filename (without path)
            filename=$(basename "$filepath")
            
            echo "linking to: $filepath from ${image_dir}/$filename" 

            # Create symbolic link in the prefix directory
            ln -s "$filepath" "${image_dir}/$filename"
        done
    done
#done
    
    out_dir=$main_out_dir/$prefix
    mkdir $out_dir

    echo "submitting clusterjob"


   #singularity exec --bind /g/schwab --bind /scratch $container python3 $wrkdir/prep_img.py --image_dir $image_dir --output_dir $out_dir --channels "${channels}" --dapi_scaling "${dapi_scaling}" --bf_scaling "${bf_scaling}"
    sbatch \
        -J "osFISH_$prefix" \
        -t 0:30:00 \
        --mem 32000 \
        -e "$logdir/log_osFISH_$prefix.txt" \
        -o "$logdir/out_osFISH_$prefix.txt" \
        --wrap="singularity exec --bind /g/schwab --bind /scratch $container python3 $wrkdir/prep_img.py --image_dir $image_dir --output_dir $out_dir --channels '$channels' --dapi_scaling '$dapi_scaling' --bf_scaling '$bf_scaling'"



done


## run manual:


prefix="mix_sub01"
image_dir="${raw_img_dir}/$prefix"
out_dir=$main_out_dir/$prefix
mkdir $out_dir
channels="Bright, DAPI, Cy3, Cy5, At590"
dapi_scaling="99.9"
bf_scaling="99.9"


sbatch \
    -J "osFISH_$prefix" \
    -t 0:30:00 \
    --mem 32000 \
    -e "$logdir/log_osFISH_$prefix.txt" \
    -o "$logdir/out_osFISH_$prefix.txt" \
    --wrap="singularity exec --bind /g/schwab --bind /scratch $container python3 $wrkdir/prep_img.py --image_dir $image_dir --output_dir $out_dir --channels '$channels' --dapi_scaling '$dapi_scaling' --bf_scaling '$bf_scaling'"






#singularity exec --bind /g/schwab $container python3 $wrkdir/main.py
squeue -u rheinnec