

wrkdir="/g/schwab/Marco/repos/osFISH"
container="/g/schwab/rheinnec/container_legacy/python_latest.sif"
logdir="/scratch/rheinnec/logs"


sf_run="SF02"


main_out_dir="/scratch/rheinnec/osFISH/$sf_run"
raw_img_dir="/g/schwab/Marco/projects/osFISH/SF/13052025/raw"




run_template="/g/schwab/Marco/projects/osFISH/image_run_templates/SF02.tsv"


mkdir $main_out_dir

cat "$run_template" | while IFS=$'\t' read -r prefix numbers channels species
do

    echo $prefix
    image_dir="${raw_img_dir}/$prefix"
    IFS=',' read -ra species_array <<< "$numbers"

    # Create the target directory
    mkdir -p "${raw_img_dir}/$prefix"

    # Loop over each species number
    for num_raw in "${species_array[@]}"; do

        num=$(printf "%02d" "$num_raw")

        # Find matching files
        pattern="${sf_run}_${num}*"
        
        # Use globbing to match files
        for filepath in "${raw_img_dir}"/$pattern; do
            # Check if file actually exists (to avoid issues if no match)
            [ -e "$filepath" ] || continue
            
            # Get just the filename (without path)
            filename=$(basename "$filepath")
            
            # Create symbolic link in the prefix directory
            ln -s "$filepath" "${image_dir}/$filename"
        done
    done

    
    out_dir=$main_out_dir/$prefix
    mkdir $out_dir

    echo "submitting clusterjob"


   # singularity exec --bind /g/schwab --bind /scratch $container python3 $wrkdir/prep_img.py --image_dir $image_dir --output_dir $out_dir --channels '$channels'

    sbatch \
        -J "osFISH_$prefix" \
        -t 1:00:00 \
        --mem 16000 \
        -e "$logdir/log_osFISH_$prefix.txt" \
        -o "$logdir/out_osFISH_$prefix.txt" \
        --wrap="singularity exec --bind /g/schwab --bind /scratch $container python3 $wrkdir/prep_img.py --image_dir $image_dir --output_dir $out_dir --channels '$channels'"



done





#singularity exec --bind /g/schwab $container python3 $wrkdir/main.py
squeue -u rheinnec