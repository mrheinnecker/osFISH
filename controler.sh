

wrkdir="/g/schwab/Marco/repos/osFISH"
container="/g/schwab/rheinnec/container_legacy/python_latest.sif"
logdir="/scratch/rheinnec/logs"


sf_run="SF02"


main_out_dir="/scratch/rheinnec/osFISH/$sf_run"
raw_img_dir="/g/schwab/Marco/projects/osFISH/SF/13052025/raw"




raw_string="M:02,05,08;mix:01;HR:03,06,09;Kmiki:04,07,10"



IFS=';' read -ra pair_array <<< "$raw_string"


for same_species in "${pair_array[@]}"; do

    echo $same_species

    prefix="${same_species%%:*}"
    numbers="${same_species#*:}"
    IFS=',' read -ra species_array <<< "$numbers"

    # Create the target directory
    mkdir -p "${raw_img_dir}/$prefix"

    # Loop over each species number
    for num in "${species_array[@]}"; do

        echo $num

        # Find matching files
        pattern="${sf_run}_${num}*"
        
        # Use globbing to match files
        for filepath in "${raw_img_dir}"/$pattern; do
            # Check if file actually exists (to avoid issues if no match)
            [ -e "$filepath" ] || continue
            
            # Get just the filename (without path)
            filename=$(basename "$filepath")
            
            # Create symbolic link in the prefix directory
            ln -s "$filepath" "${raw_img_dir}/$prefix/$filename"
        done
    done
done





mkdir $main_out_dir

ls $raw_img_dir -F | grep '/$' | while read subdir
do

input_dir=$raw_img_dir/$subdir

echo $input_dir

out_dir=$main_out_dir/$subdir
mkdir $out_dir

echo $out_dir

sbatch \
    -J "osFISH_$subdir" \
    -t 1:00:00 \
    --mem 32000 \
    -e "$logdir/log_osFISH_$subdir.txt" \
    -o "$logdir/out_osFISH_$subdir.txt" \
    --wrap="singularity exec --bind /g/schwab --bind /scratch $container python3 $wrkdir/prep_img.py --image_dir $input_dir --output_dir $out_dir"


done




#singularity exec --bind /g/schwab $container python3 $wrkdir/main.py
squeue -u rheinnec