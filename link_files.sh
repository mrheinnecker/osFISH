

raw_img_dir="/g/schwab/Marco/projects/osFISH/SF/13052025/raw"

#raw_img_dir="/scratch/rheinnec/osFISH/image_analysis/SF02_01_01"


sf_run="SF02"


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



