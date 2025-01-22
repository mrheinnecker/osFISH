#!/bin/bash

## run container 
# singularity shell --bind /media/rheinnec/OS ~/container/test.sif


# Input files
full_database="/home/rheinnec/projects/osCLEM/pr2_version_5.0.0_SSU_taxo_long.fasta"


species_file="/home/rheinnec/projects/osCLEM/RCC_cultures_ordered.txt"

outdir="/home/rheinnec/projects/osCLEM/species_seq"


cat $species_file | while read spec 
do

echo $spec

sequenceSelect.py -f $full_database -o "${outdir}/${spec}.fasta" -p $spec -a k -v


done



# # Output file
# filtered_fasta="filtered_sequences.fasta"

# # Create a temporary regex file from the species list (accounting for minor spelling differences)
# # Convert spaces to match any whitespace and add case-insensitivity.
# awk '{print tolower($0)}' "$species_file" | \
# sed -E 's/ /\\s+/g' | \
# sed -E 's/(.*)/\(^|[>|_])\1(\\s|\$)/' > temp_species_patterns.txt

# # Filter FASTA file for matching headers using the generated regex
# grep -i -A 1 -f temp_species_patterns.txt "$fasta_file" | \
# grep -v -- "^--$" > "$filtered_fasta"

# # Clean up temporary file
# rm temp_species_patterns.txt

# # Completion message
# echo "Filtered sequences saved to: $filtered_fasta"

