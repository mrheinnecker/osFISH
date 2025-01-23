#!/bin/bash

## run container 
# singularity shell --bind /media/rheinnec/OS ~/container/test.sif


# Input files
full_database="/home/rheinnec/projects/osCLEM/pr2_version_5.0.0_SSU_taxo_long.fasta"


species_file="/home/rheinnec/projects/osCLEM/RCC_cultures_ordered.txt"

projdir="/home/rheinnec/projects/osCLEM"

outdir="/home/rheinnec/projects/osCLEM/species_seq"


#cat $species_file | while read spec 
head -n 2 $species_file | while read spec 
do

spec="Akashiwo_sanguinea"

echo $spec

target_fasta="${outdir}/${spec}_target.fasta"
reference_fasta="${outdir}/${spec}_reference.fasta"

target_fasta_sl="${outdir}/${spec}_target_sl.fasta"
reference_fasta_sl="${outdir}/${spec}_reference_sl.fasta"

## filter full database for relevant species
sequenceSelect.py -f $full_database -o $target_fasta -p $spec -a k -v
sequenceSelect.py -f $full_database -o $reference_fasta -p $spec -a r -v
## make single line fasta file instead of multiline sequence
multi2linefasta.py -f $target_fasta -o $target_fasta_sl
multi2linefasta.py -f $reference_fasta -o $reference_fasta_sl


findOligo -t $target_fasta_sl -r $reference_fasta_sl -o probes -l '18-22' -m 0.8 -s 0.001  


done



filtered_fasta="~/projects/osCLEM/species_seq/Akashiwo_sanguinea.fasta"

single_line_fasta="~/projects/osCLEM/test/Akashiwo_sanguinea_single_line.fasta"


multi2linefasta.py -f $filtered_fasta -o $single_line_fasta




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

