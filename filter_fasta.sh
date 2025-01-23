#!/bin/bash

## run container 
# singularity shell --bind /media/rheinnec/OS ~/container/test.sif



full_database="/home/rheinnec/projects/osCLEM/pr2_version_5.0.0_SSU_taxo_long.fasta"
species_file="/home/rheinnec/projects/osCLEM/RCC_cultures_ordered.txt"
projdir="/home/rheinnec/projects/osCLEM"
outdir="/home/rheinnec/projects/osCLEM/species_seq"

wrkdir="/home/rheinnec/repos/osCLEM"

container="/mnt/schwab/Marco/container_legacy/oligoN_design.sif"

# Loop through each species and call the processing script
head -n 3 $species_file | while read spec
do
    echo "Processing species: $spec"
    singularity exec --bind /media/rheinnec/OS $container \
      $wrkdir/process_species.sh "$full_database" "$outdir" "$spec"
done






