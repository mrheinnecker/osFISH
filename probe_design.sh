#!/bin/bash

## run container 
# singularity shell --bind /media/rheinnec/OS ~/container/test.sif



full_database="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta"
species_file="/g/schwab/Marco/projects/osFISH/RCC_cultures_ordered.txt"
projdir="/g/schwab/Marco/projects/osFISH"
outdir="/scratch/rheinnec/osFISH/species_seq_40bp"

wrkdir="/g/schwab/Marco/repos/osFISH"

container="/g/schwab/Marco/container_legacy/oligoN_design.sif"

mkdir $outdir

# Loop through each species and call the processing script
cat $species_file | while read spec
do
    echo "Processing species: $spec"

sbatch \
    -J "probe_design_$spec" \
    -t 30:00:00 \
    --mem 5000 \
    -e "/scratch/rheinnec/logs/log_probedesign_$spec.txt" \
    -o "/scratch/rheinnec/logs/out_probedesign_$spec.txt" \
    $wrkdir/container.sh "$full_database" "$outdir" "$spec" "$wrkdir" "$container" 0.8 0.001
done


squeue -u rheinnec



