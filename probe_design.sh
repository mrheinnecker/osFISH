#!/bin/bash

## run container 
# singularity shell --bind /media/rheinnec/OS ~/container/test.sif



full_database="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta"
species_file="/g/schwab/Marco/projects/osFISH/RCC_cultures_ordered_secondary.txt"
projdir="/g/schwab/Marco/projects/osFISH"
outdir="/scratch/rheinnec/osFISH/species"

wrkdir="/g/schwab/Marco/repos/osFISH"

container="/g/schwab/Marco/container_legacy/oligoN_design.sif"

mkdir $outdir

# Loop through each species and call the processing script
cat $species_file | while read spec
do
echo "Processing species: $spec"

sbatch \
    -J "probe_design_$spec" \
    -t 40:00:00 \
    --mem 2000 \
    -e "/scratch/rheinnec/logs/log_probedesign_$spec.txt" \
    -o "/scratch/rheinnec/logs/out_probedesign_$spec.txt" \
    $wrkdir/container.sh "$full_database" "$outdir" "$spec" "$wrkdir" "$container" 0.5 0.04
done


squeue -u rheinnec






for bp in 35 36 37 38 39 40 41 42
do
    echo "$bp"

bpdir="$outdir/bp$bp"

mkdir $bpdir

cat $species_file | while read spec
do
echo "Processing species: $spec"

sbatch \
    -J "$bp_$spec" \
    -t 20:00:00 \
    --mem 2000 \
    -e "/scratch/rheinnec/logs/log_probedesign_$bp_$spec.txt" \
    -o "/scratch/rheinnec/logs/out_probedesign_$bp_$spec.txt" \
    $wrkdir/container.sh "$full_database" "$bpdir" "$spec" "$wrkdir" "$container" 0.5 0.04 $bp
done



done







