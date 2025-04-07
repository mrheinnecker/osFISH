#!/bin/bash


### run my r pipeline


full_database="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta"
sabeRprobes_container="/g/schwab/Marco/container_legacy/probeDesign_rtool.sif"


ls /g/schwab/Marco/projects/osFISH/test3 | while read spec
do
echo "Processing species: $spec"

target_file="/g/schwab/Marco/projects/osFISH/test2/$spec"

sbatch --job-name="$spec" \
    --time=10:00:00 \
    --mem=20000 \
    --cpus-per-task=10 \
    --error="/scratch/rheinnec/logs/log__paraprobedesign_$spec.txt" \
    --output="/scratch/rheinnec/logs/out_paraprobedesign_$spec.txt" \
    --wrap="singularity exec --bind /g/schwab --bind /scratch $sabeRprobes_container Rscript /g/schwab/Marco/repos/osFISH/full_probe_design.R -r '$full_database' -t $target_file -a '55' -b '75' -l '40' -m '3'"

done


squeue -u rheinnec


## run oligo-N-design

# ## run container 
# # singularity shell --bind /media/rheinnec/OS ~/container/test.sif



full_database="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta"
species_file="/g/schwab/Marco/projects/osFISH/flora_species.txt"
projdir="/g/schwab/Marco/projects/osFISH"
#outdir="/scratch/rheinnec/osFISH/flora_oligoNpipe"
outdir="/scratch/rheinnec/flora_oligo2"
logdir="$outdir/log"
wrkdir="/g/schwab/Marco/repos/osFISH"

container="/g/schwab/Marco/container_legacy/oligoN_design.sif"

mkdir $outdir

# Loop through each species and call the processing script
cat $species_file | while read spec
do
echo "Processing species: $spec"

sbatch \
    -J "$spec" \
    -t 00:01:00 \
    --mem 2000 \
    -e "$logdir/log_$spec.txt" \
    -o "$logdir/out_$spec.txt" \
    --wrap="singularity exec --bind /g/schwab --bind /scratch $container $wrkdir/process_species.sh "$full_database" "$outdir" "$spec" 0.5 0.05 "38-40""



done


squeue -u rheinnec




# ## oligo-N-design in length chunks

# for bp in 35 36 37 38 39 40 41 42
# do
#     echo "$bp"

# bpdir="$outdir/bp$bp"

# mkdir $bpdir

# cat $species_file | while read spec
# do
# echo "Processing species: $spec"

# sbatch \
#     -J "$bp_$spec" \
#     -t 20:00:00 \
#     --mem 2000 \
#     -e "/scratch/rheinnec/logs/log_probedesign_$bp_$spec.txt" \
#     -o "/scratch/rheinnec/logs/out_probedesign_$bp_$spec.txt" \
#     $wrkdir/container.sh "$full_database" "$bpdir" "$spec" "$wrkdir" "$container" 0.5 0.04 $bp
# done



# done




# ## probe design using sabeRprobes

# sabeRprobes_container="/g/schwab/Marco/container_legacy/probeDesign_rtool.sif"


# full_database="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta"
# species_file="/g/schwab/Marco/projects/osFISH/RCC_cultures_ordered.txt"
# projdir="/g/schwab/Marco/projects/osFISH"
# outdir="/scratch/rheinnec/osFISH/species_sabeRprobes"

# wrkdir="/g/schwab/Marco/repos/osFISH"

# container="/g/schwab/Marco/container_legacy/oligoN_design.sif"

# mkdir $outdir

# # Loop through each species and call the processing script
# cat $species_file | while read spec
# do
# echo "Processing species: $spec"

# sbatch \
#     -J "probe_design_$spec" \
#     -t 40:00:00 \
#     --mem 2000 \
#     -e "/scratch/rheinnec/logs/log_probedesign_$spec.txt" \
#     -o "/scratch/rheinnec/logs/out_probedesign_$spec.txt" \
#     singularity exec --bind /g/schwab --bind /scratch $container $wrkdir/process_species.sh "$full_database" "$outdir_raw" "$spec" "$abundance_target" "$abundance_ref" "$bp"

# done








### make speed test

full_database="/g/schwab/Marco/projects/osFISH/pr2_version_5.0.0_SSU_taxo_long.fasta"
sabeRprobes_container="/g/schwab/Marco/container_legacy/probeDesign_rtool.sif"

target_file="/g/schwab/Marco/projects/osFISH/test/Heterocapsa_rotundata_seq_all.fasta"

sbatch --job-name="mclapply" \
    --time=00:30:00 \
    --mem=20000 \
    --cpus-per-task=10 \
    --error="/home/rheinnec/logs/log_mcapply.txt" \
    --output="/home/rheinnec/logs/out_mcapply.txt" \
    --wrap="singularity exec --bind /g/schwab --bind /scratch $sabeRprobes_container Rscript /g/schwab/Marco/repos/osFISH/test_para_mclapply.R -r '$full_database' -t $target_file -a '55' -b '75' -l '40' -m '3'"



sbatch --job-name="lapply" \
    --time=00:30:00 \
    --mem=2000 \
    --cpus-per-task=10 \
    --error="/home/rheinnec/logs/log_lapply.txt" \
    --output="/home/rheinnec/logs/out_lapply.txt" \
    --wrap="singularity exec --bind /g/schwab --bind /scratch $sabeRprobes_container Rscript /g/schwab/Marco/repos/osFISH/test_para_lapply.R -r '$full_database' -t $target_file -a '55' -b '75' -l '40' -m '3'"











