container="/g/schwab/marco/container_legacy/probeDesign_rtool.sif"  
script="/g/schwab/marco/repos/osFISH/one_time_jobs/acarterae_i28.R"  
sbatch \
        -J "osFISH_aln" \
        -t 10:00:00 \
        --mem 5000 \
        --ntasks-per-node 8 \
        -e "/home/rheinnec/aln.log" \
        -o "/home/rheinnec/aln.txt" \
        --wrap="singularity exec --bind /g/schwab --bind /scratch -B /home/rheinnec/R/x86_64-pc-linux-gnu-library/4.4 $container Rscript $script"

