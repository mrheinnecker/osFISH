#!/bin/bash

## run container 
# singularity shell --bind /media/rheinnec/OS ~/container/test.sif


# Input files
full_database="/home/rheinnec/projects/osCLEM/pr2_version_5.0.0_SSU_taxo_long.fasta"


species_file="/home/rheinnec/projects/osCLEM/RCC_cultures_ordered.txt"

projdir="/home/rheinnec/projects/osCLEM"

outdir="/home/rheinnec/projects/osCLEM/species_seq"


cat $species_file | while read spec 
#head -n 2 $species_file | while read spec 
do

#spec="Akashiwo_sanguinea"

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

## find oligos... not parallelized: time intensive
findOligo -t $target_fasta_sl -r $reference_fasta_sl -o "${outdir}/${spec}_probes" -l '18-22' -m 0.8 -s 0.001  

## test oligos again agiants whole reference database if the match other sequences if 1 or two mismatches are allowed... 
## i dont fuly understand what the output tells me
testOligo -r $reference_fasta_sl -p "${outdir}/${spec}_probes.fasta" -o "${outdir}/${spec}_probes_tested.tsv"  

## now the accesability is checked... two step process... first is fast
alignOligo -t $target_fasta_sl -p "${outdir}/${spec}_probes.fasta" -o "${outdir}/${spec}_probes_aligned.fasta"

## secod  step... fast... output also not fully understood but i get what it is trying to do
rateAccess -f "${outdir}/${spec}_probes_aligned.fasta" -o "${outdir}/${spec}_probes_access.tsv" 


## luckily there is a script which automatically merges all the info created before and lets you select the probes
## does basically: left_join by oligoN id
bindLogs -f "${outdir}/${spec}_probes.tsv" "${outdir}/${spec}_probes_tested.tsv" "${outdir}/${spec}_probes_access.tsv" -o "${outdir}/${spec}_probes_log.tsv" -r  

## this filters based on desired criteria
filterLog -l "${outdir}/${spec}_probes_log.tsv" -s "0.4" -M "0.005" -b "0.4"  

selectLog -l "${outdir}/${spec}_probes_log_filtered.tsv" -N "4"  

done

