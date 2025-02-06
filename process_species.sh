#!/bin/bash

# Arguments
full_database="$1"
outdir_raw="$2"
spec="$3"
abundance_target="$4"
abundance_ref="$5"

outdir="${outdir_raw}/${spec}"

echo $outdir

echo $outdir_raw

mkdir $outdir

export TMPDIR="/scratch/rheinnec/tmp"

echo "Processing: $spec"

target_fasta="${outdir}/${spec}_target.fasta"
reference_fasta="${outdir}/${spec}_reference.fasta"

target_fasta_sl="${outdir}/${spec}_target_sl.fasta"
reference_fasta_sl="${outdir}/${spec}_reference_sl.fasta"

## Filter full database for relevant species
sequenceSelect.py -f "$full_database" -o "$target_fasta" -p "$spec" -a k -v
sequenceSelect.py -f "$full_database" -o "$reference_fasta" -p "$spec" -a r -v

## Make single-line FASTA file instead of multiline sequence
multi2linefasta.py -f "$target_fasta" -o "$target_fasta_sl"
multi2linefasta.py -f "$reference_fasta" -o "$reference_fasta_sl"

## Find oligos (time-intensive)
findOligo -t "$target_fasta_sl" -r "$reference_fasta_sl" -o "${outdir}/${spec}_probes" -l '18-22' -m $abundance_target -s $abundance_ref  

## Test oligos against the whole reference database
testOligo -r "$reference_fasta_sl" -p "${outdir}/${spec}_probes.fasta" -o "${outdir}/${spec}_probes_tested.tsv"  

## Check accessibility (two-step process)
alignOligo -t "$target_fasta_sl" -p "${outdir}/${spec}_probes.fasta" -o "${outdir}/${spec}_probes_aligned.fasta"


rateAccess -f "${outdir}/${spec}_probes_aligned.fasta" -o "${outdir}/${spec}_probes_access.tsv" 

## Merge logs and select probes
bindLogs -f "${outdir}/${spec}_probes.tsv" "${outdir}/${spec}_probes_tested.tsv" "${outdir}/${spec}_probes_access.tsv" -o "${outdir}/${spec}_probes_log.tsv" -r  
## -s = GC content minimum; -M = percentage of hits agains reference database allowing 2 mismatches; -b accesability brightness of 0.4
filterLog -l "${outdir}/${spec}_probes_log.tsv" -s "0.4" -M "0.005" -b "0.4"  
selectLog -l "${outdir}/${spec}_probes_log_filtered.tsv" -N "4"
