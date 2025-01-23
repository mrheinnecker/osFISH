#!/bin/bash

# Arguments
full_database="$1"
outdir_raw="$2"
spec="$3"
wrkdir="$4"

singularity exec --bind /g/schwab/rheinnec/ $container \
    $wrkdir/process_species.sh "$full_database" "$outdir" "$spec"