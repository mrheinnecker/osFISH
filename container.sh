#!/bin/bash

# Arguments
full_database="$1"
outdir_raw="$2"
spec="$3"
wrkdir="$4"
container="$5"
abundance_target="$6"
abundance_ref="$7"

echo "running container"
echo 


singularity exec --bind /g/schwab --bind /scratch $container \
    $wrkdir/process_species.sh "$full_database" "$outdir_raw" "$spec" "$abundance_target" "$abundance_ref"