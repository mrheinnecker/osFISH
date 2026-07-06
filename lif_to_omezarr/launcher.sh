  cd /mnt/c/repos/osFISH

  bash lif_to_omezarr/osfish_lif_main.sh interactive \
    --input_lif /path/to/image.lif \
    --dataset_name test_osfish \
    --main_dir /scratch/rheinnec/osFISH/lif_to_omezarr \
    --s3_bucket "s3embl/imatrec/central_data_processing/osfish" \
    --collection_table_url "https://docs.google.com/spreadsheets/d/YOUR_SHEET_ID/edit?gid=0#gid=0" \
    --collection_table_sheet "osfish_collection_table" \
    --google_key "/g/schwab/marco/repos/osFISH/trec-tem-screen-e98a2e03f58b.json" \
    --upload TRUE