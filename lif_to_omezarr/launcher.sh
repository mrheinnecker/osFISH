  cd /g/schwab/marco/repos/osFISH

  bash lif_to_omezarr/osfish_lif_main.sh interactive \
    --input_lif /g/schwab/marco/projects/taxseq/runs/confocal01/26-07-06/treatment.lif \
    --dataset_name block_screens_pond_26-06-26 \
    --main_dir /scratch/rheinnec/osFISH/lif_to_omezarr \
    --s3_bucket "s3embl/temscreen/blocks" \
    --collection_table_url "https://docs.google.com/spreadsheets/d/1vFMQKq8MDs3nURapyu6odc58IlkWsWXj0sKeb06NLnE/edit?gid=0#gid=0" \
    --collection_table_sheet "block_screens" \
    --google_key "/g/schwab/marco/repos/tem_classification/trec-tem-screen-e98a2e03f58b.json" \
    --upload TRUE \
    --resume FALSE