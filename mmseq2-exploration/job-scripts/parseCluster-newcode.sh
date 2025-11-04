#!/usr/bin/env bash

export BASE_DIR=$HOME/BioHack2025/4-haploblocks/Haploblock_Clusters_ElixirBH25/mmseq2-exploration
export POP_ID=PUR

export OLD_TSV=$BASE_DIR/data/TNF/data/${POP_ID}/all_params_${POP_ID}.tsv

export RESULT_DIR=$BASE_DIR/results/$POP_ID/
export NEW_TSV=$RESULT_DIR/cluster_info_${POP_ID}.tsv

# read param_set min_seq cov cov_mode <<< $(awk -F'\t' -v r="$row" 'NR==r {print $1, $2, $3, $4}' "$POP_PARAM_FILE")

# Add header to new TSV
awk 'NR==1 {print $0"\tNUM_CLUSTER"}' "$OLD_TSV" > "$NEW_TSV"

# Iterate over data rows
tail -n +2 "$OLD_TSV" | while IFS=$'\t' read -r param_set min_seq cov cov_mode; do
    # Construct the folder path (example: POP_ID / PARAM0_INFO)
    folder=$RESULT_DIR/PARAM${param_set}_TNF_${POP_ID}_chr6_region_31480875-31598421
    
    # Adjust pattern if needed
    n=$(grep -c "^>" $folder/*_rep_seq.fasta)
    
    # Append the row + num_cluster to new TSV
    echo -e "$param_set\t$min_seq\t$cov\t$cov_mode\t$n" >> "$NEW_TSV"
done
