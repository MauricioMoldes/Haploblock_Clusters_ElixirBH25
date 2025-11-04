#!/bin/sh
#PBS -N mmseq_easyClust_CHB
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

set -e
set -o pipefail

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate biohack2025_haplo
echo "start at $(date)"

export BASE_DIR=$HOME/BioHack2025/4-haploblocks/Haploblock_Clusters_ElixirBH25/mmseq2-exploration
export POP_ID=PUR

export INPUT_PARAMS_FILE=$BASE_DIR/data/TNF/data/${POP_ID}/haploblock_phased_seq_TNFa/TNFa_${POP_ID}_variant_counts.tsv

# Read the TSV
read start end mean stdev < <(awk -F'\t' 'NR==2 {print $1, $2, $3, $4}' "$INPUT_PARAMS_FILE")

# Compute end - start
haplo_length=$((end - start))

lower=$(echo "scale=6; 1 - ($mean - $stdev) / $haplo_length" | bc -l)
mid=$(echo "scale=6; 1 - ($mean) / $haplo_length" | bc -l)
upper=$(echo "scale=6; 1 - ($mean + $stdev) / $haplo_length" | bc -l)

# Store in an array
seq_id_array=($lower $mid $upper)
cov_array=(0.999 0.9942 0.9999)

# READ FROM A POPULATION SPECIFIC PARAMETER FILE HERE:
export POP_PARAM_FILE=$BASE_DIR/data/TNF/data/${POP_ID}/all_params_${POP_ID}.tsv

# Convert arrays to space-separated strings for awk
seq_str="${seq_id_array[*]}"
cov_str="${cov_array[*]}"

# Use awk to generate the Cartesian product and add header
awk -v seq="$seq_str" -v cov="$cov_str" 'BEGIN{
    print "PARAM_SET\tMIN_SEQ_ID\tCOV_FRACTION\tCOV_MODE"   # header
    n = split(seq, s_arr, " ")
    m = split(cov, c_arr, " ")
    count = 0
    for(i=1;i<=n;i++){
        for(j=1;j<=m;j++){
            print count "\t" s_arr[i] "\t" c_arr[j] "\t" 0
            count += 1
        }
    }
}' > "$POP_PARAM_FILE"

echo "end at $(date)"
echo "DONE!"