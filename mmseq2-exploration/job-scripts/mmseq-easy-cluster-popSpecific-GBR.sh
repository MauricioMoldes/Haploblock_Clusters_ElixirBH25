#!/bin/sh
#PBS -N mmseq_easyClust_GBR
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -J 1-9

set -e
set -o pipefail

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate biohack2025_haplo

echo "start at $(date)"

export BASE_DIR=$HOME/BioHack2025/4-haploblocks/Haploblock_Clusters_ElixirBH25/mmseq2-exploration
export POP_ID=GBR

export INPUT_FASTA=$BASE_DIR/data/TNF/data/${POP_ID}/haploblock_phased_seq_TNFa/haploblock_phased_seq_merged/chr6_region_31480875-31598421.fa

# READ FROM A POPULATION SPECIFIC PARAMETER FILE HERE:
export POP_PARAM_FILE=$BASE_DIR/data/TNF/data/${POP_ID}/all_params_${POP_ID}.tsv

row=$((1 + PBS_ARRAY_INDEX)) #-- STARTS FROM 2 as NR=1 is header
read param_set min_seq cov cov_mode <<< $(awk -F'\t' -v r="$row" 'NR==r {print $1, $2, $3, $4}' "$POP_PARAM_FILE")

export PARAM_SET=$param_set
export MIN_SEQ_ID=$min_seq
export COV_FRACTION=$cov
export COV_MODE=$cov_mode

echo "min_seq_id ${MIN_SEQ_ID}"
echo "C ${COV_FRACTION}"
echo "cov-mode ${COV_MODE}"

export RESULTS=PARAM${PARAM_SET}_TNF_${POP_ID}_chr6_region_31480875-31598421
export TEMP_FOLDER=$EPHEMERAL/BioHack2025_tempfiles/$RESULTS

export RESULTS_DIR=$BASE_DIR/results/${POP_ID}/$RESULTS

mkdir -p $TEMP_FOLDER

mkdir -p $RESULTS_DIR
cd $RESULTS_DIR

mmseqs easy-cluster $INPUT_FASTA $RESULTS $TEMP_FOLDER \
            --min-seq-id $MIN_SEQ_ID -c $COV_FRACTION --cov-mode $COV_MODE

# ADD NUMBER OF CLUSTERS TO THE POPULATION SPECIFIC PARAMETER FILE DIRECTLY -- DIFFERENT SCRIPT

echo "end at $(date)"
echo "DONE!"