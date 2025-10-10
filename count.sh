#!/bin/bash

set -euo pipefail

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
conda activate /home/glbrc.org/mplace/miniforge3/envs/plantrnaseq

IN_FOLDER=/home/glbrc.org/mplace/data/RNA-seq/sang_jin_kim/Analysis/map
SAMPLE_LIST=/mnt/bigdata/linuxhome/mplace/data/RNA-seq/sang_jin_kim/Analysis/map/sample_list.txt
OUT_FOLDER=/mnt/bigdata/linuxhome/mplace/data/RNA-seq/sang_jin_kim/Analysis
ANNOTATION=/mnt/bigdata/linuxhome/mplace/data/RNA-seq/sang_jin_kim/Analysis/map/Athaliana_167_gene.gff3

echo "[$(date)] Starting count job for samples: ${SAMPLE_LIST}"

# Filenames are used as matrix column headers
file_list=()
while IFS="" read -r sample || [ -n "$sample" ]
do
  file_list+=( "$sample" )
done < "${SAMPLE_LIST}"

# run featureCounts to get sample gene counts
featureCounts -s 2 -p --countReadPairs -B -C -T 32 -t gene -g ID -a "${ANNOTATION}" -o "${OUT_FOLDER}/counts.txt" "${file_list[@]}"

cut -f 1,6 "${OUT_FOLDER}/counts.txt" | csvtk tab2csv > gene_lengths.csv

cut --complement -f 2,3,4,5,6 "${OUT_FOLDER}/counts.txt" | csvtk transpose -t | csvtk tab2csv | rnanorm tpm --gene-lengths gene_lengths.csv | csvtk transpose > "${OUT_FOLDER}/tpm_counts.txt"
cut --complement -f 2,3,4,5,6 "${OUT_FOLDER}/counts.txt" | csvtk transpose -t | csvtk tab2csv | rnanorm fpkm --gene-lengths gene_lengths.csv | csvtk transpose > "${OUT_FOLDER}/fpkm_counts.txt"

