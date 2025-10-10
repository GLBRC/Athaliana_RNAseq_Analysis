#!/bin/bash

set -euo pipefail

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
conda activate /home/glbrc.org/mplace/miniforge3/envs/plantrnaseq

IN_FOLDER="$1"
SAMPLE="$2"
OUT_FOLDER="$3"
INDEX_FOLDER="$4"
INDEX_D="$5"
INDEX_C="$6"

echo "[$(date)] Starting map job for: $SAMPLE"

mkdir -p "${OUT_FOLDER}/BAMs/"

# gsnap args
# -t num threads
# --split-output Basename for multiple-file output, separately for nomapping, paired_uniq, unpaired_uniq etc...
# --gunzip Uncompress gzipped input files
# -D Genome Directory (index dir) /mnt/bigdata/linuxhome/mplace/data/RNA-seq/sang_jin_kim/Reference/Athaliana_167/Athaliana_167_index
# -d Genome database /mnt/bigdata/linuxhome/mplace/data/RNA-seq/sang_jin_kim/Reference/Athaliana_167/Athaliana_167_index/Athaliana_167_fa
# -c /mnt/bigdata/linuxhome/mplace/data/RNA-seq/sang_jin_kim/Reference/Athaliana_167/Athaliana_167_index
# -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)
# -n -n, --npaths=INT               Maximum number of paths to print (default 100).
# -A, --format=STRING            Another format type, other than default. Currently implemented: sam
#

gsnap -t 32 --split-output "${SAMPLE}" --gunzip --transcriptome-mode=assist -D "${INDEX_FOLDER}" -d "${INDEX_D}" -c "${INDEX_C}" -N 1 -n 1 -A sam "${IN_FOLDER}/${SAMPLE}_R1.filter.fastq.gz" "${IN_FOLDER}/${SAMPLE}_R2.filter.fastq.gz"

samtools sort -O bam -o "${OUT_FOLDER}/BAMs/${SAMPLE}.bam" "${SAMPLE}.concordant_uniq"
samtools index --threads 8 -o "${OUT_FOLDER}/BAMs/${SAMPLE}.bam.bai"  "${OUT_FOLDER}/BAMs/${SAMPLE}.bam"
