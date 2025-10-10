#!/bin/bash

set -euo pipefail

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
conda activate /home/glbrc.org/mplace/miniforge3/envs/deeptools

IN_FOLDER="$1"
SAMPLE="$2"
OUT_FOLDER="$3"

echo "[$(date)] Starting wig conversion job for: $SAMPLE"

bamCoverage -p 8 -b "${IN_FOLDER}/${SAMPLE}.bam" -o "${OUT_FOLDER}/${SAMPLE}.bw"

