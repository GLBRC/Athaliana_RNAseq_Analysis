#!/bin/bash

BASE_DIR="$1"
SAMPLE_LIST="$2"

echo "Summarizing samples from: ${SAMPLE_LIST}" >&2

echo -e "Sample\tRawFragments\tPassFilterFragments\tMappedFragments\tGeneAligned\t%PassFilter\t%Mapped\t%MappedRaw\t%GeneAligned\t%GeneAlignedRaw"

while IFS="" read -r SAMPLE || [ -n "$SAMPLE" ]
do
echo "- ${SAMPLE}" >&2

BBDUK_LOG="${BASE_DIR}/trim/logs/${SAMPLE}.trim.err"

RAW_READS=$(grep -E "Input:" "$BBDUK_LOG" | head -n 1 | awk '{print $2}')
FILTERED_READS=$(grep -E "Result:" "$BBDUK_LOG" | tail -n 1 | awk '{print $2}')

RAW_READS=$(( ${RAW_READS:-0} / 2 ))
FILTERED_READS=$(( ${FILTERED_READS:-0} / 2 ))


GSNAP_LOG="${BASE_DIR}/trim/logs/${SAMPLE}.map.err"
ALIGNED_INPUT=$(grep -E "Aligned" "$GSNAP_LOG" | cut -d " " -f 2)

if [ $ALIGNED_INPUT -ne $FILTERED_READS ]; then
  echo "${SAMPLE} counts differ between FILTERED READS and ALIGNMENT INPUT" >&2
fi

MAPPED_READS=$(grep -A 5 "file ${SAMPLE}" ${BASE_DIR}/count/logs/all.count.err | grep "Total" | cut -d ":" -f 2 | cut -f 2 -d " ")
COUNTED_READS=$(grep -A 5 "file ${SAMPLE}" ${BASE_DIR}/count/logs/all.count.err | grep "Success" | cut -d ":" -f 2 | cut -f 2 -d " ")


if (( $RAW_READS > 0 )); then
  PCT_FILTERED=$(awk "BEGIN {printf \"%.2f\", ($FILTERED_READS / $RAW_READS) * 100}")
  PCT_MAPPED=$(awk "BEGIN {printf \"%.2f\", ($MAPPED_READS / $FILTERED_READS) * 100}")
  PCT_MAPPED_RAW=$(awk "BEGIN {printf \"%.2f\", ($MAPPED_READS / $RAW_READS) * 100}")
  PCT_GENE=$(awk "BEGIN {printf \"%.2f\", ($COUNTED_READS / $MAPPED_READS) * 100}")
  PCT_GENE_RAW=$(awk "BEGIN {printf \"%.2f\", ($COUNTED_READS / $RAW_READS) * 100}")
else
  echo "ERROR: No RAW READS found for sample ${SAMPLE}" >&2
  PCT_FILTERED="0.0"
  PCT_MAPPED="0.0"
  PCT_MAPPED_RAW="0.0"
  PCT_GENE="0.0"
  PCT_GENE_RAW="0.0"
fi

echo -e "${SAMPLE}\t${RAW_READS}\t${FILTERED_READS}\t${MAPPED_READS}\t${COUNTED_READS}\t${PCT_FILTERED}\t${PCT_MAPPED}\t${PCT_MAPPED_RAW}\t${PCT_GENE}\t${PCT_GENE_RAW}"

done < "${SAMPLE_LIST}"


