#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/qtlseq_out_$(date +%Y%m%d_%H%M%S)"

qtlseq -o "${OUT_DIR}" \
       -n1 20 \
       -n2 20 \
       -w 100 \
       -s 20 \
       -r "${SCRIPT_DIR}/qtlseq_ref.fasta" \
       -p "${SCRIPT_DIR}/qtlseq_parent.1.fastq.gz,${SCRIPT_DIR}/qtlseq_parent.2.fastq.gz" \
       -b1 "${SCRIPT_DIR}/qtlseq_bulk1.1.fastq.gz,${SCRIPT_DIR}/qtlseq_bulk1.2.fastq.gz" \
       -b2 "${SCRIPT_DIR}/qtlseq_bulk2.1.fastq.gz,${SCRIPT_DIR}/qtlseq_bulk2.2.fastq.gz"

echo "QTL-seq test completed. Output: ${OUT_DIR}"
