#!/usr/bin/env bash
# Remove entries with no CDS from an accessions file.
# Usage: remove_no_cds.sh ACCESSIONS_FILE
set -euo pipefail

ACCESSIONS_FILE="$1"
PREFETCH_DIR="$2"

: "${PREFETCH_DIR:?PREFETCH_DIR not set}"
mkdir -p "$PREFETCH_DIR"

# Determine accession column index based on header
HEADER=$(head -n 1 "$ACCESSIONS_FILE")
IFS=$'\t' read -r -a HEADERS <<< "$HEADER"
ACC_IDX=-1
for i in "${!HEADERS[@]}"; do
  case "${HEADERS[$i]}" in
    "ACCESSION"|"Accession")
      ACC_IDX=$i
      break
      ;;
  esac
done

if [[ $ACC_IDX -lt 0 ]]; then
  echo "ACCESSION column not found" >&2
  exit 1
fi

TMP_FILE=$(mktemp)
echo "$HEADER" > "$TMP_FILE"

tail -n +2 "$ACCESSIONS_FILE" | while IFS=$'\t' read -r -a FIELDS; do
  ACC="${FIELDS[$ACC_IDX]}"
  cds_na="$PREFETCH_DIR/${ACC}_cds_na.fasta"
  cds_aa="$PREFETCH_DIR/${ACC}_cds_aa.fasta"
  genome="$PREFETCH_DIR/${ACC}.fasta"

  # Prefetch sequences if missing
  if [[ ! -s "$cds_na" || ! -s "$cds_aa" ]]; then
    efetch -db nucleotide -id "$ACC" -format fasta_cds_na > "$cds_na"
    efetch -db nucleotide -id "$ACC" -format fasta_cds_aa > "$cds_aa"
  fi
  if [[ ! -s "$genome" ]]; then
    efetch -db nucleotide -id "$ACC" -format fasta > "$genome"
  fi

  CDS_COUNT=$(grep -c '^>' "$cds_na" || true)
  if [[ $CDS_COUNT -gt 0 ]]; then
    (IFS=$'\t'; echo "${FIELDS[*]}") >> "$TMP_FILE"
  else
    echo "Removing $ACC - no CDS found" >&2
  fi
done

mv "$TMP_FILE" "$ACCESSIONS_FILE"
