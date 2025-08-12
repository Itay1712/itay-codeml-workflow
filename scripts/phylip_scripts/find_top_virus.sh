#!/usr/bin/env bash
# Usage:   ./find_top_virus.sh ACCESSIONS_FILE ACCESSION_COL_INDEX
# Example: ./find_top_virus.sh all_viruses.tsv 0
#
# Dependencies: Entrez Direct (edirect: efetch, xtract), awk, grep
set -euo pipefail

ACCESSIONS_FILE="$1"
ACCESSION_COL_INDEX="$2"
PREFETCH_DIR="$3"

: "${PREFETCH_DIR:?PREFETCH_DIR not set}"
mkdir -p "$PREFETCH_DIR"

# Trackers for “best” virus
max_cds=0
best_len=0
best_acc=""

# We use process-substitution so variables aren’t in a subshell
while IFS=$'\t' read -ra FIELDS; do
    acc="${FIELDS[$ACCESSION_COL_INDEX]}"

    cds_file="$PREFETCH_DIR/${acc}_cds_na.fasta"
    genome_file="$PREFETCH_DIR/${acc}.fasta"

    # Retrieve sequences if not already cached
    if [[ ! -s "$cds_file" ]]; then
        efetch -db nucleotide -id "$acc" -format fasta_cds_na > "$cds_file"
    fi
    if [[ ! -s "$genome_file" ]]; then
        efetch -db nucleotide -id "$acc" -format fasta > "$genome_file"
    fi

    # 1) count CDS entries from cached file
    cds_count=$(grep -c '^>' "$cds_file")

    # 2) compute total genome length
    seq_len=$(awk '/^>/ {next} {total+=length($0)} END {print total}' "$genome_file")

    # 3) compare/update best
    if (( cds_count > max_cds )) || { (( cds_count == max_cds )) && (( seq_len > best_len )); }; then
        max_cds=$cds_count
        best_len=$seq_len
        best_acc=$acc
    fi

    # optional progress indicator
    printf "  %-12s  CDS=%3d  Len=%7d\n" "$acc" "$cds_count" "$seq_len" >&2
done < <(tail -n +2 "$ACCESSIONS_FILE")

echo "$best_acc"
