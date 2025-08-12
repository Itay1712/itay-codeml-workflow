#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


GLOBALS="$1"
source "$GLOBALS"
FASTA="$2"
OUTDIR="$3"
BASE="${4:-$(basename "$FASTA" .fasta)}"

if [[ ! -f "$FASTA" ]]; then
    echo "!!! ERROR: FASTA file not found: $FASTA"
    exit 1
fi


run_3seq_on_fasta() {
    local fasta_file="$1"
    local outdir="$2"
    local base_label="$3"
    local global_mask_file="$outdir/recombination_regions.mask.tsv"

    # Remove identical sequences and ensure at least 3 unique sequences remain
    local dedup_fasta
    dedup_fasta=$(mktemp)
    seqkit rmdup -s "$fasta_file" > "$dedup_fasta"
    local num_unique
    num_unique=$(grep -c '^>' "$dedup_fasta")
    if (( num_unique < 3 )); then
        echo "Skipping $base_label — only $num_unique unique sequences (need ≥ 3)"
        rm -f "$dedup_fasta"
        return 0
    fi

    rm -f 3s.log 3s.pvalHist 3s.rec.csv 3s.longRec

    # Run 3SEQ quietly; when it fails to generate outputs (e.g., <3 unique sequences)
    # we skip downstream processing instead of erroring out
    3seq -y -c "$PV_TABLE_FILE" >/dev/null 2>&1 || true
    3seq -f "$dedup_fasta" -p "$PV_TABLE_FILE" -d >/dev/null 2>&1 <<< 'Y' || true

    if [[ ! -f 3s.pvalHist ]]; then
        echo "Skipping $base_label — 3SEQ output not generated (likely <3 unique sequences)"
        [[ -f 3s.log ]] && mv "3s.log" "$outdir/3s-${base_label}.log"
        rm -f "$dedup_fasta"
        return 0
    fi

    [[ -f 3s.log ]] && mv "3s.log" "$outdir/3s-${base_label}.log"
    mv "3s.pvalHist" "$outdir/3s-${base_label}.pvalHist"

    [[ -f 3s.rec.csv ]] && mv "3s.rec.csv" "$outdir/3s-${base_label}.rec.csv"

    if [[ -f 3s.longRec ]]; then
        mv "3s.longRec" "$outdir/3s-${base_label}.longRec"
        echo "✅ Recombination found in: $base_label"

        # Extract recombinant regions and append to global mask file
        awk '
          /\[.*\]/ {
              match($0, /\[([0-9]+)-([0-9]+)\]/, coords)
              if (coords[1] && coords[2]) {
                  header = $1
                  print header "\t" coords[1] "\t" coords[2]
              }
          }
        ' "$outdir/3s-${base_label}.longRec" >> "$global_mask_file"

    else
        echo "ℹ️ No recombination detected in: $base_label"
    fi

    rm -f "$dedup_fasta"
}

mkdir -p "$OUTDIR"

# === Detect max sequence length ===
MAX_LEN=$(seqkit fx2tab "$FASTA" | awk -F'\t' '{ print length($2) }' | sort -nr | head -n 1)
echo "Max sequence length: $MAX_LEN nt"



# If the max aligned length is < PV_DIM, run 3SEQ once on full FASTA
if (( MAX_LEN < PV_DIM )); then
    run_3seq_on_fasta "$FASTA" "$OUTDIR" "$BASE"
    exit 0
else
    WINDOW_DIR="${OUTDIR}/3seq_windows_${BASE}"
    mkdir -p ${WINDOW_DIR}

    python3 "${SCRIPT_DIR}/sliding_window.py" -s "$STEP" -W "$PV_DIM" "$WINDOW_DIR/window_{start}-{end}.fasta" "$FASTA"

    # === Run 3SEQ on each window ===
    for win_fasta in "$WINDOW_DIR"/*.fasta; do
        base_win=$(basename "$win_fasta" .fasta)

        num_seqs=$(grep -c '^>' "$win_fasta")
        if (( num_seqs < 3 )); then
            echo "Skipping $base_win — only $num_seqs sequences (need ≥ 3)"
            continue
        fi


        echo "Running 3SEQ on window: $base_win"
        run_3seq_on_fasta "$win_fasta" "$OUTDIR" "${BASE}_${base_win}"
    done
fi