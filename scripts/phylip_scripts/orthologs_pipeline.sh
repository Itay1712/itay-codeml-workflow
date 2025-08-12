#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ACCESSIONS_FILE="$1"
# REF_ACC="$2"
# REF_CDS_ID="$3"
# PHY_FILE_TEMPLATE="$4"
# TREE_FILE="$5"
# GROUP="$6"
# RESULTS_DIR="${SCRIPT_DIR}/../../pipeline_output"

GLOBALS="$1"
source "$GLOBALS"

for f in \
  "$SCRIPT_DIR/find_top_virus.sh" \
  "$SCRIPT_DIR/find_orthologs.sh"
do
  sed -i 's/\r$//' "$f"
  chmod +x "$f"
done


print_box() {
  local lines=()
  local maxlen=0

  # Read input lines
  while IFS= read -r line; do
    lines+=("$line")
    (( ${#line} > maxlen )) && maxlen=${#line}
  done

  local border="+-$(printf '%*s' "$maxlen" '' | tr ' ' '-')-+"
  echo "$border"
  for line in "${lines[@]}"; do
    printf "| %-*s |\n" "$maxlen" "$line"
  done
  echo "$border"
}

generate_unique_filename() {
    local base="$1"
    local ext="$2"
    local filepath="${base}.${ext}"

    if [[ ! -e "$filepath" ]]; then
        echo "$filepath"
        return
    fi

    local counter=1
    while [[ -e "${base} (${counter}).${ext}" ]]; do
        ((counter++))
    done
    echo "${base} (${counter}).${ext}"
}

HEADER_LINE=$(head -n 1 "$ACCESSIONS_FILE")
IFS=$'\t' read -ra HEADERS <<< "$HEADER_LINE"

TREE_FILE="$FINAL_TREE_FILE_TEMPLATE"

ACCESSION_COL_INDEX=-1
HOST_COL_INDEX=-1

for i in "${!HEADERS[@]}"; do
    case "${HEADERS[$i]}" in
        "ACCESSION"|"Accession")
            ACCESSION_COL_INDEX="$i"
            ;;
        "HOST"|"Host")
            HOST_COL_INDEX="$i"
            ;;
    esac
done

if [[ -z "$ACCESSION_COL_INDEX" || -z "$HOST_COL_INDEX" ]]; then
    echo "ERROR: Unable to locate ACCESSION or HOST column in accessions.txt"
    exit 1
fi

echo "[Stage 1] Getting reference and targets..."
if [[ -z "$REF_ACC" ]]; then
	REF_ACC=$(bash $SCRIPT_DIR/find_top_virus.sh "$ACCESSIONS_FILE" "$ACCESSION_COL_INDEX")
fi

TARGETS=($(tail -n +2 "$ACCESSIONS_FILE" | awk -F'\t' -v col=$((ACCESSION_COL_INDEX + 1)) -v ref="$REF_ACC" '$col != ref {print $col}'))
{
  echo "Ref:"
  echo "- $REF_ACC"
  echo "Targets:"
  for target in "${TARGETS[@]}"; do
    echo "- $target"
  done
} | print_box
ORTHOLOGS_RESULTS_DIR="${RESULTS_DIR}"

echo "[Stage 2] Finding orthologs..."

## If ACCESSION_FILE is test.txt, we process the target regardless of the file's existence.
## If ACCESSION_FILE is not test.txt, we only process the target if the file does not exist.
filtered_targets=()
for tgt in "${TARGETS[@]}"; do
  if [[ "$(basename "${ACCESSIONS_FILE:-}")" == "test.txt" || ! -f "${ORTHOLOGS_RESULTS_DIR}/reciprocal_pairs/${tgt}_reciprocal_pairs.tsv" ]]; then
    filtered_targets+=("$tgt")
  fi
done
TARGETS=("${filtered_targets[@]}")
${SCRIPT_DIR}/find_orthologs.sh "$REF_ACC" "$REF_CDS_ID" "$ORTHOLOGS_RESULTS_DIR" "$PREFETCH_DIR" "${TARGETS[@]}"

echo "[Stage 3] Processing each ref CDS..."

for global in "${ORTHOLOGS_RESULTS_DIR}/globals/"*.fasta; do
  awk '/^>/ {sub(/\|.*/, "", $0)} {print}' "$global" > tmp && mv tmp "$global"

  base="$(basename "$global" .fasta)"
  workg="${ORTHOLOGS_RESULTS_DIR}/globals/${base}"
  mkdir -p "$workg"
  IFS="_" read -r param1 CDS <<< "$base"
  echo -e "\n\n\n\n========================================== Processing CDS: [${base}] ==========================================\n\n\n\n"

  echo -e "\n\n--------------[3.1] Codon alignment\n\n"
  if [ "$ALIGN_CODONS_WITH" = "prank" ]; then
    time prank -d="$global" -o="${workg}/aligned_temp" -t="$FINAL_TREE_FILE_TEMPLATE" -once -f=fasta +F -codon
    mv "${workg}/aligned_temp.best.fas" "${workg}/aligned_codons.fasta"
  else
    python3 "$SCRIPT_DIR/align_codons.py" "$global" "${workg}/aligned_codons.fasta"
  fi
  echo "Initial aligned FASTA file created at: ${workg}/aligned_codons.fasta"


  echo -e "\n\n--------------[3.2] Masking poorly aligned regions\n\n"
  python3 ${SCRIPT_DIR}/mask_alignment.py -i "${workg}/aligned_codons.fasta" \
	--aa-threshold 0.85 \
	--blosum62-threshold 0.2 \
	--o "${workg}/aligned_codons_masked_poor.fasta"

  echo -e "\n\n--------------[3.3] Recombination filtering report\n\n"

  FASTA_INPUT="${workg}/aligned_codons_masked_poor.fasta"
  RECOMB_OUTPUT_DIR="${RECOMB_OUTPUT_DIR//<CDS>/$CDS}"
  BASE_NAME="${base}"
  f=$SCRIPT_DIR/recombination_filtering.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && "$f" "$GLOBALS" "$FASTA_INPUT" "$RECOMB_OUTPUT_DIR" "$BASE_NAME"


  echo -e "\n\n--------------[3.4] Masking recombination regions\n\n"

  MASK_FILE="${RECOMB_OUT}/recombination_regions.mask.tsv"
  if [[ -f "$MASK_FILE" ]]; then
    python3 ${SCRIPT_DIR}/mask_recomb_regions.py \
      "${workg}/aligned_codons_masked_poor.fasta" \
      "$MASK_FILE" \
      "${workg}/aligned_codons_masked.fasta" \
      --mask-char N
  else
      echo "No mask file for $BASE — skipping masking recom regions"
      cp "${workg}/aligned_codons_masked_poor.fasta" "${workg}/aligned_codons_masked.fasta"
  fi



  echo -e "\n\n--------------[3.5] Convert to Phylip\n\n"
 
  mkdir -p "$(dirname "$PHY_FILE_TEMPLATE")"
  PHY_OUT="${PHY_FILE_TEMPLATE//<CDS>/$CDS}"
  # PHY_EXT="${PHY_FILE_TEMPLATE##*.}"
  # PHY_BASE="${PHY_FILE_TEMPLATE%.*}"
  # PHY_OUT="${PHY_BASE}_${base}.${PHY_EXT}"
  touch $PHY_OUT
  python3 "${SCRIPT_DIR}/fasta_to_phylip.py" "${workg}/aligned_codons_masked.fasta" "$PHY_OUT"
  echo "Phylip file created at: $PHY_OUT"

echo "=== Pipeline complete ==="
done



rm -rf tmpdirprankmsa*
