#!/bin/bash
# run_tree_pipeline.sh: Download sequences, merge FASTA files, align codons, and convert to Phylip format.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GLOBALS="$1"
source "$GLOBALS"

## [Constants] ##
OUTPUT_INTR_DIR="${PROCESSED_DIR}/trees_output"
ERRORS_FILE="$OUTPUT_INTR_DIR/tree_pipeline_logs.txt"

IFS=',' read -r -a REMOVE_SAMPLES <<< "$2"
TREE_INTERMIDIATE="$OUTPUT_INTR_DIR/tree_tmp.tree"

mkdir -p "$OUTPUT_INTR_DIR"

cat "$ROOTED_TREE_PATH" > "$TREE_INTERMIDIATE"

echo "=== Step 1: Cleaning tree labels ==="
EXT="${TREE_INTERMIDIATE##*.}"
BASE="${TREE_INTERMIDIATE%.*}"
sed -E 's/\|[^():,;]*//g' "${TREE_INTERMIDIATE}" > "${BASE}_clean.${EXT}"
TREE_INTERMIDIATE="${BASE}_clean.${EXT}"
BASE="${BASE}_clean"

if [[ -z "${TARGET_LABEL:-}" ]]; then
  TARGET_LABEL=$(python3 "${SCRIPT_DIR}/find_latest_diverged_group.py" "$TREE_INTERMIDIATE" "$ACCESSIONS_FILE")
  echo "Auto-selected TARGET_LABEL=$TARGET_LABEL"
  if [[ -n "$TARGET_LABEL" ]]; then
    sed -i "s/^TARGET_LABEL=.*/TARGET_LABEL=\"$TARGET_LABEL\"/" "$GLOBALS"
  fi
fi

if [[ -n "${TARGET_LABEL:-}" ]]; then
  echo "=== Step 2: Mark tree ==="
  python3 "${SCRIPT_DIR}/mark_foreground.py" "$TREE_INTERMIDIATE" "$ACCESSIONS_FILE" "$TARGET_LABEL" > "${BASE}_marked.${EXT}"
  TREE_INTERMIDIATE="${BASE}_marked.${EXT}"
  BASE="${TREE_INTERMIDIATE%.*}"
  echo "Marked tree saved as $TREE_INTERMIDIATE"
else
  echo "Skipping marking and reducing tree (step 2)"
fi
echo "$TARGET_LABEL"

echo "=== Step 3: Creating leaf only tree (not used) ==="
sed -E '
  s/\)([^():,;#]+)(#[0-9]+)?(:[0-9.eE+-]+)?/\)\2\3/g;
  s/\)([^():,;#]+)(#[0-9]+)?;/)\2;/g
' "$TREE_INTERMIDIATE" > "${BASE}_leaf_only.${EXT}"
echo "leaf only tree saved in ${BASE}_leaf_only.${EXT}"

echo "=== Step 4: Removing leaves that are not found in ACCESSIONS.txt ==="
if [ ${#REMOVE_SAMPLES[@]} -ne 0 ]; then
    python3 "${SCRIPT_DIR}/prune_leaves_by_name.py" "$TREE_INTERMIDIATE" "${REMOVE_SAMPLES[@]}"
    TREE_INTERMIDIATE="${BASE}_pruned.${EXT}"
    BASE="${BASE}_pruned"
fi

echo "=== Step 5: Limiting tree to ${MAX_TREE_LEAVES} leaves ==="
python3 "$SCRIPT_DIR/prune_random_leaves.py" "$TREE_INTERMIDIATE" "${BASE}_limited.${EXT}" --max-leaves "$MAX_TREE_LEAVES"
TREE_INTERMIDIATE="${BASE}_limited.${EXT}"
BASE="${TREE_INTERMIDIATE%.*}"

mkdir -p "$(dirname "$FINAL_TREE_FILE_PATH")"
cp "$TREE_INTERMIDIATE" "$FINAL_TREE_FILE_PATH"
echo "Final tree saved as $FINAL_TREE_FILE_PATH"

echo "=== Pipeline complete ==="
