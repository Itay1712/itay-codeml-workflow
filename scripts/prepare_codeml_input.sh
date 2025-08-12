#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ACCESSIONS_FILE="$1"
# GROUP="$2" # Hepeviridae_6
# REF_ACC="${3:-}"
# ROOTED_TREE_PATH="${SCRIPT_DIR}/../input/rooted_trees/${GROUP}.rooted.anc_recon.tree"

GLOBALS="$1"
source "$GLOBALS"

# f=tree_scripts/run_tree_pipeline.sh && sed -i 's/\r$//' "$f" && chmod +x "$f"
# f=phylip_scripts/orthologs_pipeline.sh && sed -i 's/\r$//' "$f" && chmod +x "$f"
# f=prepare_codeml_input.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && ./"$f" ./accessions.txt "Hepeviridae_6"

for f in \
  "$SCRIPT_DIR/tree_scripts/run_tree_pipeline.sh"
do
  sed -i 's/\r$//' "$f"
  chmod +x "$f"
done

echo "########################################### PREPARING SAMPLES ###########################################"

f=$SCRIPT_DIR/remove_no_cds_samples.sh && \
  sed -i 's/\r$//' "$f" && \
  chmod +x "$f" && \
  "$f" "$ACCESSIONS_FILE" "$PREFETCH_DIR"

TREE_TMP=$(mktemp)
ACCESSIONS_TMP=$(mktemp)
MISSING_TMP=$(mktemp)
grep -oE '[^(),:;|]+' "$ROOTED_TREE_PATH" | grep -E '^[A-Z]{2}[0-9]+\.[0-9]+$' | sort -u > "$TREE_TMP"
#tail -n +2 "$ACCESSIONS_FILE" | sort -u > "$ACCESSIONS_TMP"
tail -n +2 "$ACCESSIONS_FILE" | cut -f2 | sort -u > "$ACCESSIONS_TMP"


comm -23 "$TREE_TMP" "$ACCESSIONS_TMP" > "$MISSING_TMP"


if [[ -s "$MISSING_TMP" ]]; then
  MISSING=$(paste -sd, "$MISSING_TMP")
  echo "Missing samples at $ROOTED_TREE_PATH:\n${MISSING}"
else
  echo "All tree accessions are present in ${ACCESSIONS_FILE}."
fi
MISSING=$(echo "$MISSING" | tr -d '\r')

rm -f "$TREE_TMP" "$ACCESSIONS_TMP" "$MISSING_TMP"



# ================ PIPELINE WORKFLOW =================
echo "########################################### TREE PIPELINE ###########################################"
f=$SCRIPT_DIR/tree_scripts/run_tree_pipeline.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && "$f" "$GLOBALS" "$MISSING"

HEADER_LINE=$(head -n 1 "$ACCESSIONS_FILE")
IFS=$'\t' read -ra HEADERS <<< "$HEADER_LINE"

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

# Step 1: Create set of accessions in reduced tree
mapfile -t KEPT_LEAVES < <(python3 -c "
from ete3 import Tree
t = Tree('$FINAL_TREE_FILE_PATH', format=1)
print('\n'.join(leaf.name for leaf in t.get_leaves()))
")
declare -A LEAF_SET
for acc in "${KEPT_LEAVES[@]}"; do
	cleaned_acc=$(echo "$acc" | sed 's/#\([1-9]\+\)//g')
	LEAF_SET["$cleaned_acc"]=1
done

# Print the number of leaves
echo "Number of leaves in the reduced tree: ${#LEAF_SET[@]}"

# Step 2: Filter ACCESSIONS_FILE by these accessions
HEADER=$(head -n 1 "$ACCESSIONS_FILE")
echo "$HEADER" > "${ACCESSIONS_FILE}.filtered"

tail -n +2 "$ACCESSIONS_FILE" | while IFS=$'\t' read -ra FIELDS; do
        ACCESSION="${FIELDS[$ACCESSION_COL_INDEX]}"
        echo "Comparing: '$ACCESSION' vs '${cleaned_acc}'"
        if [[ -n "${LEAF_SET[$ACCESSION]}" ]]; then
                (IFS=$'\t'; echo -e "${FIELDS[*]}") >> "${ACCESSIONS_FILE}.filtered"
        fi
done

#Replace original
mv "${ACCESSIONS_FILE}.filtered" "$ACCESSIONS_FILE"

# After pruning, determine a valid REF_ACC and update globals
if ! awk -v col="$((ACCESSION_COL_INDEX + 1))" -F'\t' 'NR > 1 { print $col }' "$ACCESSIONS_FILE" | grep -qxF "$REF_ACC"; then
    echo "<-> Selecting REF_ACC from filtered ACCESSIONS_FILE"
    REF_ACC=$(bash "$SCRIPT_DIR/phylip_scripts/find_top_virus.sh" "$ACCESSIONS_FILE" "$ACCESSION_COL_INDEX" "$PREFETCH_DIR")
fi

RESULTS_DIR="${PROCESSED_DIR}/results_${REF_ACC}_${GROUP}"
RECOMB_OUTPUT_DIR="${RESULTS_DIR}/3seq_report_${REF_ACC}_<CDS>"
PHY_FILE_TEMPLATE="${CODEML_DIR}/input/${GROUP}_${REF_ACC}_<CDS>.phy"

update_global() {
  local var="$1" value="$2"
  if grep -q "^${var}=" "$GLOBALS"; then
    sed -i "s|^${var}=.*|${var}=\"${value}\"|" "$GLOBALS"
  else
    echo "${var}=\"${value}\"" >> "$GLOBALS"
  fi
}

update_global REF_ACC "$REF_ACC"
update_global RESULTS_DIR "$RESULTS_DIR"
update_global RECOMB_OUTPUT_DIR "$RECOMB_OUTPUT_DIR"
update_global PHY_FILE_TEMPLATE "$PHY_FILE_TEMPLATE"

export REF_ACC RESULTS_DIR RECOMB_OUTPUT_DIR PHY_FILE_TEMPLATE

echo "########################################### PHYLIP PIPELINE ###########################################"

f=$SCRIPT_DIR/phylip_scripts/orthologs_pipeline.sh && \
 sed -i 's/\r$//' "$f" && \
 chmod +x "$f" && \
 "$f" "$GLOBALS"

