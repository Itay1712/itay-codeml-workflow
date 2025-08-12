#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="$1"
GLOBALS_FILE="${2:-$(dirname "$CONFIG_FILE")/globals.sh}"

CONFIG_DIR="$(cd "$(dirname "$CONFIG_FILE")" && pwd)"
SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse YAML into key\tvalue lines using python
parse_yaml() {
  python3 - "$CONFIG_FILE" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    data = yaml.safe_load(f)
for k,v in data.items():
    if isinstance(v, bool):
        print(f"{k}\t{str(v).lower()}")
    else:
        print(f"{k}\t{v}")
PY
}

declare -A VARS
# ---- Config variables ----
while IFS=$'\t' read -r key value; do
  VARS[$key]="$value"
done < <(parse_yaml)

VARS[CONFIG_DIR]="$CONFIG_DIR"

# === Resolve OUTPUT_DIR ===
BASE_OUTPUT="${VARS[OUTPUT_DIR]:-$CONFIG_DIR}"
if [[ "$BASE_OUTPUT" != /* ]]; then
  BASE_OUTPUT="${CONFIG_DIR}/${BASE_OUTPUT}"
fi
mkdir -p "$BASE_OUTPUT"
VARS[OUTPUT_DIR]="$BASE_OUTPUT"
VARS[PROCESSED_DIR]="${BASE_OUTPUT}/processed"
VARS[CODEML_DIR]="${BASE_OUTPUT}/codeml"
mkdir -p "${VARS[PROCESSED_DIR]}" "${VARS[CODEML_DIR]}/input" "${VARS[CODEML_DIR]}/output"

# Shared cache for downloaded sequences
VARS[PREFETCH_DIR]="${BASE_OUTPUT}/prefetch"
mkdir -p "${VARS[PREFETCH_DIR]}"

# Ensure codeml template exists
TEMPLATE_DEST="${VARS[CODEML_DIR]}/input/codeml_template.ctl"
TEMPLATE_URL="https://raw.githubusercontent.com/abacus-gene/paml-tutorial/main/positive-selection/templates/template_CODEML.ctl"
if [[ ! -f "$TEMPLATE_DEST" ]]; then
  curl -L -o "$TEMPLATE_DEST" "$TEMPLATE_URL"
fi

# === Validate ACCESSIONS_FILE ===
if [[ -z "${VARS[ACCESSIONS_FILE]:-}" || "${VARS[ACCESSIONS_FILE]}" == "-" ]]; then
  VARS[ACCESSIONS_FILE]="${VARS[PROCESSED_DIR]}/accessions_${VARS[GROUP]}.txt"
fi
if [[ "${VARS[ACCESSIONS_FILE]}" != /* ]]; then
  VARS[ACCESSIONS_FILE]="${CONFIG_DIR}/${VARS[ACCESSIONS_FILE]}"
fi
if [[ ! -f "${VARS[ACCESSIONS_FILE]}" ]]; then
  echo "<-> ACCESSIONS_FILE not found — auto-generating from GROUP=${VARS[GROUP]}"
  mkdir -p "$(dirname "${VARS[ACCESSIONS_FILE]}")"
  : > "${VARS[ACCESSIONS_FILE]}"
  python3 "$SCRIPT_ROOT/../generate_accessions/get_accession_file.py" "${VARS[GROUP]}" "${VARS[ACCESSIONS_FILE]}"
fi

# === Resolve input files ===
if [[ -n "${VARS[ROOTED_TREE_PATH]:-}" && "${VARS[ROOTED_TREE_PATH]}" != /* ]]; then
  VARS[ROOTED_TREE_PATH]="${CONFIG_DIR}/${VARS[ROOTED_TREE_PATH]}"
fi

# === Check REF_ACC (selection deferred) ===
# The reference accession must correspond to a sample that remains after
# pruning the tree and filtering the accession list.  Selecting it now may
# pick a sequence that will later be removed.  We therefore postpone the
# automatic selection until after sample pruning inside `prepare_codeml_input.sh`.
# If a REF_ACC is provided in the config we keep it; otherwise leave it empty
# so that it can be determined later.
if [[ -z "${VARS[REF_ACC]:-}" || "${VARS[REF_ACC]}" == "-" ]]; then
  VARS[REF_ACC]=""
  echo "<-> REF_ACC will be selected after sample pruning"
fi

# === Check PV_DIM ===
if (( VARS[PV_DIM] > 400 )); then
  echo "!!! PV_DIM too large (${VARS[PV_DIM]}) — resetting to 400"
  VARS[PV_DIM]=400
fi

# === Default MAX_TREE_LEAVES ===
VARS[MAX_TREE_LEAVES]="${VARS[MAX_TREE_LEAVES]:-150}"

# === Validate ANALYSIS ===
if [[ "${VARS[ANALYSIS]}" != "site-model" && "${VARS[ANALYSIS]}" != "branch-site" ]]; then
  echo "!!! ANALYSIS '${VARS[ANALYSIS]}' is invalid — resetting to 'site-model'"
  VARS[ANALYSIS]="site-model"
fi

# ---- Default variables ----
VARS[RESULTS_DIR]="${VARS[PROCESSED_DIR]}/results_${VARS[REF_ACC]}_${VARS[GROUP]}"
VARS[MASKING_OUTPUT_DIR]="${VARS[RESULTS_DIR]}/masked_alignments"
if [[ "${VARS[PV_TABLE_FILE]}" != /* ]]; then
  VARS[PV_TABLE_FILE]="${CONFIG_DIR}/${VARS[PV_TABLE_FILE]}"
fi
VARS[RECOMB_OUTPUT_DIR]="${VARS[RESULTS_DIR]}/3seq_report_${VARS[REF_ACC]}_<CDS>"
VARS[FINAL_TREE_FILE_PATH]="${VARS[CODEML_DIR]}/input/${VARS[GROUP]}.tree"
VARS[FINAL_TREE_FILE_TEMPLATE]="${VARS[FINAL_TREE_FILE_PATH]}"
VARS[PHY_FILE_TEMPLATE]="${VARS[CODEML_DIR]}/input/${VARS[GROUP]}_${VARS[REF_ACC]}_<CDS>.phy"
VARS[CODEML_INPUT_DIR]="${VARS[CODEML_DIR]}/input"
VARS[CODEML_RESULTS_DIR]="${VARS[CODEML_RESULTS_DIR]:-${VARS[CODEML_DIR]}/output/${VARS[GROUP]}}" ## allow config overwrite with :-
VARS[LOG_FILE]="${VARS[OUTPUT_DIR]}/log.txt"

# Expand ${VAR} patterns
# helper to expand ${VAR} patterns
# export variables for envsubst
for key in "${!VARS[@]}"; do
  export "$key"="${VARS[$key]}"
done

for key in "${!VARS[@]}"; do
  VARS[$key]="$(echo "${VARS[$key]}" | envsubst)" # Handles variable and command substitution
done

# Write globals file
{
  echo "# Generated globals"
  for key in "${!VARS[@]}"; do
    echo "${key}=\"${VARS[$key]}\""
  done
} > "$GLOBALS_FILE"

echo "Globals written to $GLOBALS_FILE"
