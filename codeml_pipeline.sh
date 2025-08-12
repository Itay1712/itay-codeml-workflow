#!/bin/bash

# f=codeml_pipeline.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && ./"$f" ./input/accessions.txt "Hepeviridae_10"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

f=$SCRIPT_DIR/scripts/check_deps.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && "$f" "python>=3.10 python<3.11 biopython==1.85 numpy scipy mafft gawk 3seq seqkit blast entrez-direct ete3 paml lxml pyqt pip packaging setuptools wheel TreeCluster"

CONFIG="$1"
GLOBALS="$SCRIPT_DIR/scripts/globals.sh"
f="$SCRIPT_DIR/scripts/generate_globals.sh" && sed -i 's/\r$//' "$f" && chmod +x "$f" && "$f" "$CONFIG" "$GLOBALS"
source "$GLOBALS"

source  ~/miniconda3/etc/profile.d/conda.sh
conda activate codeml_env

for f in \
  "$SCRIPT_DIR/scripts/prepare_codeml_input.sh"
do
  sed -i 's/\r$//' "$f"
  chmod +x "$f"
done

# ACCESSIONS_FILE="$1"
# GROUP="$2" # Hepeviridae_6
# ANALYSIS="site-model"  # "branch-site" or "site-model"


echo "*************************** GETTING CODEML INPUT ***************************"
echo $CODEML_GET_INPUT
echo $CODEML_RESULTS_DIR
# JF915746.1 - Hepeviridae_6
if [ "$CODEML_GET_INPUT" = "true" ]; then
  f=$SCRIPT_DIR/scripts/prepare_codeml_input.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && "$f" "$GLOBALS"
  ## reload globals in case REF_ACC or related paths were updated during input preparation
  source "$GLOBALS"
fi

# place inside $SCRIPT_DIR/codeml/input directory:
# ${GROUP}_${REF}_${CDS}.phy
# ${GROUP}_leaf_only_reduced.tree

echo "*************************** RUNNING CODEML ***************************"
RESULTS_DIR="$CODEML_RESULTS_DIR"
if ! grep -qE '^[0-9]+$' "$FINAL_TREE_FILE_PATH"; then
    sed -i '1i1' "$FINAL_TREE_FILE_PATH"
fi

if [ "$CODEML_RUN" = "true" ]; then
  CTL_TEMPLATE="$CODEML_INPUT_DIR/codeml_template.ctl"
  mkdir -p "$RESULTS_DIR"

  if [[ "$ANALYSIS" == "site-model" ]]; then
    # Loop through each .phy file in codeml/input
    for PHY_FILE in "$CODEML_INPUT_DIR"/${GROUP}_*.phy; do
      BASENAME="$(basename "$PHY_FILE" .phy)"
      echo "Running codeml on $BASENAME..."
      for model in M1a M2a; do
        NSsites=1; [[ "$model" == "M2a" ]] && NSsites=2
        OUT_DIR="$RESULTS_DIR/$BASENAME/$model"
        MLC_FILE="$OUT_DIR/mlc"

        if [[ -f "$MLC_FILE" ]]; then
          echo "Skipping $model: $MLC_FILE already exists."
          continue
        fi

        mkdir -p "$OUT_DIR"

        # Prepare a unique ctl file
        CTL_FILE="$OUT_DIR/codeml.ctl"
        cp "$CTL_TEMPLATE" "$CTL_FILE"

        # Create short symlinks for codeml to avoid long paths
        ln -sf "$PHY_FILE" "$OUT_DIR/aln.phy"
        ln -sf "$FINAL_TREE_FILE_PATH" "$OUT_DIR/tree.tre"
        ln -sf "$OUT_DIR/mlc" "$OUT_DIR/mlc_link"

        # Replace variables in ctl file using symlink paths
        sed -i "s|^[[:space:]]*seqfile.*|seqfile = aln.phy|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*treefile.*|treefile = tree.tre|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*outfile.*|outfile = mlc_link|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*model.*|model = 0|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*fix_omega.*|fix_omega = 0|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*omega.*|omega = 1|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*NSsites.*|NSsites = $NSsites|" "$CTL_FILE"

        # Run codeml
        cd "$OUT_DIR"
        codeml codeml.ctl
        rm -f aln.phy tree.tre mlc_link
        cd - > /dev/null
      done
    done
  fi

  if [[ "$ANALYSIS" == "branch-site" ]]; then
    for PHY_FILE in "$CODEML_INPUT_DIR"/${GROUP}_*.phy; do
      BASENAME="$(basename "$PHY_FILE" .phy)"
      echo "Running branch-site models on $BASENAME..."

      for model in Null Positive; do
        OUT_DIR="$RESULTS_DIR/$BASENAME/$model"
        MLC_FILE="$OUT_DIR/mlc"

        if [[ -f "$MLC_FILE" ]]; then
          echo "Skipping $model: $MLC_FILE already exists."
          continue
        fi

        mkdir -p "$OUT_DIR"

        CTL_FILE="$OUT_DIR/codeml.ctl"
        cp "$CTL_TEMPLATE" "$CTL_FILE"

        # Create short symlinks for codeml to avoid long paths
        ln -sf "$PHY_FILE" "$OUT_DIR/aln.phy"
        ln -sf "$FINAL_TREE_FILE_PATH" "$OUT_DIR/tree.tre"
        ln -sf "$OUT_DIR/mlc" "$OUT_DIR/mlc_link"

        # Shared settings for branch-site models using symlink paths
        sed -i "s|^[[:space:]]*seqfile.*|seqfile = aln.phy|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*treefile.*|treefile = tree.tre|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*outfile.*|outfile = mlc_link|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*model.*|model = 2|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*NSsites.*|NSsites = 2|" "$CTL_FILE"
        sed -i "s|^[[:space:]]*clock.*|clock = 0|" "$CTL_FILE"

        # Differentiate between null and positive model
        if [[ "$model" == "Null" ]]; then
            sed -i "s|^[[:space:]]*fix_omega.*|fix_omega = 1|" "$CTL_FILE"
            sed -i "s|^[[:space:]]*omega.*|omega = 1|" "$CTL_FILE"
        else
            sed -i "s|^[[:space:]]*fix_omega.*|fix_omega = 0|" "$CTL_FILE"
            sed -i "s|^[[:space:]]*omega.*|omega = 1|" "$CTL_FILE"
        fi

        cd "$OUT_DIR"
        codeml codeml.ctl
        rm -f aln.phy tree.tre mlc_link
        cd - > /dev/null
      done
    done
  fi
fi

echo "*************************** ANALYZING RESULTS ***************************"

if [ "$CODEML_ANALYSIS" = "true" ]; then
  SUMMARY_FILE="$RESULTS_DIR/summary_${ANALYSIS}_${GROUP}.tsv"
  echo -e "CDS\tLRT\tp-value\tSelectedSites" > "$SUMMARY_FILE"

  # ---------------------- SITE MODEL ANALYSIS ----------------------
  if [[ "$ANALYSIS" == "site-model" ]]; then
    for CDS_DIR in "$RESULTS_DIR"/${GROUP}_*/; do
      CDS_NAME="$(basename "$CDS_DIR")"
      M1A_MLC="${CDS_DIR}/M1a/mlc"
      M2A_MLC="${CDS_DIR}/M2a/mlc"

      if [[ ! -f "$M1A_MLC" || ! -f "$M2A_MLC" ]]; then
        echo "Skipping $CDS_NAME — missing M1a or M2a results"
        continue
      fi

      lnL1=$(grep -m1 "lnL" "$M1A_MLC" | awk '{print $(NF-1)}')
      lnL2=$(grep -m1 "lnL" "$M2A_MLC" | awk '{print $(NF-1)}')

      if [[ -z "$lnL1" || -z "$lnL2" ]]; then
        echo "Could not extract log-likelihoods for $CDS_NAME"
        continue
      fi

      LRT=$(echo "scale=5; 2 * ($lnL2 - $lnL1)" | bc)
      PVALUE=$(python3 -c "from scipy.stats import chi2; print(round(chi2.sf($LRT, 2), 6))")

      echo "$lnL1"
      echo "$lnL2"
      echo "$LRT"
      echo "$PVALUE"

      echo -e "$CDS_NAME\t$LRT\t$PVALUE\t" >> "$SUMMARY_FILE"
    done
  fi

  # ------------------- BRANCH-SITE MODEL ANALYSIS -------------------
  if [[ "$ANALYSIS" == "branch-site" ]]; then
    for CDS_DIR in "$RESULTS_DIR"/${GROUP}_*/; do
      CDS_NAME="$(basename "$CDS_DIR")"
      NULL_MLC="${CDS_DIR}/Null/mlc"
      POS_MLC="${CDS_DIR}/Positive/mlc"

      if [[ ! -f "$NULL_MLC" || ! -f "$POS_MLC" ]]; then
        echo "Skipping $CDS_NAME — missing Null or Positive results"
        continue
      fi

      lnL1=$(grep -m1 "lnL" "$NULL_MLC" | awk '{print $(NF-1)}')
      lnL2=$(grep -m1 "lnL" "$POS_MLC" | awk '{print $(NF-1)}')

      if [[ -z "$lnL1" || -z "$lnL2" ]]; then
        echo "Could not extract log-likelihoods for $CDS_NAME"
        continue
      fi

      LRT=$(echo "scale=5; 2 * ($lnL2 - $lnL1)" | bc)
      PVALUE=$(python3 -c "from scipy.stats import chi2; print(round(chi2.sf($LRT, 1), 6))")

      SELECTED=""
      if (( $(echo "$PVALUE < 0.05" | bc -l) )); then
        SELECTED=$(awk '/Bayes Empirical Bayes/,/^$/ {
          if ($0 ~ /^[ ]*[0-9]+/ && $(NF) ~ /\*/) {
            pos=$1; aa=$2; prob=$(NF-1)
            if (prob > 0.8) {
              printf "%s(%s,%.3f);", pos, aa, prob
            }
          }
        }' "$POS_MLC")
      fi

      echo -e "$CDS_NAME\t$LRT\t$PVALUE\t$SELECTED" >> "$SUMMARY_FILE"
    done
  fi
  echo "Summary saved to $SUMMARY_FILE"
fi


echo "*************************** CODEML PIPELINE COMPLETE ***************************"
