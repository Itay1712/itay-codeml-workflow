#!/usr/bin/env bash
set -euo pipefail


# process_target.sh
# f=find_orthologs.sh && sed -i 's/\r$//' "$f" && chmod +x "$f" && ./"$f" AB740222.1 "" AB740220.1
# Usage:
#   ./process_targets_multi.sh REF_ACCESSION [REF_PROTEIN_ID] TARGET1 [TARGET2 ...]
# Example:
#   ./process_targets_multi.sh AB740220.1 BAM68890.1 JX565469.1 AB740222.1

if (( $# < 3 )); then
  echo "Usage: $0 REF_ACC REF_CDS_ID TARGET1 [TARGET2 ...]"
  exit 1
fi

# --- arguments & setup ---
ref_acc="$1"                 # e.g. AB740220.1
ref_cds_id="$2"              # e.g. BAM68890.1
refdir="$3"            ## for output
PREFETCH_DIR="$4"      ## data
shift 4
targets=("$@")               # remaining args

: "${PREFETCH_DIR:?PREFETCH_DIR not set}"
mkdir -p "$PREFETCH_DIR"

if (( ${#targets[@]} == 0 )); then
  echo "find_orthologs.sh | Error: need at least one TARGET accession"
  exit 1
fi

# PARAMETERS
[[ -n "$ref_cds_id" ]] && refdir+="_${ref_cds_id}"
mkdir -p "$refdir"
mkdir -p "${refdir}/globals"

MAX_JOBS=4     # adjust to your cores
EVALUE=1e-3
PIDENT=30

### 1) PREPARE REFERENCE & GLOBAL FASTAs ###

if [[ -n "$ref_cds_id" ]]; then
  echo "find_orthologs.sh | [REF] preparing full proteome & CDS…"
else
  echo "find_orthologs.sh | [REF] extracting only $ref_cds_id from $ref_acc…"
fi

# 1) Retrieve the raw AA FASTA so we can grab the real header
raw_aa="${refdir}/ref_proteins_raw.fasta"
prefetch_aa="$PREFETCH_DIR/${ref_acc}_cds_aa.fasta"
if [[ ! -s "$prefetch_aa" ]]; then
  efetch -db nucleotide -id "$ref_acc" -format fasta_cds_aa > "$prefetch_aa"
fi
cp "$prefetch_aa" "$raw_aa"

headers_list="${refdir}/headers_list.txt"  # or whatever file you want

# prepare/empty the header list
: > "$headers_list"

# extract headers based on ref_cds_id
if [[ -z "$ref_cds_id" ]]; then
    # if ref_cds_id is empty, keep all headers
    grep '^>' "$raw_aa" >> "$headers_list"
else
    # if ref_cds_id is given, keep only matching header
    grep -E "^>.*$ref_cds_id" "$raw_aa" >> "$headers_list"
fi

# 3) Now simplify the headers *from that same raw file*
sed -E 's/^>.*\[protein_id=([^]]+)\].*/>'"$ref_acc"'|\1/' \
"$raw_aa" \
> "${refdir}/ref_proteins.fasta"
rm -f "$raw_aa"

if [[ -n "$ref_cds_id" ]]; then
	# 4) Extract just the one block (with our new simple header) 
	awk -v id=">${ref_acc}|${ref_cds_id}" '
	  $0 == id {print; found=1; next}
	  found && /^>/ {exit}
	  found {print}
	' "${refdir}/ref_proteins.fasta" > "${refdir}/ref_proteins.fasta.tmp" \
	&& mv "${refdir}/ref_proteins.fasta.tmp" "${refdir}/ref_proteins.fasta"
fi

# repeat *exactly* the same pattern for the nucleotide CDS:
raw_na="${refdir}/ref_cds_raw.fasta"
prefetch_na="$PREFETCH_DIR/${ref_acc}_cds_na.fasta"
if [[ ! -s "$prefetch_na" ]]; then
  efetch -db nucleotide -id "$ref_acc" -format fasta_cds_na > "$prefetch_na"
fi
cp "$prefetch_na" "$raw_na"

# (we already have ORIGINAL_HEADER, so no need to recalc)
sed -E 's/^>.*\[protein_id=([^]]+)\].*/>'"$ref_acc"'|\1/' \
"$raw_na" \
> "${refdir}/ref_cds.fasta"
rm -f "$raw_na"

if [[ -n "$ref_cds_id" ]]; then
	awk -v id=">${ref_acc}|${ref_cds_id}" '
	  $0 == id {print; found=1; next}
	  found && /^>/ {exit}
	  found {print}
	' "${refdir}/ref_cds.fasta" > "${refdir}/ref_cds.tmp" && \
	mv "${refdir}/ref_cds.tmp" "${refdir}/ref_cds.fasta"
fi

echo "find_orthologs.sh | [REF] building reference BLAST DB…"
makeblastdb -in "${refdir}/ref_proteins.fasta" -dbtype prot \
            -out "${refdir}/ref_prot_db"

echo "find_orthologs.sh | [REF] preparing global FASTAs for each ref CDS…"
# iterate over each header in ref_cds.fasta:
grep '^>' "${refdir}/ref_cds.fasta" | sed 's/^>//' \
| while IFS= read -r qid; do
  # sanitize filename
  file="${refdir}/globals/${qid//|/_}.fasta"
  # empty or overwrite
  : > "$file"
  # append the reference sequence block
  # sed -e '$a\' "${refdir}/ref_cds.fasta" | \
  # awk -v id=">${qid}" '
    # $0 == id {print; found=1; next}
    # found && /^>/ {exit}
    # found {print}
  # ' >> "$file"
  
	{
	  printf "\n"   # make sure we start on a fresh line
	  awk -v id="$qid" '
		BEGIN {
		  RS = ">"      # read whole FASTA records at a time (splits on ">")
                  ORS = ""      # do not auto-reinsert RS
		}
		NR > 1 {        # skip the very first empty chunk before the first ">"
		  # extract the header (up to first newline)
		  split($0, rec, "\n")
		  if (rec[1] == id) {
			# print ">" plus the entire record (header + all its seq lines)
			print ">" $0
		  }
		}
	  ' "${refdir}/ref_cds.fasta"
	} >> "$file"
done

### 2) DOWNLOAD ALL TARGETS IN PARALLEL ###
echo "find_orthologs.sh | [ALL] Stage 2: fetching all target FASTAs…"
for tgt in "${targets[@]}"; do
  (
    mkdir -p "${refdir}/reciprocal_pairs/${tgt}"

    prot_out="${refdir}/reciprocal_pairs/${tgt}/${tgt}_proteins.fasta"
    cds_out="${refdir}/reciprocal_pairs/${tgt}/${tgt}_cds.fasta"

    prot_cache="$PREFETCH_DIR/${tgt}_cds_aa.fasta"
    cds_cache="$PREFETCH_DIR/${tgt}_cds_na.fasta"

    [[ -s "$prot_cache" ]] || efetch -db nucleotide -id "$tgt" -format fasta_cds_aa > "$prot_cache"
    [[ -s "$cds_cache" ]] || efetch -db nucleotide -id "$tgt" -format fasta_cds_na > "$cds_cache"

    sed -E 's/^>.*\[protein_id=([^]]+)\].*/>'"$tgt"'|\1/' "$prot_cache" > "$prot_out"
    sed -E 's/^>.*\[protein_id=([^]]+)\].*/>'"$tgt"'|\1/' "$cds_cache" > "$cds_out"
  ) &
  while (( $(jobs -rp|wc -l) >= 1 )); do sleep 1; done
done
wait

### 3) BUILD BLAST DBs FOR ALL TARGETS IN PARALLEL ###
echo "find_orthologs.sh | [ALL] Stage 3: building target BLAST DBs…"
for tgt in "${targets[@]}"; do
  (
    makeblastdb -in "${refdir}/reciprocal_pairs/${tgt}/${tgt}_proteins.fasta" \
                -dbtype prot \
                -out "${refdir}/reciprocal_pairs/${tgt}/${tgt}_prot_db"
  ) &
  while (( $(jobs -rp|wc -l) >= MAX_JOBS )); do sleep 1; done
done
wait

### 4) PROCESS EACH TARGET, APPENDING TO GLOBAL FASTAs ###
echo "find_orthologs.sh | [ALL] Stage 4: reciprocal‐BLAST + append to globals…"
for tgt in "${targets[@]}"; do
  (
    work="${refdir}/reciprocal_pairs"

    # --- Reciprocal BLASTp (as before) ---
    blastp -num_threads 1 \
      -query "${work}/${tgt}/${tgt}_proteins.fasta" \
      -db "${refdir}/ref_prot_db" \
      -evalue "$EVALUE" \
      -outfmt "6 qseqid sseqid pident evalue bitscore" \
      > "${work}/${tgt}/ref_vs_${tgt}.blast"

    blastp -num_threads 1 \
      -query "${refdir}/ref_proteins.fasta" \
      -db "${work}/${tgt}/${tgt}_prot_db" \
      -evalue "$EVALUE" \
      -outfmt "6 qseqid sseqid pident evalue bitscore" \
      > "${work}/${tgt}/${tgt}_vs_ref.blast"

    # --- pick best and find true 1:1 reciprocals ---
    sort -k1,1 -k4,4g -k5,5nr "${work}/${tgt}/ref_vs_${tgt}.blast" \
      | awk '!h1[$1]++ && $3>= '"$PIDENT"'' > "${work}/${tgt}/best1.tsv"
    sort -k1,1 -k4,4g -k5,5nr "${work}/${tgt}/${tgt}_vs_ref.blast" \
      | awk '!h2[$1]++ && $3>= '"$PIDENT"'' > "${work}/${tgt}/best2.tsv"

    awk '
      NR==FNR { seen[$2 "|" $1]; next }
      ($1 "|" $2) in seen { print $1 "\t" $2 }
    ' "${work}/${tgt}/best2.tsv" "${work}/${tgt}/best1.tsv" \
      > "${work}/${tgt}_reciprocal_pairs.tsv"

    # --- append each target CDS to the correct global FASTA ---

    while IFS=$'\t' read -r tgt_id ref_id; do
      global="${refdir}/globals/${ref_id//|/_}.fasta"
		{
		  awk -v id="$tgt_id" '
			BEGIN {
			  RS  = ">"     # read whole FASTA entries at once
                            ORS = ""      # suppress awk record-separator
			}
			NR > 1 {
			  # split record into lines: rec[1]=header, rec[2..n]=seq (or blanks)
			  n = split($0, rec, "\n")
			  # match the ID exactly (allowing optional desc after it)
			  if ( rec[1] == id ) {
				# print exactly one blank line, then header+newline
				printf "\n>%s\n", rec[1]
                                  # print only non-empty sequence lines
				for (i = 2; i <= n; i++) {
				  if (rec[i] ~ /[^[:space:]]/) {
					printf "%s\n", rec[i]
				  }
				}
			  }
			}
		  ' "${work}/${tgt}/${tgt}_cds.fasta"
		} >> "$global"
    done < "${work}/${tgt}_reciprocal_pairs.tsv"



    # while IFS=$'\t' read -r tgt_id ref_id; do
      # global="${refdir}/globals/${ref_id//|/_}.fasta"
      # # extract the target block from this target’s CDS.fasta
	  # sed -e '$a\' "${work}/${tgt}/${tgt}_cds.fasta" | \
	  # awk -v id=">${tgt_id}" '
	    # $0 == id {print; found=1; next}
	    # found && /^>/ {exit}
	    # found {print}
	  # ' >> "$global"

    # done < "${work}/${tgt}_reciprocal_pairs.tsv"

    # cleanup intermediates
        rm -r -- "${work:?}/${tgt}"
    #rm -f "${work}"/*.blast "${work}/best1.tsv" "${work}/best2.tsv"
    #rm -f "${work}/${tgt}_proteins.fasta" "${work}/${tgt}_cds.fasta"
    #rm -f "${work}/${tgt}_prot_db".*
  ) &
  # throttle parallel jobs
  while (( $(jobs -rp|wc -l) >= MAX_JOBS )); do sleep 1; done
done
wait

rm -f "${refdir}/ref_prot_db".*
