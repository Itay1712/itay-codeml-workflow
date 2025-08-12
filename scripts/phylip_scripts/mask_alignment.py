#!/usr/bin/env python3
"""
mark_poorly_aligned_regions.py

Process PRANK codon-aligned FASTA output, identify poorly aligned regions using one or more metrics:
  1. Simple codon identity
  2. Weighted nucleotide identity (user-defined weights per codon position)
  3. Amino acid substitution score (e.g., BLOSUM62)

Aggregates contiguous low-quality codons into regions and optionally outputs:
  - GFF3 annotation of poor regions
  - Masked FASTA with those codons replaced by 'NNN'

Usage examples:
  # Original: simple codon identity threshold
  python mark_poorly_aligned_regions.py -i aln.fasta --threshold 0.8

  # Weighted nucleotide identity (e.g., third position half-weight)
  python mark_poorly_aligned_regions.py -i aln.fasta --use-weighted-nuc --weights 1,1,0.5 --threshold 0.75

  # Amino acid substitution scoring (BLOSUM62) with threshold
  python mark_poorly_aligned_regions.py -i aln.fasta --use-aa --aa-threshold 0.6

  # Combine both metrics into composite score
  python mark_poorly_aligned_regions.py -i aln.fasta --use-weighted-nuc --use-aa --weights 1,1,0.5 --threshold 0.7

"""
import argparse
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

# Precompute BLOSUM62 score range
_BLOSUM_SCORES = list(blosum62.values())
MIN_BLOSUM, MAX_BLOSUM = min(_BLOSUM_SCORES), max(_BLOSUM_SCORES)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify poorly aligned codon regions in a PRANK codon alignment"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to codon-aligned FASTA file (PRANK output)"
    )
    parser.add_argument(
        "--codons-threshold", type=float,
        help="General threshold for chosen metric(s) (0-1)"
    )
    parser.add_argument(
        "--use-weighted-nuc", action="store_true",
        help="Use weighted nucleotide identity metric for codons threshold"
    )
    parser.add_argument(
        "--weights", default="1,1,0.5",
        help="Comma-separated weights for codon positions (first,second,third)"
    )
    parser.add_argument(
        "--aa-threshold", type=float,
        help="Threshold for consensus amino acid % (0-1), any aa above avoid masking"
    )
    parser.add_argument(
        "--blosum62-threshold", type=float,
        help="Threshold for normlized blosum62 score for consensus aa (0-1), any aa avoid masking"
    )
    parser.add_argument(
        "--o", default=None,
        help="Output FASTA where poorly aligned codons are masked as 'NNN'"
    )
    parser.add_argument(
        "--gff-out", default=None,
        help="Output GFF3 file listing poorly aligned regions"
    )
    return parser.parse_args()


def load_alignment(path):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError("No sequences found in input file.")
    length = len(records[0].seq)
    if any(len(r.seq) != length for r in records):
        raise ValueError("All sequences must be same length in alignment.")
    if length % 3 != 0:
        raise ValueError("Alignment length not a multiple of 3 (codon-aligned?).")
    return records, length


def compute_codon_identity(records, length):
    n = len(records)
    identities = []
    for start in range(0, length, 3):
        codons = [str(r.seq[start:start+3]) for r in records]
        most_common, count = Counter(codons).most_common(1)[0]
        identities.append(count / n)
    return identities


def compute_weighted_nuc_identity(records, length, weights):
    n = len(records)
    w = list(weights)
    total_w = sum(w)
    identities = []
    for start in range(0, length, 3):
        # consensus base at each position
        cols = [ [r.seq[start+p] for r in records] for p in range(3) ]
        consensus_base = [Counter(col).most_common(1)[0][0] for col in cols]
        # compute per-codon weighted identity then average
        codon_scores = []
        for r in records:
            cod = str(r.seq[start:start+3])
            score = 0.0
            for p, base in enumerate(cod):
                if base == consensus_base[p]:
                    score += w[p]
            codon_scores.append(score / total_w)
        identities.append(sum(codon_scores)/n)
    return identities


def compute_aa_identity(records, length):
    n = len(records)
    identities = []
    for start in range(0, length, 3):
        codons = [str(r.seq[start:start+3]) for r in records]
        # translate, handling gaps
        aas = []
        for c in codons:
            if '-' in c:
                aas.append('X')
            else:
                try:
                    aas.append(str(Seq(c).translate()))
                except:
                    aas.append('X')
        # consensus aa (ignore X)
        aa_counts = Counter([a for a in aas if a != 'X'])
        consensus_aa = aa_counts.most_common(1)[0][0] if aa_counts else 'X'
        # compute % consensus aa
        match_count = sum(1 for a in aas if a == consensus_aa and a != 'X')
        valid_count = sum(1 for a in aas if a != 'X')
        identity = match_count / valid_count if valid_count > 0 else 0
        identities.append(identity)
    return identities

def compute_blosum62_identity(records, length):
    n = len(records)
    identities = []
    for start in range(0, length, 3):
        codons = [str(r.seq[start:start+3]) for r in records]
        # translate, handling gaps
        aas = []
        for c in codons:
            if '-' in c:
                aas.append('X')
            else:
                try:
                    aas.append(str(Seq(c).translate()))
                except:
                    aas.append('X')
        # consensus aa (ignore X)
        aa_counts = Counter([a for a in aas if a != 'X'])
        consensus_aa = aa_counts.most_common(1)[0][0] if aa_counts else 'X'
        # compute normalized substitution scores
        scores = []
        for a in aas:
            if a == 'X' or consensus_aa == 'X':
                raw = MIN_BLOSUM
            else:
                raw = blosum62.get((a, consensus_aa), blosum62.get((consensus_aa, a), MIN_BLOSUM))
            norm = (raw - MIN_BLOSUM) / (MAX_BLOSUM - MIN_BLOSUM)
            scores.append(norm)
        identities.append(sum(scores)/n)
    return identities

def find_poor_codons(scores, threshold):
    return [i for i, s in enumerate(scores) if s < threshold]

def find_poor_regions(poor_codons):
    regions = []
    if not poor_codons:
        return regions
    start = prev = poor_codons[0]
    for idx in poor_codons[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            regions.append((start, prev))
            start = prev = idx
    regions.append((start, prev))
    return regions


def write_gff(regions, outfile):
    with open(outfile, 'w') as gf:
        gf.write("##gff-version 3\n")
        for idx, (s, e) in enumerate(regions, 1):
            nt_start = s*3 + 1
            nt_end   = (e+1)*3
            gf.write(f"chr1\tmark_poor\tlow_quality_region\t{nt_start}\t{nt_end}\t.\t+\t.\tID=poor{idx}\n")


def mask_alignment(records, regions):
    masked = []
    bad = set()
    for s, e in regions:
        for cod in range(s, e+1):
            for p in range(cod*3, cod*3+3):
                bad.add(p)
    for r in records:
        seq = list(str(r.seq))
        for p in bad:
            seq[p] = 'N'
        r.seq = Seq(''.join(seq))
        masked.append(r)
    return masked


def main():
    args = parse_args()
    recs, length = load_alignment(args.input)
    poor_codon_groups = []
    regions = []

    if args.aa_threshold is not None:
        aa_scores = compute_aa_identity(recs, length)
        poor_codon_groups.append(set(find_poor_codons(aa_scores, args.aa_threshold)))
    if args.blosum62_threshold is not None:
        blosum_scores = compute_blosum62_identity(recs, length)
        poor_codon_groups.append(set(find_poor_codons(blosum_scores, args.blosum62_threshold)))
    if args.codons_threshold is not None:
        if args.use_weighted_nuc:
            w = [float(x) for x in args.weights.split(',')]
            scores = compute_weighted_nuc_identity(recs, length, w)
        else:
            scores = compute_codon_identity(recs, length)
        poor_codon_groups.append(set(find_poor_codons(scores, args.codons_threshold)))

    ## Compute intersection (AND logic) for all groups
    while len(poor_codon_groups) > 1:
        poor_codon_groups[0] &= poor_codon_groups.pop(1)
    poor_codons = poor_codon_groups[0] if poor_codon_groups else set()

    regions = find_poor_regions(sorted(poor_codons))

    if not regions:
        print(f"No poorly aligned regions found")
    else:
        print(f"Found {len(regions)} poor region(s):")
        for s, e in regions:
            print(f"  Codons {s+1}-{e+1} (nt {s*3+1}-{(e+1)*3})")

    if args.gff_out:
        write_gff(regions, args.gff_out)
        print(f"GFF3 annotations written to {args.gff_out}")
    if args.o:
        masked = mask_alignment(recs, regions)
        SeqIO.write(masked, args.o, 'fasta')
        print(f"Masked FASTA written to {args.o}")

if __name__ == '__main__':
    main()
