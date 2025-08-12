#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import sys

def parse_mask_file(mask_file):
    """Parse tab-delimited recombination regions: header \t start \t end."""
    mask_regions = {}
    with open(mask_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue
            seq_id, start, end = parts
            start, end = int(start), int(end)
            if seq_id not in mask_regions:
                mask_regions[seq_id] = []
            mask_regions[seq_id].append((start, end))
    return mask_regions

def apply_mask(seq, regions, mask_char='N'):
    """Mask 1-based inclusive regions in sequence string."""
    seq = list(seq)
    for start, end in regions:
        for i in range(start - 1, end):
            if 0 <= i < len(seq):
                seq[i] = mask_char
    return ''.join(seq)

def main():
    parser = argparse.ArgumentParser(description="Mask recombination regions in FASTA")
    parser.add_argument("fasta", help="Input FASTA alignment file")
    parser.add_argument("mask_file", help="TSV file with: header TAB start TAB end")
    parser.add_argument("output", help="Masked FASTA output path")
    parser.add_argument("--mask-char", default='N', help="Character to use for masking (default: N)")

    args = parser.parse_args()

    mask_regions = parse_mask_file(args.mask_file)
    masked_records = []

    for record in SeqIO.parse(args.fasta, "fasta"):
        regions = mask_regions.get(record.id, [])
        if regions:
            print(f"Masking {len(regions)} region(s) in {record.id}")
            record.seq = apply_mask(str(record.seq), regions, args.mask_char)
        masked_records.append(record)

    SeqIO.write(masked_records, args.output, "fasta")
    print(f"Masked alignment saved to: {args.output}")

if __name__ == "__main__":
    main()