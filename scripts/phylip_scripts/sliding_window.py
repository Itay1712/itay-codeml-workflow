#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Split aligned FASTA into sliding windows.")
    parser.add_argument("output_pattern", help="Output pattern, e.g. windows/window_{start}-{end}.fasta")
    parser.add_argument("fasta", help="Aligned input FASTA file")
    parser.add_argument("-W", "--window", type=int, default=700, help="Window size (default: 700)")
    parser.add_argument("-s", "--step", type=int, default=500, help="Step size (default: 500)")
    return parser.parse_args()

def main():
    args = parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        raise ValueError("No sequences found.")

    aln_len = len(records[0].seq)
    for r in records:
        if len(r.seq) != aln_len:
            raise ValueError(f"All sequences must be aligned to the same length (got {len(r.seq)} vs {aln_len})")

    for start in range(0, aln_len, args.step):
        end = start + args.window
        if end > aln_len:
            break  # Skip incomplete window

        sliced_records = []
        for rec in records:
            new_rec = rec[start:end]
            new_rec.id = rec.id
            new_rec.description = f"{rec.description} [{start+1}-{end}]"
            sliced_records.append(new_rec)

        out_path = args.output_pattern.replace("{start}", str(start + 1)).replace("{end}", str(end))
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        SeqIO.write(sliced_records, out_path, "fasta")
        print(f"Wrote: {out_path}")

if __name__ == "__main__":
    main()