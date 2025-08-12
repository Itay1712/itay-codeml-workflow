#!/usr/bin/env python3
import sys

def parse_fasta(filename):
    """
    Parses a FASTA file and returns a dictionary with sequence identifiers as keys 
    and sequences as values.
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save the previous sequence if it exists
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                # Use the first token (up to the first whitespace) as the identifier
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        # Save the last sequence
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    return sequences

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 fasta_to_phylip.py inputfasta output.phy")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_phylip = sys.argv[2]
    
    sequences = parse_fasta(input_fasta)
    
    if not sequences:
        print("Error: No sequences found in the FASTA file.")
        sys.exit(1)
    
    # Check that all sequences have equal lengths (they must be aligned)
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        print("Error: Sequences have different lengths. They must be aligned.")
        sys.exit(1)
    
    seq_length = lengths.pop()
    num_sequences = len(sequences)
    
    with open(output_phylip, "w") as out:
        # Write the header: number of sequences and sequence length
        out.write(f"{num_sequences} {seq_length}\n")
        for header, seq in sequences.items():
            # Truncate or pad the identifier to 10 characters
            label = header[:10].ljust(10)
            out.write(f"{label}  {seq}\n")

if __name__ == "__main__":
    main()