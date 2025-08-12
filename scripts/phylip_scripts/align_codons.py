from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os
import sys


def translate_sequences(nuc_records):
    """
    Translate nucleotide sequences into amino acids.
    Assuming sequence is in-frame and no introns
    """
    prot_records = []
    for record in nuc_records:
        full_prot_len = len(record.seq) // 3
        prot_seq = record.seq.translate(to_stop=True)
        if len(prot_seq) == full_prot_len-1:
            print(f"[Warning] Sequence ended with a stop codon - translated {len(prot_seq)} aa from {full_prot_len} codons.")
        elif len(prot_seq) < full_prot_len-1:
            print(f"[Warning] Early stop codon in sequence '{record.id}' — translated only {len(prot_seq)} aa of expected {full_prot_len}.")

        prot_record = SeqRecord(prot_seq, id=record.id, description="Translated")
        prot_records.append(prot_record)
    return prot_records

def run_mafft(input_fasta, output_fasta):
    """
    Run MAFFT alignment on the input FASTA file and write the alignment to output_fasta.
    """
    cmd = ["mafft", "--auto", input_fasta]
    with open(output_fasta, "w") as outf:
        subprocess.run(cmd, stdout=outf, check=True)

def back_translate(aligned_prot_records, orig_nuc_dict):
    """
    Convert protein alignment back to codon alignment.
    For each gap in the protein alignment, insert three gaps (---) in the codon alignment.
    """
    aligned_nuc_records = []
    for prot_record in aligned_prot_records:
        nuc_seq = str(orig_nuc_dict[prot_record.id].seq)
        codon_index = 0  ## pointer for original codon
        aligned_nuc_seq = ""
        for aa in str(prot_record.seq):
            if aa == "-":
                aligned_nuc_seq += "---"  ## gap for a missing codon
            else:
                ## extract corresponding codon from original nucleotide sequence
                codon = nuc_seq[codon_index*3 : codon_index*3 + 3]
                aligned_nuc_seq += codon
                codon_index += 1
        aligned_nuc_records.append(SeqRecord(Seq.Seq(aligned_nuc_seq), 
                                               id=prot_record.id,
                                               description="Codon alignment"))
    return aligned_nuc_records

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 fasta_to_phylip.py inputfasta output.fasta")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_aligned_fasta = sys.argv[2]
    
    nuc_records = list(SeqIO.parse(input_fasta, "fasta"))
    orig_nuc_dict = {record.id: record for record in nuc_records}

    prot_records = translate_sequences(nuc_records)
    temp_prot_fasta = "temp_prot.fasta"
    SeqIO.write(prot_records, temp_prot_fasta, "fasta")

    aligned_prot_fasta = "aligned_prot.fasta"
    run_mafft(temp_prot_fasta, aligned_prot_fasta)
    aligned_prot_records = list(SeqIO.parse(aligned_prot_fasta, "fasta"))

    aligned_nuc_records = back_translate(aligned_prot_records, orig_nuc_dict)
    SeqIO.write(aligned_nuc_records, output_aligned_fasta, "fasta") ## final output

    os.remove(temp_prot_fasta)
    os.remove(aligned_prot_fasta)

if __name__ == "__main__":
    main()
