#!/usr/bin/env python3
import argparse
import re
from Bio import SeqIO

def parse_ltr_passlist(input_file):
    entries = []
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11 or parts[9].strip() == "unknown":
                continue
            chrom_loc = parts[0]
            strand = parts[8].replace("+", "+").replace("-", "-")
            match = re.match(r"(.*):(\d+)\.\.(\d+)", chrom_loc)
            if match:
                chrom, start, end = match.groups()
                start, end = int(start), int(end)
                # No swapping, coordinates as-is
                entries.append((chrom, start, end, strand))
    return entries

def write_bed(entries, bed_file):
    with open(bed_file, 'w') as out:
        for idx, (chrom, start, end, strand) in enumerate(entries, 1):
            # Write start and end as-is, no -1 adjustment
            out.write(f"{chrom}\t{start}\t{end}\tLTR_{idx}\t0\t{strand}\n")

def write_fasta(entries, fasta_file, output_file):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    with open(output_file, "w") as out_fasta:
        seen = set()
        for idx, (chrom, start, end, strand) in enumerate(entries, 1):
            if chrom not in fasta_dict:
                continue
            # Extract sequence using coordinates as-is, no start-1 here
            seq = fasta_dict[chrom].seq[start:end]
            header = f">LTR_{idx}_{chrom}:{start}-{end}"
            if header in seen:
                continue
            seen.add(header)
            out_fasta.write(f"{header}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Parse LTR .mod.pass.list to BED and FASTA")
    parser.add_argument("--input", required=True, help="Input .mod.pass.list file")
    parser.add_argument("--fasta", required=True, help="Reference FASTA genome")
    parser.add_argument("--bed", required=True, help="Output BED file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file")

    args = parser.parse_args()

    entries = parse_ltr_passlist(args.input)
    write_bed(entries, args.bed)
    write_fasta(entries, args.fasta, args.output_fasta)

if __name__ == "__main__":
    main()
