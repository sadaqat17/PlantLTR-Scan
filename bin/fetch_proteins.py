#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import re

# -----------------------
# Argument parser
# -----------------------
parser = argparse.ArgumentParser(description="Fetch protein sequences based on input IDs and annotation mapping")
parser.add_argument("-i", "--input", required=True, help="Input file with gene/CDS/mRNA IDs (tab-delimited, first column)")
parser.add_argument("-p", "--proteome", required=True, help="Proteome FASTA file")
parser.add_argument("-a", "--annotation", required=True, help="Annotation file (GFF3 or GTF)")
parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
args = parser.parse_args()

# -----------------------
# Step 1: Read annotation mapping
# -----------------------
mapping = {} 

with open(args.annotation) as ann_file:
    for line in ann_file:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        info = fields[8]
        # parse key=value pairs
        info_dict = {}
        for kv in re.split(";|,", info):
            if "=" in kv:
                k, v = kv.split("=", 1)
                info_dict[k.strip()] = v.strip()
            elif " " in kv: 
                parts = kv.strip().split()
                if len(parts) == 2:
                    info_dict[parts[0]] = parts[1].strip('"')

        # Try to link input IDs to protein_id
        if "protein_id" in info_dict:
            prot_id = info_dict["protein_id"]
            # ID or Parent can be input
            if "ID" in info_dict:
                mapping[info_dict["ID"]] = prot_id
            if "Parent" in info_dict:
                mapping[info_dict["Parent"]] = prot_id
        elif "ID" in info_dict:
            # sometimes gene->mRNA->CDS chain only
            mapping[info_dict["ID"]] = info_dict.get("ID") 

# -----------------------
# Step 2: Read input IDs
# -----------------------
input_ids = set()
with open(args.input) as f:
    next(f)  
    for line in f:
        if line.strip():
            gene_id = line.split("\t")[0].strip()
            for acc in gene_id.split("|"):
                input_ids.add(acc.strip())

# Map input IDs to protein IDs
protein_ids = set()
for inp in input_ids:
    if inp in mapping:
        protein_ids.add(mapping[inp])

# -----------------------
# Step 3: Fetch sequences
# -----------------------
found_records = []
seen_ids = set()

for record in SeqIO.parse(args.proteome, "fasta"):
    for pid in protein_ids:
        if pid in record.id and pid not in seen_ids:
            found_records.append(record)
            seen_ids.add(pid)
            break

# -----------------------
# Step 4: Save output
# -----------------------
SeqIO.write(found_records, args.output, "fasta")
print(f"{len(found_records)} protein sequences saved to {args.output}")
