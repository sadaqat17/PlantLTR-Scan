import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Filter gffcompare tracking file and output unique Gene/mRNA accession + LTR coordinates + class code")
parser.add_argument("-i", "--input", required=True, help="Input .tracking file")
parser.add_argument("-o", "--output", required=True, help="Output filtered file")
args = parser.parse_args()

# Allowed class codes
allowed_codes = {"=", "c", "m", "n", "j", "o", "e", "i"}

seen = set() 

with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    # Write header
    outfile.write("Gene_ID\tLTR_coordinates\tClass_code\n")
    
    for line in infile:
        if line.strip():
            parts = line.split("\t")
            if len(parts) > 4:
                class_code = parts[3].strip()
                if class_code in allowed_codes:
                    gene_acc = parts[2].strip() 
                    ltr_info = parts[4].split("|")[1] if "|" in parts[4] else parts[4]
                    triple = (gene_acc, ltr_info, class_code)
                    if triple not in seen:
                        seen.add(triple)
                        outfile.write(f"{gene_acc}\t{ltr_info}\t{class_code}\n")

print(f"Filtered file saved to: {args.output}")
