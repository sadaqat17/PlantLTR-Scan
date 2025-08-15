import argparse

def extract_chr_id(input_file, output_file):
    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if line.startswith("#") or not line.strip():
                f_out.write(line)
                continue
            cols = line.strip().split("\t")
            # Extract chromosome ID from column 1
            chr_id = "_".join(cols[0].split("_")[2:])
            cols[0] = chr_id
            f_out.write("\t".join(cols) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace column 1 with chromosome ID from LTR index")
    parser.add_argument("-i", "--input", required=True, help="Input GFF3-like file")
    parser.add_argument("-o", "--output", required=True, help="Output file with chr IDs in column 1")
    args = parser.parse_args()

    extract_chr_id(args.input, args.output)
