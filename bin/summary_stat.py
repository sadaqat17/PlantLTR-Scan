import argparse
from collections import defaultdict


def parse_fasta_headers(file_path):
    """Parse fasta headers and return clade info with sequence lengths."""
    data = {"Copia": defaultdict(list), "Gypsy": defaultdict(list)}

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                try:
                    coords = header.split(":")[1].split("#")[0]
                    start, end = map(int, coords.split("-")[-2:])
                    length = abs(end - start) + 1

                    family_info = header.split("#")[1].split("/")
                    family = family_info[1]  # Copia or Gypsy
                    clade = family_info[2]   # e.g., Ivana, Reina

                    if family in data:
                        data[family][clade].append(length)
                except Exception as e:
                    print(f"Skipping malformed header: {header} ({e})")
    return data


def compute_stats(lengths):
    if not lengths:
        return (0, 0, 0, 0)
    count = len(lengths)
    min_len = min(lengths)
    max_len = max(lengths)
    avg_len = sum(lengths) / count
    return (count, min_len, max_len, avg_len)


def write_stats(data, out_prefix):
    out_file = f"{out_prefix}_stat.txt"
    with open(out_file, "w") as out:
        out.write("Family\tClade\tCount\tMin_seq_length\tMax_seq_length\tAverage_seq_length\n")
        for family in ["Copia", "Gypsy"]:
            for clade, lengths in data[family].items():
                count, min_len, max_len, avg_len = compute_stats(lengths)
                out.write(f"{family}\t{clade}\t{count}\t{min_len}\t{max_len}\t{avg_len:.2f}\n")

        # Family-level totals
        out.write("\nSummary:\n")
        family_lengths = {}
        for family in ["Copia", "Gypsy"]:
            all_lengths = []
            for lengths in data[family].values():
                all_lengths.extend(lengths)
            total_count, min_len, max_len, avg_len = compute_stats(all_lengths)
            family_lengths[family] = all_lengths
            out.write(f"{family}\tTotal Sequences: {total_count}\tMin_seq_length: {min_len}\tMax_seq_length: {max_len}\tAverage_seq_length: {avg_len:.2f}\n")

        # Overall summary
        combined_lengths = family_lengths["Copia"] + family_lengths["Gypsy"]
        total_count, min_len, max_len, avg_len = compute_stats(combined_lengths)
        out.write(f"Overall\tTotal Sequences: {total_count}\tMin_seq_length: {min_len}\tMax_seq_length: {max_len}\tAverage_seq_length: {avg_len:.2f}\n")

    print(f"Saved combined stats to {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Compute Copia and Gypsy clade statistics from fasta headers.")
    parser.add_argument("-i", "--input", required=True, help="Input fasta file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output file prefix")
    args = parser.parse_args()

    data = parse_fasta_headers(args.input)
    write_stats(data, args.out_prefix)


if __name__ == "__main__":
    main()
