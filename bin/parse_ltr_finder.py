import re
import argparse

def parse_ltr_finder(input_file, output_file):
    with open(input_file, 'r') as f:
        data = f.read()

    # Split file into sections starting with ">Sequence:"
    sections = re.split(r'(?=>Sequence:\s*[^\s]+)', data)

    total_ltrs = 0
    with open(output_file, 'w') as out:
        out.write("Seq_ID\tStart-End\tLength\tStrand\tScore\tSimilarity\n")

        for section in sections:
            if not section.strip():
                continue

            # Extract sequence ID from header
            seq_id_match = re.search(r'>Sequence:\s*([^\s]+)', section)
            seq_id = seq_id_match.group(1) if seq_id_match else "NA"

            # Find all LTR entries for this sequence
            entries = re.findall(r'\[\d+\](.*?)Details of PPT', section, re.DOTALL)
            total_ltrs += len(entries)

            for entry in entries:
                # Location and Strand
                loc_match = re.search(r'Location\s*:\s*(\d+)\s*-\s*(\d+).*Strand:([+-])', entry)
                if loc_match:
                    start, end, strand = loc_match.groups()
                    length = int(end) - int(start) + 1
                    location = f"{start}-{end}"
                else:
                    location = length = strand = "NA"

                # Score and Similarity
                score_match = re.search(r'Score\s*:\s*(\d+)\s+\[LTR region similarity:([\d\.]+)\]', entry)
                if score_match:
                    score, sim = score_match.groups()
                else:
                    score = sim = "NA"

                out.write(f"{seq_id}\t{location}\t{length}\t{strand}\t{score}\t{sim}\n")

    print(f"Total LTR found by LTR_Finder: {total_ltrs}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse LTR_FINDER output to tabular format.")
    parser.add_argument("-i", "--input", required=True, help="LTR_FINDER output file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    parse_ltr_finder(args.input, args.output)
