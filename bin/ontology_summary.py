#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

# -----------------------
# Argument parser
# -----------------------
parser = argparse.ArgumentParser(
    description="Classify input file by ontology (BP, CC, MF), count 'desc', and save sorted TSV files"
)
parser.add_argument(
    "-i", "--input", required=True, help="Input TSV file containing 'ontology' and 'desc' columns"
)
parser.add_argument(
    "-o", "--output_prefix", required=True, help="Output prefix for generated TSV files"
)
args = parser.parse_args()

# -----------------------
# Load input file
# -----------------------
def detect_delimiter(file_path, n_lines=5):
    with open(file_path) as f:
        lines = [next(f) for _ in range(n_lines)]
    delimiters = ['\t', ',', ';', ' ']
    scores = {d: sum(line.count(d) for line in lines) for d in delimiters}
    return max(scores, key=scores.get)

delimiter = detect_delimiter(args.input)
print(f"Detected delimiter: '{delimiter}'")

try:
    df = pd.read_csv(args.input, delimiter=delimiter)
except Exception as e:
    print(f"Failed to read {args.input}: {e}")
    sys.exit(1)

if 'ontology' not in df.columns or 'desc' not in df.columns:
    print("Input file must contain 'ontology' and 'desc' columns")
    sys.exit(1)

# -----------------------
# Function to filter, count, and save
# -----------------------
def process_ontology(df, ontology, output_file):
    filtered_df = df[df['ontology'] == ontology]
    if filtered_df.empty:
        print(f"No rows found for ontology '{ontology}'")
        return
    counts = filtered_df['desc'].value_counts().reset_index()
    counts.columns = ['desc', 'count']
    counts = counts.sort_values(by='count', ascending=False)
    counts.to_csv(output_file, sep='\t', index=False)
    print(f"{output_file} written successfully")

# -----------------------
# Process each ontology
# -----------------------
ontologies = ['BP', 'CC', 'MF']
for ont in ontologies:
    out_file = f"{args.output_prefix}_{ont}.tsv" 
    process_ontology(df, ont, out_file)

print("All done!")
