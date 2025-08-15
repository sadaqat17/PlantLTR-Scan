#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import subprocess
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from pathlib import Path

# -----------------------------
# Argument Parser
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Parallel LTR_FINDER for multi-FASTA genomes.")
    parser.add_argument("-i", "--input", required=True, help="Input multi-FASTA genome file")
    parser.add_argument("-o", "--output", required=True, help="Output combined LTR_FINDER result file (*.finder.scn)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Total threads to use (default: 4)")
    return parser.parse_args()

# -----------------------------
# Sequence Splitter
# -----------------------------
def split_fasta(input_fasta, output_dir):
    input_path = Path(input_fasta)
    os.makedirs(output_dir, exist_ok=True)

    split_files = []

    with open(input_path, "r") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            seq_id = record.id.replace("/", "_")
            out_file = Path(output_dir) / f"{seq_id}.fa"
            SeqIO.write(record, out_file, "fasta")
            split_files.append(str(out_file))
    return split_files

# -----------------------------
# LTR_FINDER Wrapper
# -----------------------------
def run_ltr_finder(args):
    fasta_file, output_file = args
    cmd = [
        "ltr_finder",
        "-D", "15000", "-d", "1000", "-L", "7000", "-l", "100", "-p", "20",
        "-C", "-M", "0.9", fasta_file
    ]
    try:
        with open(output_file, "w") as out_f:
            subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        print(f"Error running LTR_FINDER on {fasta_file}", file=sys.stderr)

# -----------------------------
# Main
# -----------------------------
def main():
    args = parse_args()

    input_fasta = args.input
    output_file = args.output
    total_threads = args.threads

    splitter_dir = "tmp/ltr_finder/split"
    results_dir = "tmp/ltr_finder/results"
    os.makedirs(splitter_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    print(f"Splitting input FASTA: {input_fasta}")
    fasta_chunks = split_fasta(input_fasta, splitter_dir)

    print(f"Running LTR_FINDER using {total_threads} threads "
          f"with 2 threads per job (running {total_threads//2} jobs in parallel)\n")
    
    # Prepare arguments for multiprocessing
    finder_tasks = []
    for fasta_path in fasta_chunks:
        out_filename = Path(results_dir) / (Path(fasta_path).stem + ".finder.scn")
        finder_tasks.append((fasta_path, str(out_filename)))

    # Run ltr_finder in parallel (2 threads per job)
    max_jobs = total_threads // 2 if total_threads >= 2 else 1
    with Pool(processes=max_jobs) as pool:
        pool.map(run_ltr_finder, finder_tasks)

    # Combine results
    print(f"\nCombining output files into {output_file}")
    with open(output_file, "w") as wout:
        for _, result_path in finder_tasks:
            if os.path.exists(result_path):
                with open(result_path, "r") as infile:
                    shutil.copyfileobj(infile, wout)
            else:
                print(f"Missing result for: {result_path}", file=sys.stderr)

    print("Done!")

if __name__ == "__main__":
    main()

