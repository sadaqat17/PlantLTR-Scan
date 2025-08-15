import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import math

def generate_multi_panel_bar_plot_single_svg(input_file, output_file):
    # Read the input file, skipping header
    df = pd.read_csv(input_file, sep="\t", skiprows=1, names=["Seq_ID", "Start-End", "Length", "Strand", "Score", "Similarity"])

    # Count occurrences of each chromosome
    chromosome_counts = Counter(df['Seq_ID'])
    chromosomes = list(chromosome_counts.keys())
    counts = list(chromosome_counts.values())

    panel_size = 10
    num_panels = math.ceil(len(chromosomes) / panel_size)

    # Create figure with subplots: vertical stack of num_panels rows, 1 col
    fig, axes = plt.subplots(num_panels, 1, figsize=(12, 4 * num_panels), constrained_layout=True)

    # If only one panel, axes is not array, make it iterable
    if num_panels == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        start_idx = i * panel_size
        end_idx = start_idx + panel_size
        panel_chromosomes = chromosomes[start_idx:end_idx]
        panel_counts = counts[start_idx:end_idx]

        bars = ax.bar(panel_chromosomes, panel_counts, color='skyblue')

        # Add value labels on bars
        for bar, count in zip(bars, panel_counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
                    str(count), ha='center', va='center', fontsize=11, color='black', weight='bold')

        ax.set_xlabel('Chromosome', fontsize=12)
        ax.set_ylabel('Count of TEs', fontsize=12)
        ax.set_title(f'Count of TEs per Chromosome (Panel {i+1} of {num_panels})', fontsize=14)

        # Fix for the warning: set ticks before setting tick labels
        ax.set_xticks(range(len(panel_chromosomes)))
        ax.set_xticklabels(panel_chromosomes, rotation=45, ha='right')

    plt.savefig(output_file, format='svg')
    print(f"Multi-panel bar plot saved as {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate multi-panel bar plots of TEs per chromosome in one SVG")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input file with TE data")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output SVG file for the multi-panel plot")
    args = parser.parse_args()

    generate_multi_panel_bar_plot_single_svg(args.input, args.output)

if __name__ == '__main__':
    main()
