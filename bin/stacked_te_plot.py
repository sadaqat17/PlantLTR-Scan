import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

def generate_multipanel_stacked_bar_plot(input_file, output_file, panel_size=10):
    # Read input, skip commented header lines starting with '#'
    df = pd.read_csv(input_file, sep="\t", comment='#', header=None,
                     names=["LTR_loc", "Category", "Motif", "TSD", "5_TSD", "3_TSD",
                            "Internal", "Identity", "Strand", "SuperFamily", "TE_type", "Insertion_Time"])

    # Extract chromosome from LTR_loc
    df['Chromosome'] = df['LTR_loc'].str.extract(r'^(.*?):')

    # Group by Chromosome and SuperFamily
    grouped = df.groupby(['Chromosome', 'SuperFamily']).size().unstack(fill_value=0)

    if grouped.empty:
        print("No data to plot. Check input file contents.")
        return

    chromosomes = grouped.index.tolist()
    total_chromosomes = len(chromosomes)
    num_panels = math.ceil(total_chromosomes / panel_size)

    # Setup figure size and subplots grid
    ncols = 1
    nrows = num_panels
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 5 * nrows), squeeze=False)

    cmap = cm.get_cmap('Pastel1', len(grouped.columns))
    colors = [cmap(i) for i in range(len(grouped.columns))]

    for panel_idx in range(num_panels):
        ax = axs[panel_idx, 0]

        start_idx = panel_idx * panel_size
        end_idx = min(start_idx + panel_size, total_chromosomes)
        panel_chromosomes = chromosomes[start_idx:end_idx]
        panel_data = grouped.loc[panel_chromosomes]

        bottom = np.zeros(len(panel_data))
        x_labels = panel_data.index.tolist()
        x = np.arange(len(x_labels))

        for idx, column in enumerate(grouped.columns):
            heights = panel_data[column].values
            bars = ax.bar(x, heights, bottom=bottom, label=column if panel_idx == 0 else None, color=colors[idx])

            # Add value labels in center of bars
            for i, bar in enumerate(bars):
                if heights[i] > 0:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bottom[i] + heights[i] / 2,
                        str(int(heights[i])),
                        ha='center', va='center', fontsize=8, color='black'
                    )

            bottom += heights

        ax.set_xlabel("Chromosome", fontsize=12)
        ax.set_ylabel("Count of TEs", fontsize=12)
        ax.set_title(f"Chromosomes {start_idx+1} to {end_idx}", fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=45, ha='right')

    # Only add one legend for all panels
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, title="SuperFamily", fontsize=10, title_fontsize=11,
               loc='upper right')

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave space for legend on right
    plt.savefig(output_file, format='svg', bbox_inches='tight', transparent=True)
    print(f"Multi-panel stacked bar plot saved as {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate multi-panel stacked bar plot of TE types per chromosome")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input LTR summary file (TSV)")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output SVG file for the bar plot")
    parser.add_argument('-p', '--panel_size', type=int, default=10, help="Number of chromosomes per panel (default 10)")
    args = parser.parse_args()

    generate_multipanel_stacked_bar_plot(args.input, args.output, args.panel_size)

if __name__ == "__main__":
    main()
