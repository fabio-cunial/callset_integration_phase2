# I have six histogram files with names histogram_1_0.csv, histogram_1_1.csv, histogram_1_2.csv, 
# histogram_2_0.csv, histogram_2_1.csv, histogram_2_2.csv, each of them containing the same number of rows and one integer per row. 
# Write a Python script that loads all such files and plots them 
# in a single stacked barchart, where each X coordinate in the barchart corresponds to row X in every file.

import matplotlib.pyplot as plt
import numpy as np
import os

input_dir = "/Users/fcunial/Downloads/TestKanpigLength/histograms"

x_labels = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]

series_colors = {
    'HET->HET': '#2e7d32',  # dark green
    'HOM->HOM': '#81c784',  # light green
    'HET->REF': '#c62828',  # distinct red
    'HET->HOM': '#b71c1c',  # dark red
    'HOM->REF': '#e53935',  # medium red
    'HOM->HET': '#ef9a9a',  # light red
}

plot_order = [
    (1, 1), (2, 2),          # green series first and consecutive
    (1, 0), (1, 2), (2, 0), (2, 1),  # red series next and consecutive
]

subplots = [
    ('all', 'All records'),
    ('ins', 'INS'),
    ('del', 'DEL'),
]

bar_width = 1.0
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

for ax, (prefix, subplot_title) in zip(axes, subplots):
    # Load data for this prefix
    data = {}
    for i in range(1, 3):
        for j in range(0, 3):
            filename = f'{input_dir}/histogram_{prefix}_{i}_{j}.csv'
            with open(filename, 'r') as f:
                data[(i, j)] = [float(line.strip()) for line in f if line.strip()]

    x = np.arange(len(data[(1, 0)]))
    if len(x) != len(x_labels):
        raise ValueError(f"Expected {len(x_labels)} bins, found {len(x)} bins in {prefix} histograms")

    # One cumulative baseline for all six files creates a true stacked bar chart.
    bottom = np.zeros(len(x))
    for i, j in plot_order:
        label = {
            (1, 1): 'HET->HET',
            (2, 2): 'HOM->HOM',
            (1, 0): 'HET->REF',
            (1, 2): 'HET->HOM',
            (2, 0): 'HOM->REF',
            (2, 1): 'HOM->HET',
        }[(i, j)]
        values = np.array(data[(i, j)], dtype=float)
        ax.bar(x, values, bar_width, label=label, bottom=bottom, color=series_colors.get(label))
        bottom += values

    ax.set_xlabel('SVLEN bin')
    ax.set_title(subplot_title)
    ax.margins(x=0)
    ax.set_xticks(x)
    ax.set_xticklabels([str(v) for v in x_labels], rotation=45, ha='right')
    ax.legend()

axes[0].set_ylabel('Fraction of records')
fig.suptitle('Kanpig on dipcall records. Union of all calls from 10 HPRC Y2 samples.')
plt.tight_layout()
plt.show()
