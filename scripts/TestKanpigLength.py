# I have six histogram files with names histogram_1_0.csv, histogram_1_1.csv, histogram_1_2.csv, 
# histogram_2_0.csv, histogram_2_1.csv, histogram_2_2.csv, each of them containing the same number of rows and one integer per row. 
# Write a Python script that loads all such files and plots them 
# in a single stacked barchart, where each X coordinate in the barchart corresponds to row X in every file.

import matplotlib.pyplot as plt
import numpy as np
import os

input_dir = "/Users/fcunial/Downloads/TestKanpigLength/histograms"

# Load data from files
data = {}
for i in range(1, 3):   # Loop through 1 and 2
    for j in range(0, 3):   # Loop through 0, 1, and 2
        filename = f'{input_dir}/histogram_{i}_{j}.csv'
        with open(filename, 'r') as f:
            data[(i, j)] = [float(line.strip()) for line in f if line.strip()]

# Plot bars
x = np.arange(len(data[(1, 0)]))
x_labels = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
if len(x) != len(x_labels):
    raise ValueError(f"Expected {len(x_labels)} bins, found {len(x)} bins in input histograms")

bar_width = 1.0
fig, ax = plt.subplots()

series_colors = {
    '1->1': '#2e7d32',  # dark green
    '2->2': '#81c784',  # light green
    '1->0': '#c62828',  # distinct red
    '1->2': '#b71c1c',  # dark red
    '2->0': '#e53935',  # medium red
    '2->1': '#ef9a9a',  # light red
}

plot_order = [
    (1, 1), (2, 2),          # green series first and consecutive
    (1, 0), (1, 2), (2, 0), (2, 1),  # red series next and consecutive
]

# One cumulative baseline for all six files creates a true stacked bar chart.
bottom = np.zeros(len(x))
for i, j in plot_order:
    label = f'{i}->{j}'
    values = np.array(data[(i, j)], dtype=float)
    ax.bar(x, values, bar_width, label=label, bottom=bottom, color=series_colors.get(label))
    bottom += values
ax.set_xlabel('SVLEN bin')
ax.set_ylabel('Fraction of records')
ax.set_title('Kanpig on dipcall SVs. INS and DEL from 10 HPRC Y2 samples')
ax.legend()
ax.margins(x=0)
ax.set_xticks(x)
ax.set_xticklabels([str(v) for v in x_labels], rotation=45, ha='right')
plt.tight_layout()
plt.show()
