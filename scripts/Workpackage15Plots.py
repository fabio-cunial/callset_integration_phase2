# Warning: this has been vibe-coded very quickly.
#
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Usage: Workpackage15Plots.py <matrix.csv> <length_quantum> <title>
#
# Input is a precomputed heatmap-count matrix (as produced by Workpackage15Plots.java):
# rows are FREQUENCY (row i = frequency i) and columns are LENGTH bins (column j =
# lengths in [j*length_quantum, (j+1)*length_quantum)). Each cell holds the number of
# records with that (FREQUENCY, LENGTH-bin) combination.
#
# Shows a heatmap whose horizontal axis is LENGTH, vertical axis is FREQUENCY, and color
# is the log10 of the count in each cell.

MATRIX_CSV = sys.argv[1]
LENGTH_QUANTUM = float(sys.argv[2])
TITLE = sys.argv[3]

if LENGTH_QUANTUM <= 0:
    sys.exit("ERROR: length quantum must be positive.")

# The matrix is written with a trailing comma per row, which yields a spurious all-NaN
# last column; drop any fully-empty columns.
matrix = pd.read_csv(MATRIX_CSV, header=None).dropna(axis=1, how="all").to_numpy(dtype=float)

# The matrix has up to ~2M length-bin columns but is overwhelmingly empty. Crop to the
# bounding box of the non-zero cells so we only render the region that has data; this is
# what makes plotting fast (and keeps memory down).
nz_rows = np.flatnonzero(matrix.any(axis=1))
nz_cols = np.flatnonzero(matrix.any(axis=0))
if nz_rows.size == 0:
    sys.exit("Matrix is empty: nothing to plot.")
r0, r1 = nz_rows[0], nz_rows[-1] + 1
c0, c1 = nz_cols[0], nz_cols[-1] + 1
sub = matrix[r0:r1, c0:c1]

# log10 of the count; empty cells (count==0) are masked so they stay blank.
log_counts = np.log10(sub, out=np.full_like(sub, np.nan), where=sub > 0)
log_counts = np.ma.masked_invalid(log_counts)

fig, ax = plt.subplots(figsize=(10, 8))
# Cell edges in real units. Length bins span [j*quantum, (j+1)*quantum); each frequency
# row r holds "r samples", so it is centered on r with integer-wide [r-0.5, r+0.5) edges.
length_edges = (np.arange(c0, c1 + 1) * LENGTH_QUANTUM).astype(float)
freq_edges = np.arange(r0, r1 + 1) - 0.5
# Both axes are log10, which cannot include 0; clamp the bottom edges to small positives.
freq_edges[0] = max(freq_edges[0], 0.5)
length_edges[0] = max(length_edges[0], LENGTH_QUANTUM / 2)

mesh = ax.pcolormesh(length_edges, freq_edges, log_counts, cmap="viridis")
fig.colorbar(mesh, ax=ax, label="Number of records (log10)")
ax.set_xlabel("SVLEN")
ax.set_ylabel("Number of samples")
ax.set_title(TITLE)

# Both axes on a log10 scale (X = SVLEN, Y = number of samples).
ax.set_xscale("log")
ax.set_yscale("log")

# Clamp the view to the non-empty bounding box (no surrounding margins).
ax.set_xlim(length_edges[0], length_edges[-1])
ax.set_ylim(freq_edges[0], freq_edges[-1])

# Trim the surrounding figure margins.
fig.tight_layout()

plt.show()
