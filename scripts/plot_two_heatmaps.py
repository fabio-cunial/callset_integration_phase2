#!/usr/bin/env python3
"""Load two TSV matrices, sort their rows by vector similarity, and show them
as two heatmap panels in a single figure."""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list


def load_matrix(path):
    """Load a TSV file of floats into a 2D numpy array."""
    return np.loadtxt(path, delimiter="\t", dtype=float)


def sort_rows_by_similarity(matrix):
    """Order rows by hierarchical clustering.

    Rows are clustered with Ward linkage on Euclidean distances, and the
    dendrogram leaf order is used to arrange them so that similar rows sit
    next to each other.
    """
    if matrix.shape[0] < 2:
        return matrix
    order = leaves_list(linkage(matrix, method="ward"))
    return matrix[order]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv1", help="first TSV matrix file")
    parser.add_argument("tsv2", help="second TSV matrix file")
    args = parser.parse_args()

    paths = [args.tsv1, args.tsv2]
    row_limits = [7000, 1600] #[10000, 10000] #[7000, 1600]
    matrices = [
        sort_rows_by_similarity(load_matrix(p)[:, 18:38])[:n]
        for p, n in zip(paths, row_limits)
    ]

    vmin = min(matrix.min() for matrix in matrices)
    vmax = 100

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, matrix, path in zip(axes, matrices, paths):
        im = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
        ax.set_title(path)
        ax.set_xlabel("column")
        ax.set_ylabel("row (sorted)")
        fig.colorbar(im, ax=ax)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
