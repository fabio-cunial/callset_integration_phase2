#!/usr/bin/env python3
"""
Usage: plot_af_histograms.py <file1> <file2> [output.png]

Loads two text files, each containing one allele frequency (0-1) per line,
and plots both count histograms on the same axes with 100 fixed bins.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


def load_af(path):
    values = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                values.append(float(line))
    return np.array(values)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    outfile = sys.argv[3] if len(sys.argv) > 3 else None

    af1 = load_af(file1)
    af2 = load_af(file2)

    bins = np.linspace(0.0, 1.0, 101)  # 100 equal-width bins in [0, 1]

    def reverse_cdf(values, bins):
        counts, edges = np.histogram(values, bins=bins)
        # cumulative sum from the highest-frequency bin downward, normalized
        cumulative = np.cumsum(counts[::-1])[::-1] / len(values)
        return edges[:-1], cumulative  # left edge of each bin

    fig, ax = plt.subplots(figsize=(10, 5))
    for af, label in [(af1, file1), (af2, file2)]:
        x, y = reverse_cdf(af, bins)
        ax.step(x, y, where="post", label=label)
    ax.set_xlabel("Allele frequency")
    ax.set_ylabel("Fraction (AF >= bin)")
    ax.set_title("Cumulative allele frequency distribution (high to low)")
    ax.legend()
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150)
        print(f"Saved to {outfile}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
