import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

# Configuration
COHORT_NAME = 'BI'
CALLER_ID = 2  # 0=pav, 1=pbsv, 2=sniffles
SELECTED_SV_TYPE = 'DEL'  # DEL, INS, DUP, INV, BND, SUB, UNK

SV_TYPES_PRIMARY = ['DEL', 'INS', 'DUP', 'INV']
SV_TYPES_SECONDARY = ['BND', 'SUB', 'UNK']
MAX_COVERAGE = 50
MAX_NRECORDS = 40000
CALLER_NAMES = ['PAV', 'PBSV', 'Sniffles']
HEATMAP_BINS_X = 50
HEATMAP_BINS_Y = 50

if SELECTED_SV_TYPE not in SV_TYPES_PRIMARY + SV_TYPES_SECONDARY:
    raise ValueError(
        f'Unsupported SV type: {SELECTED_SV_TYPE}. '
        f'Valid values: {SV_TYPES_PRIMARY + SV_TYPES_SECONDARY}'
    )

# Load data
A = pd.read_csv(f'/Users/fcunial/Downloads/svqc/plot/{COHORT_NAME}/counts.csv', header=None)

# Meaning of the columns of `counts.csv` (numbers are offsets from the first 
# column, which contains the coverage):
#
# - >=20
#   - whole genome
#     - DEL
#       1 pav
#       2 pbsv
#       3 sniffles
#     - INS
#       4 pav
#       5 pbsv
#       6 sniffles
#     - DUP
#       7 pav
#       8 pbsv
#       9 sniffles
#     - INV
#       10 pav
#       11 pbsv
#       12 sniffles
#   - inside TR
#     - DEL
#       13 pav
#       14 pbsv
#       15 sniffles
#     - INS
#       16 pav
#       17 pbsv
#       18 sniffles
#     - DUP
#       19 pav
#       20 pbsv
#       21 sniffles
#     - INV
#       22 pav
#       23 pbsv
#       24 sniffles
#   - outside TR
#     - DEL
#       25 pav
#       26 pbsv
#       27 sniffles
#     - INS
#       28 pav
#       29 pbsv
#       30 sniffles
#     - DUP
#       31 pav
#       32 pbsv
#       33 sniffles
#     - INV
#       34 pav
#       35 pbsv
#       36 sniffles
#
# - >=50
#   - whole genome
#     - DEL
#       37 pav
#       38 pbsv
#       39 sniffles
#     - INS
#       40 pav
#       41 pbsv
#       42 sniffles
#     - DUP
#       43 pav
#       44 pbsv
#       45 sniffles
#     - INV
#       46 pav
#       47 pbsv
#       48 sniffles
#   - inside TR
#     - DEL
#       49 pav
#       50 pbsv
#       51 sniffles
#     - INS
#       52 pav
#       53 pbsv
#       54 sniffles
#     - DUP
#       55 pav
#       56 pbsv
#       57 sniffles
#     - INV
#       58 pav
#       59 pbsv
#       60 sniffles
#   - outside TR
#     - DEL
#       61 pav
#       62 pbsv
#       63 sniffles
#     - INS
#       64 pav
#       65 pbsv
#       66 sniffles
#     - DUP
#       67 pav
#       68 pbsv
#       69 sniffles
#     - INV
#       70 pav
#       71 pbsv
#       72 sniffles
#
# - Others
#   - whole genome
#    - BND
#       73 pav
#       74 pbsv
#       75 sniffles
#    - SUB
#       76 pav
#       77 pbsv
#       78 sniffles
#    - UNK
#       79 pav
#       80 pbsv
#       81 sniffles
#   - inside TR
#     - BND
#       82 pav
#       83 pbsv 
#       84 sniffles
#     - SUB
#       85 pav
#       86 pbsv
#       87 sniffles
#     - UNK
#       88 pav
#       89 pbsv
#       90 sniffles
#   - outside TR
#     - BND
#       91 pav
#       92 pbsv
#       93 sniffles
#     - SUB
#       94 pav  
#       95 pbsv
#       96 sniffles
#     - UNK
#       97 pav
#       98 pbsv
#       99 sniffles

fig = plt.figure(figsize=(15, 12))
x = A.iloc[:, 0].values  # First column is coverage


def _make_bin_edges(values, nbins):
    if values.size == 0:
        return np.linspace(0.0, 1.0, nbins + 1)
    vmin = np.min(values)
    vmax = np.max(values)
    if vmin == vmax:
        pad = max(abs(vmin) * 0.01, 0.5)
        vmin -= pad
        vmax += pad
    return np.linspace(vmin, vmax, nbins + 1)


def plot_count_heatmap(fig_obj, ax_obj, x_values, y_values):
    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)

    valid_mask = np.isfinite(x_values) & np.isfinite(y_values)
    x_valid = x_values[valid_mask]
    y_valid = y_values[valid_mask]

    x_edges = _make_bin_edges(x_valid, HEATMAP_BINS_X)
    y_edges = _make_bin_edges(y_valid, HEATMAP_BINS_Y)
    counts, _, _ = np.histogram2d(x_valid, y_valid, bins=[x_edges, y_edges])

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_under('white')
    vmax = max(1.0, float(np.max(counts)))
    norm = colors.Normalize(vmin=0.5, vmax=vmax)

    mesh = ax_obj.pcolormesh(
        x_edges,
        y_edges,
        counts.T,
        cmap=cmap,
        norm=norm,
        shading='auto',
    )
    fig_obj.colorbar(mesh, ax=ax_obj, label='Number of samples')







# >=20bp
# Subplot 1: whole genome
ax1 = plt.subplot(3, 3, 1)
y1 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'DEL':
    y1 = A.iloc[:, 1 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INS':
    y1 = A.iloc[:, 4 + CALLER_ID].values
if SELECTED_SV_TYPE == 'DUP':
    y1 = A.iloc[:, 7 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INV':
    y1 = A.iloc[:, 10 + CALLER_ID].values
plot_count_heatmap(fig, ax1, x, y1)
ax1.set_xlabel('Coverage')
ax1.set_ylabel('Number of calls')
ax1.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}>=20bp, whole genome')
ax1.set_box_aspect(1)
ax1.set_axisbelow(False)
ax1.grid(True, color='lightgray')
ax1.tick_params(labelsize=10)

# Subplot 4: inside TR
ax4 = plt.subplot(3, 3, 4)
y4 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'DEL':
    y4 = A.iloc[:, 13 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INS':
    y4 = A.iloc[:, 16 + CALLER_ID].values
if SELECTED_SV_TYPE == 'DUP':
    y4 = A.iloc[:, 19 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INV':
    y4 = A.iloc[:, 22 + CALLER_ID].values
plot_count_heatmap(fig, ax4, x, y4)
ax4.set_xlabel('Coverage')
ax4.set_ylabel('Number of calls')
ax4.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}>=20bp, inside TR')
ax4.set_box_aspect(1)
ax4.set_axisbelow(False)
ax4.grid(True, color='lightgray')
ax4.tick_params(labelsize=10)

# Subplot 7: outside TR
ax7 = plt.subplot(3, 3, 7)
y7 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'DEL':
    y7 = A.iloc[:, 25 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INS':
    y7 = A.iloc[:, 28 + CALLER_ID].values
if SELECTED_SV_TYPE == 'DUP':
    y7 = A.iloc[:, 31 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INV':
    y7 = A.iloc[:, 34 + CALLER_ID].values
plot_count_heatmap(fig, ax7, x, y7)
ax7.set_xlabel('Coverage')
ax7.set_ylabel('Number of calls')
ax7.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}>=20bp, outside TR')
ax7.set_box_aspect(1)
ax7.set_axisbelow(False)
ax7.grid(True, color='lightgray')
ax7.tick_params(labelsize=10)

# >=50bp
# Subplot 2: whole genome
ax2 = plt.subplot(3, 3, 2)
y2 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'DEL':
    y2 = A.iloc[:, 37 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INS':
    y2 = A.iloc[:, 40 + CALLER_ID].values
if SELECTED_SV_TYPE == 'DUP':
    y2 = A.iloc[:, 43 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INV':
    y2 = A.iloc[:, 46 + CALLER_ID].values
plot_count_heatmap(fig, ax2, x, y2)
ax2.set_xlabel('Coverage')
ax2.set_ylabel('Number of calls')
ax2.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}>=50bp, whole genome')
ax2.set_box_aspect(1)
ax2.set_axisbelow(False)
ax2.grid(True, color='lightgray')
ax2.tick_params(labelsize=10)

# Subplot 5: inside TR
ax5 = plt.subplot(3, 3, 5)
y5 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'DEL':
    y5 = A.iloc[:, 49 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INS':
    y5 = A.iloc[:, 52 + CALLER_ID].values
if SELECTED_SV_TYPE == 'DUP':
    y5 = A.iloc[:, 55 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INV':
    y5 = A.iloc[:, 58 + CALLER_ID].values
plot_count_heatmap(fig, ax5, x, y5)
ax5.set_xlabel('Coverage')
ax5.set_ylabel('Number of calls')
ax5.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}>=50bp, inside TR')
ax5.set_box_aspect(1)
ax5.set_axisbelow(False)
ax5.grid(True, color='lightgray')
ax5.tick_params(labelsize=10)

# Subplot 8: outside TR
ax8 = plt.subplot(3, 3, 8)
y8 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'DEL':
    y8 = A.iloc[:, 61 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INS':
    y8 = A.iloc[:, 64 + CALLER_ID].values
if SELECTED_SV_TYPE == 'DUP':
    y8 = A.iloc[:, 67 + CALLER_ID].values
if SELECTED_SV_TYPE == 'INV':
    y8 = A.iloc[:, 70 + CALLER_ID].values
plot_count_heatmap(fig, ax8, x, y8)
ax8.set_xlabel('Coverage')
ax8.set_ylabel('Number of calls')
ax8.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}>=50bp, outside TR')
ax8.set_box_aspect(1)
ax8.set_axisbelow(False)
ax8.grid(True, color='lightgray')
ax8.tick_params(labelsize=10)

# Other types
# Subplot 3: whole genome
ax3 = plt.subplot(3, 3, 3)
y3 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'BND':
    y3 = A.iloc[:, 73 + CALLER_ID].values
if SELECTED_SV_TYPE == 'SUB':
    y3 = A.iloc[:, 76 + CALLER_ID].values
if SELECTED_SV_TYPE == 'UNK':
    y3 = A.iloc[:, 79 + CALLER_ID].values
plot_count_heatmap(fig, ax3, x, y3)
ax3.set_xlabel('Coverage')
ax3.set_ylabel('Number of calls')
ax3.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}, whole genome')
ax3.set_box_aspect(1)
ax3.set_axisbelow(False)
ax3.grid(True, color='lightgray')
ax3.tick_params(labelsize=10)

# Subplot 6: inside TR
ax6 = plt.subplot(3, 3, 6)
y6 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'BND':
    y6 = A.iloc[:, 82 + CALLER_ID].values
if SELECTED_SV_TYPE == 'SUB':
    y6 = A.iloc[:, 85 + CALLER_ID].values
if SELECTED_SV_TYPE == 'UNK':
    y6 = A.iloc[:, 88 + CALLER_ID].values
plot_count_heatmap(fig, ax6, x, y6)
ax6.set_xlabel('Coverage')
ax6.set_ylabel('Number of calls')
ax6.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}, inside TR')
ax6.set_box_aspect(1)
ax6.set_axisbelow(False)
ax6.grid(True, color='lightgray')
ax6.tick_params(labelsize=10)

# Subplot 9: outside TR
ax9 = plt.subplot(3, 3, 9)
y9 = np.full(x.shape, np.nan)
if SELECTED_SV_TYPE == 'BND':
    y9 = A.iloc[:, 91 + CALLER_ID].values
if SELECTED_SV_TYPE == 'SUB':
    y9 = A.iloc[:, 94 + CALLER_ID].values
if SELECTED_SV_TYPE == 'UNK':
    y9 = A.iloc[:, 97 + CALLER_ID].values
plot_count_heatmap(fig, ax9, x, y9)
ax9.set_xlabel('Coverage')
ax9.set_ylabel('Number of calls')
ax9.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}, outside TR')
ax9.set_box_aspect(1)
ax9.set_axisbelow(False)
ax9.grid(True, color='lightgray')
ax9.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig(f'/Users/fcunial/Downloads/svqc/plot/{COHORT_NAME}/counts_{SELECTED_SV_TYPE}.png', dpi=300, bbox_inches='tight')
plt.show()









# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: heatmaps with X = number of DEL, Y = number of INS, cell = number
# of samples at those (DEL, INS) counts.
# Column 1: >=20bp panels.  Column 2: >=50bp panels.
# ─────────────────────────────────────────────────────────────────────────────

fig2 = plt.figure(figsize=(10, 12))

# >=20bp
# Subplot 1: whole genome
ax2_1 = fig2.add_subplot(3, 2, 1)
plot_count_heatmap(fig2, ax2_1,
                   A.iloc[:, 1 + CALLER_ID].values,
                   A.iloc[:, 4 + CALLER_ID].values)
ax2_1.set_xlabel('Number of DEL')
ax2_1.set_ylabel('Number of INS')
ax2_1.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=20bp, whole genome')
ax2_1.set_box_aspect(1)
ax2_1.set_axisbelow(False)
ax2_1.grid(True, color='lightgray')
ax2_1.tick_params(labelsize=10)

# Subplot 3: inside TR
ax2_3 = fig2.add_subplot(3, 2, 3)
plot_count_heatmap(fig2, ax2_3,
                   A.iloc[:, 13 + CALLER_ID].values,
                   A.iloc[:, 16 + CALLER_ID].values)
ax2_3.set_xlabel('Number of DEL')
ax2_3.set_ylabel('Number of INS')
ax2_3.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=20bp, inside TR')
ax2_3.set_box_aspect(1)
ax2_3.set_axisbelow(False)
ax2_3.grid(True, color='lightgray')
ax2_3.tick_params(labelsize=10)

# Subplot 5: outside TR
ax2_5 = fig2.add_subplot(3, 2, 5)
plot_count_heatmap(fig2, ax2_5,
                   A.iloc[:, 25 + CALLER_ID].values,
                   A.iloc[:, 28 + CALLER_ID].values)
ax2_5.set_xlabel('Number of DEL')
ax2_5.set_ylabel('Number of INS')
ax2_5.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=20bp, outside TR')
ax2_5.set_box_aspect(1)
ax2_5.set_axisbelow(False)
ax2_5.grid(True, color='lightgray')
ax2_5.tick_params(labelsize=10)

# >=50bp
# Subplot 2: whole genome
ax2_2 = fig2.add_subplot(3, 2, 2)
plot_count_heatmap(fig2, ax2_2,
                   A.iloc[:, 37 + CALLER_ID].values,
                   A.iloc[:, 40 + CALLER_ID].values)
ax2_2.set_xlabel('Number of DEL')
ax2_2.set_ylabel('Number of INS')
ax2_2.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=50bp, whole genome')
ax2_2.set_box_aspect(1)
ax2_2.set_axisbelow(False)
ax2_2.grid(True, color='lightgray')
ax2_2.tick_params(labelsize=10)

# Subplot 4: inside TR
ax2_4 = fig2.add_subplot(3, 2, 4)
plot_count_heatmap(fig2, ax2_4,
                   A.iloc[:, 49 + CALLER_ID].values,
                   A.iloc[:, 52 + CALLER_ID].values)
ax2_4.set_xlabel('Number of DEL')
ax2_4.set_ylabel('Number of INS')
ax2_4.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=50bp, inside TR')
ax2_4.set_box_aspect(1)
ax2_4.set_axisbelow(False)
ax2_4.grid(True, color='lightgray')
ax2_4.tick_params(labelsize=10)

# Subplot 6: outside TR
ax2_6 = fig2.add_subplot(3, 2, 6)
plot_count_heatmap(fig2, ax2_6,
                   A.iloc[:, 61 + CALLER_ID].values,
                   A.iloc[:, 64 + CALLER_ID].values)
ax2_6.set_xlabel('Number of DEL')
ax2_6.set_ylabel('Number of INS')
ax2_6.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=50bp, outside TR')
ax2_6.set_box_aspect(1)
ax2_6.set_axisbelow(False)
ax2_6.grid(True, color='lightgray')
ax2_6.tick_params(labelsize=10)

fig2.tight_layout()
fig2.savefig(f'/Users/fcunial/Downloads/svqc/plot/{COHORT_NAME}/counts_del_ins.png', dpi=300, bbox_inches='tight')
plt.show()
