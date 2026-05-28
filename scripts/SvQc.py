import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np


DEFAULT_COHORT_NAME = 'BI'
DEFAULT_CALLER_ID = 1  # 0=pav, 1=pbsv, 2=sniffles
DEFAULT_SELECTED_SV_TYPE = 'DEL'  # DEL, INS, DUP, INV, BND, SUB, UNK
DEFAULT_SUFFIX = 'revio'

SV_TYPES_PRIMARY = ['DEL', 'INS', 'DUP', 'INSDUP', 'INV']
SV_TYPES_SECONDARY = ['BND', 'SUB', 'UNK']
CALLER_NAMES = ['PAV', 'PBSV', 'Sniffles']
HEATMAP_BINS_X = 50
HEATMAP_BINS_Y = 50


def parse_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', default=DEFAULT_COHORT_NAME)
    parser.add_argument(
        '--caller',
        type=int,
        default=DEFAULT_CALLER_ID,
        choices=range(len(CALLER_NAMES)),
        help='0=pav, 1=pbsv, 2=sniffles',
    )
    parser.add_argument(
        '--svtype',
        default=DEFAULT_SELECTED_SV_TYPE,
        choices=SV_TYPES_PRIMARY + SV_TYPES_SECONDARY,
    )
    parser.add_argument('--suffix', default=DEFAULT_SUFFIX)
    parser.add_argument('--input-dir')
    parser.add_argument(
        '--row',
        default=None,
        help='Comma-separated matrix row to highlight on each plot.',
    )
    args = parser.parse_args()
    input_dir = args.input_dir
    if not input_dir.endswith('/'):
        input_dir = f'{input_dir}/'
    return (
        args.cohort,
        args.caller,
        args.svtype,
        args.suffix,
        input_dir,
        args.row,
    )


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


def parse_highlight_row(row_str, expected_columns):
    if row_str is None:
        return None

    tokens = [token.strip() for token in row_str.split(',')]
    if len(tokens) != expected_columns:
        raise ValueError(
            f'--row has {len(tokens)} values, expected {expected_columns} values.'
        )

    row_values = np.array([float(token) for token in tokens], dtype=float)
    return row_values


def add_highlight_marker(ax_obj, row_values, x_column, y_column):
    if row_values is None:
        return
    x_value = row_values[x_column]
    y_value = row_values[y_column]
    if np.isfinite(x_value) and np.isfinite(y_value):
        ax_obj.plot(
            x_value,
            y_value,
            marker='o',
            markersize=8,
            markerfacecolor='none',
            markeredgecolor='red',
            linestyle='None',
            zorder=10,
        )




# -------------------------------- Main program --------------------------------

COHORT_NAME, CALLER_ID, SELECTED_SV_TYPE, SUFFIX, INPUT_DIR, HIGHLIGHT_ROW_STR = parse_command_line_args()
A = pd.read_csv(f'{INPUT_DIR}counts_{SUFFIX}.csv', header=None)
HIGHLIGHT_ROW = parse_highlight_row(HIGHLIGHT_ROW_STR, A.shape[1])
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




# Figure 1: heatmaps with X=coverage, Y=number of calls of the selected SV type, 
# cell = number of samples.
fig = plt.figure(figsize=(15, 12))
x = A.iloc[:, 0].values  # First column is coverage

# >=20bp
# Subplot 1: whole genome
ax1 = plt.subplot(3, 3, 1)
y1 = np.full(x.shape, np.nan)
y1_column = None
if SELECTED_SV_TYPE == 'DEL':
    y1_column = 1 + CALLER_ID
if SELECTED_SV_TYPE == 'INS':
    y1_column = 4 + CALLER_ID
if SELECTED_SV_TYPE == 'DUP':
    y1_column = 7 + CALLER_ID
if SELECTED_SV_TYPE == 'INV':
    y1_column = 10 + CALLER_ID
if y1_column is not None:
    y1 = A.iloc[:, y1_column].values
plot_count_heatmap(fig, ax1, x, y1)
if y1_column is not None:
    add_highlight_marker(ax1, HIGHLIGHT_ROW, 0, y1_column)
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
y4_column = None
if SELECTED_SV_TYPE == 'DEL':
    y4_column = 13 + CALLER_ID
if SELECTED_SV_TYPE == 'INS':
    y4_column = 16 + CALLER_ID
if SELECTED_SV_TYPE == 'DUP':
    y4_column = 19 + CALLER_ID
if SELECTED_SV_TYPE == 'INV':
    y4_column = 22 + CALLER_ID
if y4_column is not None:
    y4 = A.iloc[:, y4_column].values
plot_count_heatmap(fig, ax4, x, y4)
if y4_column is not None:
    add_highlight_marker(ax4, HIGHLIGHT_ROW, 0, y4_column)
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
y7_column = None
if SELECTED_SV_TYPE == 'DEL':
    y7_column = 25 + CALLER_ID
if SELECTED_SV_TYPE == 'INS':
    y7_column = 28 + CALLER_ID
if SELECTED_SV_TYPE == 'DUP':
    y7_column = 31 + CALLER_ID
if SELECTED_SV_TYPE == 'INV':
    y7_column = 34 + CALLER_ID
if y7_column is not None:
    y7 = A.iloc[:, y7_column].values
plot_count_heatmap(fig, ax7, x, y7)
if y7_column is not None:
    add_highlight_marker(ax7, HIGHLIGHT_ROW, 0, y7_column)
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
y2_column = None
if SELECTED_SV_TYPE == 'DEL':
    y2_column = 37 + CALLER_ID
if SELECTED_SV_TYPE == 'INS':
    y2_column = 40 + CALLER_ID
if SELECTED_SV_TYPE == 'DUP':
    y2_column = 43 + CALLER_ID
if SELECTED_SV_TYPE == 'INV':
    y2_column = 46 + CALLER_ID
if y2_column is not None:
    y2 = A.iloc[:, y2_column].values
plot_count_heatmap(fig, ax2, x, y2)
if y2_column is not None:
    add_highlight_marker(ax2, HIGHLIGHT_ROW, 0, y2_column)
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
y5_column = None
if SELECTED_SV_TYPE == 'DEL':
    y5_column = 49 + CALLER_ID
if SELECTED_SV_TYPE == 'INS':
    y5_column = 52 + CALLER_ID
if SELECTED_SV_TYPE == 'DUP':
    y5_column = 55 + CALLER_ID
if SELECTED_SV_TYPE == 'INV':
    y5_column = 58 + CALLER_ID
if y5_column is not None:
    y5 = A.iloc[:, y5_column].values
plot_count_heatmap(fig, ax5, x, y5)
if y5_column is not None:
    add_highlight_marker(ax5, HIGHLIGHT_ROW, 0, y5_column)
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
y8_column = None
if SELECTED_SV_TYPE == 'DEL':
    y8_column = 61 + CALLER_ID
if SELECTED_SV_TYPE == 'INS':
    y8_column = 64 + CALLER_ID
if SELECTED_SV_TYPE == 'DUP':
    y8_column = 67 + CALLER_ID
if SELECTED_SV_TYPE == 'INV':
    y8_column = 70 + CALLER_ID
if y8_column is not None:
    y8 = A.iloc[:, y8_column].values
plot_count_heatmap(fig, ax8, x, y8)
if y8_column is not None:
    add_highlight_marker(ax8, HIGHLIGHT_ROW, 0, y8_column)
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
y3_column = None
if SELECTED_SV_TYPE == 'BND':
    y3_column = 73 + CALLER_ID
if SELECTED_SV_TYPE == 'SUB':
    y3_column = 76 + CALLER_ID
if SELECTED_SV_TYPE == 'UNK':
    y3_column = 79 + CALLER_ID
if y3_column is not None:
    y3 = A.iloc[:, y3_column].values
plot_count_heatmap(fig, ax3, x, y3)
if y3_column is not None:
    add_highlight_marker(ax3, HIGHLIGHT_ROW, 0, y3_column)
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
y6_column = None
if SELECTED_SV_TYPE == 'BND':
    y6_column = 82 + CALLER_ID
if SELECTED_SV_TYPE == 'SUB':
    y6_column = 85 + CALLER_ID
if SELECTED_SV_TYPE == 'UNK':
    y6_column = 88 + CALLER_ID
if y6_column is not None:
    y6 = A.iloc[:, y6_column].values
plot_count_heatmap(fig, ax6, x, y6)
if y6_column is not None:
    add_highlight_marker(ax6, HIGHLIGHT_ROW, 0, y6_column)
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
y9_column = None
if SELECTED_SV_TYPE == 'BND':
    y9_column = 91 + CALLER_ID
if SELECTED_SV_TYPE == 'SUB':
    y9_column = 94 + CALLER_ID
if SELECTED_SV_TYPE == 'UNK':
    y9_column = 97 + CALLER_ID
if y9_column is not None:
    y9 = A.iloc[:, y9_column].values
plot_count_heatmap(fig, ax9, x, y9)
if y9_column is not None:
    add_highlight_marker(ax9, HIGHLIGHT_ROW, 0, y9_column)
ax9.set_xlabel('Coverage')
ax9.set_ylabel('Number of calls')
ax9.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, {SELECTED_SV_TYPE}, outside TR')
ax9.set_box_aspect(1)
ax9.set_axisbelow(False)
ax9.grid(True, color='lightgray')
ax9.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig(f'{INPUT_DIR}counts_{SUFFIX}_{CALLER_NAMES[CALLER_ID]}_{SELECTED_SV_TYPE}.png', dpi=300, bbox_inches='tight')




# Figure 2: heatmaps with X=number of DEL, Y=number of INS, cell = number of 
# samples.
fig2 = plt.figure(figsize=(10, 12))


# >=20bp
# Subplot 1: whole genome
ax2_1 = fig2.add_subplot(3, 2, 1)
plot_count_heatmap(fig2, ax2_1,
                   A.iloc[:, 1 + CALLER_ID].values,
                   A.iloc[:, 4 + CALLER_ID].values)
add_highlight_marker(ax2_1, HIGHLIGHT_ROW, 1 + CALLER_ID, 4 + CALLER_ID)
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
add_highlight_marker(ax2_3, HIGHLIGHT_ROW, 13 + CALLER_ID, 16 + CALLER_ID)
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
add_highlight_marker(ax2_5, HIGHLIGHT_ROW, 25 + CALLER_ID, 28 + CALLER_ID)
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
add_highlight_marker(ax2_2, HIGHLIGHT_ROW, 37 + CALLER_ID, 40 + CALLER_ID)
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
add_highlight_marker(ax2_4, HIGHLIGHT_ROW, 49 + CALLER_ID, 52 + CALLER_ID)
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
add_highlight_marker(ax2_6, HIGHLIGHT_ROW, 61 + CALLER_ID, 64 + CALLER_ID)
ax2_6.set_xlabel('Number of DEL')
ax2_6.set_ylabel('Number of INS')
ax2_6.set_title(f'{COHORT_NAME}, {CALLER_NAMES[CALLER_ID]}, >=50bp, outside TR')
ax2_6.set_box_aspect(1)
ax2_6.set_axisbelow(False)
ax2_6.grid(True, color='lightgray')
ax2_6.tick_params(labelsize=10)

fig2.tight_layout()
fig2.savefig(f'{INPUT_DIR}counts_{SUFFIX}_{CALLER_NAMES[CALLER_ID]}_DEL_INS.png', dpi=300, bbox_inches='tight')
