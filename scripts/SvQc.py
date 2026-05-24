import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configuration
CALLER_ID = 0  # 0=pav, 1=pbsv, 2=sniffles

SV_TYPES_PRIMARY = ['DEL', 'INS', 'DUP', 'INV']
SV_TYPES_SECONDARY = ['BND', 'SUB', 'UNK']
MAX_COVERAGE = 50
MAX_NRECORDS = 40000
CALLER_NAMES = ['PAV', 'PBSV', 'Sniffles']

# Load data
A = pd.read_csv('/Users/fcunial/Downloads/svqc/plot/counts.csv', header=None)

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

# >=20bp
# Subplot 1: whole genome
ax1 = plt.subplot(3, 3, 1)
for i in [1, 4, 7, 10]:
    ax1.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax1.set_xlabel('Coverage')
ax1.set_ylabel('Number of records')
ax1.set_title(f'{CALLER_NAMES[CALLER_ID]}, >=20bp, whole genome')
ax1.set_xlim(0, MAX_COVERAGE)
ax1.set_ylim(0, MAX_NRECORDS)
ax1.set_box_aspect(1)
ax1.grid(True, color='lightgray')
ax1.tick_params(labelsize=10)
ax1.legend(['DEL', 'INS', 'DUP', 'INV'], loc='upper left')

# Subplot 4: inside TR
ax4 = plt.subplot(3, 3, 4)
for i in [13, 16, 19, 22]:
    ax4.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax4.set_xlabel('Coverage')
ax4.set_ylabel('Number of records')
ax4.set_title(f'{CALLER_NAMES[CALLER_ID]}, >=20bp, inside TR')
ax4.set_xlim(0, MAX_COVERAGE)
ax4.set_ylim(0, MAX_NRECORDS)
ax4.set_box_aspect(1)
ax4.grid(True, color='lightgray')
ax4.tick_params(labelsize=10)
ax4.legend(['DEL', 'INS', 'DUP', 'INV'], loc='upper left')

# Subplot 7: outside TR
ax7 = plt.subplot(3, 3, 7)
for i in [25, 28, 31, 34]:
    ax7.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax7.set_xlabel('Coverage')
ax7.set_ylabel('Number of records')
ax7.set_title(f'{CALLER_NAMES[CALLER_ID]}, >=20bp, outside TR')
ax7.set_xlim(0, MAX_COVERAGE)
ax7.set_ylim(0, MAX_NRECORDS)
ax7.set_box_aspect(1)
ax7.grid(True, color='lightgray')
ax7.tick_params(labelsize=10)
ax7.legend(['DEL', 'INS', 'DUP', 'INV'], loc='upper left')

# >=50bp
# Subplot 2: whole genome
ax2 = plt.subplot(3, 3, 2)
for i in [37, 40, 43, 46]:
    ax2.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax2.set_xlabel('Coverage')
ax2.set_ylabel('Number of records')
ax2.set_title(f'{CALLER_NAMES[CALLER_ID]}, >=50bp, whole genome')
ax2.set_xlim(0, MAX_COVERAGE)
ax2.set_ylim(0, MAX_NRECORDS)
ax2.set_box_aspect(1)
ax2.grid(True, color='lightgray')
ax2.tick_params(labelsize=10)
ax2.legend(['DEL', 'INS', 'DUP', 'INV'], loc='upper left')

# Subplot 5: inside TR
ax5 = plt.subplot(3, 3, 5)
for i in [49, 52, 55, 58]:
    ax5.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax5.set_xlabel('Coverage')
ax5.set_ylabel('Number of records')
ax5.set_title(f'{CALLER_NAMES[CALLER_ID]}, >=50bp, inside TR')
ax5.set_xlim(0, MAX_COVERAGE)
ax5.set_ylim(0, MAX_NRECORDS)
ax5.set_box_aspect(1)
ax5.grid(True, color='lightgray')
ax5.tick_params(labelsize=10)
ax5.legend(['DEL', 'INS', 'DUP', 'INV'], loc='upper left')

# Subplot 8: outside TR
ax8 = plt.subplot(3, 3, 8)
for i in [61, 64, 67, 70]:
    ax8.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax8.set_xlabel('Coverage')
ax8.set_ylabel('Number of records')
ax8.set_title(f'{CALLER_NAMES[CALLER_ID]}, >=50bp, outside TR')
ax8.set_xlim(0, MAX_COVERAGE)
ax8.set_ylim(0, MAX_NRECORDS)
ax8.set_box_aspect(1)
ax8.grid(True, color='lightgray')
ax8.tick_params(labelsize=10)
ax8.legend(['DEL', 'INS', 'DUP', 'INV'], loc='upper left')

# Other types
# Subplot 3: whole genome
ax3 = plt.subplot(3, 3, 3)
for i in [73, 76, 79]:
    ax3.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax3.set_xlabel('Coverage')
ax3.set_ylabel('Number of records')
ax3.set_title(f'{CALLER_NAMES[CALLER_ID]}, Other types, whole genome')
ax3.set_xlim(0, MAX_COVERAGE)
ax3.set_ylim(0, MAX_NRECORDS)
ax3.set_box_aspect(1)
ax3.grid(True, color='lightgray')
ax3.tick_params(labelsize=10)
ax3.legend(['BND', 'SUB', 'UNK'], loc='upper left')

# Subplot 6: inside TR
ax6 = plt.subplot(3, 3, 6)
for i in [82, 85, 88]:
    ax6.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax6.set_xlabel('Coverage')
ax6.set_ylabel('Number of records')
ax6.set_title(f'{CALLER_NAMES[CALLER_ID]}, Other types, inside TR')
ax6.set_xlim(0, MAX_COVERAGE)
ax6.set_ylim(0, MAX_NRECORDS)
ax6.set_box_aspect(1)
ax6.grid(True, color='lightgray')
ax6.tick_params(labelsize=10)
ax6.legend(['BND', 'SUB', 'UNK'], loc='upper left')

# Subplot 9: outside TR
ax9 = plt.subplot(3, 3, 9)
for i in [91, 94, 97]:
    ax9.plot(x, A.iloc[:, 1 + i + CALLER_ID], '.')
ax9.set_xlabel('Coverage')
ax9.set_ylabel('Number of records')
ax9.set_title(f'{CALLER_NAMES[CALLER_ID]}, Other types, outside TR')
ax9.set_xlim(0, MAX_COVERAGE)
ax9.set_ylim(0, MAX_NRECORDS)
ax9.set_box_aspect(1)
ax9.grid(True, color='lightgray')
ax9.tick_params(labelsize=10)
ax9.legend(['BND', 'SUB', 'UNK'], loc='upper left')

plt.tight_layout()
plt.show()
