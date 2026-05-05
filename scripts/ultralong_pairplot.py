# a bit of code for reading the HDF5 file of the extracted annotations and generating a pairplot of
# a subsample of 1000 SNPs and 1000 indels

import h5py
import numpy as np
import seaborn as sns
import pandas as pd

# advanced users can use this code as a rough guide to understand the schema of the HDF5 file,
# which contains labels, intervals, etc. that aren't extracted here.
# see also https://www.hdfgroup.org/downloads/hdfview/ for a GUI tool that can be used to inspect the HDF5
def read_annotations(h5file):
    with h5py.File(h5file, 'r') as f:
        annotation_names_i = f['/annotations/names'][()].astype(str)

        # read chunked annotations
        num_chunks = int(f['/annotations/num_chunks'][()])
        num_columns = int(f['/annotations/num_columns'][()])
        num_rows = int(f['/annotations/num_rows'][()])
        X_ni = np.zeros((num_rows, num_columns))
        n = 0
        for chunk_index in range(num_chunks):
            chunk_ni = f[f'/annotations/chunk_{chunk_index}'][()]
            num_rows_in_chunk = len(chunk_ni)
            X_ni[n:n + num_rows_in_chunk, :] = chunk_ni
            n += num_rows_in_chunk
        snp_n = f['/labels/snp'][()].astype(bool)
        assert n == num_rows
        
    return annotation_names_i, X_ni, snp_n

annotation_names_i, X_ni, snp_n = read_annotations(f'{TUTORIAL_DIR}/vets.extract.annot.hdf5')
annotations_df = pd.DataFrame(X_ni, columns=annotation_names_i)
annotations_df['type'] = np.take(['indel', 'snp'], snp_n)


wanted = ["LRCALLER_AD11_left", "LRCALLER_AD12_left", "LRCALLER_AD13_left", "LRCALLER_GTCOUNT1_left", "LRCALLER_PL11_left", "LRCALLER_PL12_left", "LRCALLER_PL13_left", "LRCALLER_VA11_left", "LRCALLER_VA12_left", "LRCALLER_VA13_left", "SVLEN" ]
selected = [c for c in wanted if c in annotations_df.columns]
plot_df = pd.concat(
    [annotations_df[snp_n].sample(400, random_state=1),
     annotations_df[~snp_n].sample(400, random_state=1)],
    ignore_index=True
)
sns.pairplot(
    plot_df,
    vars=selected,
    hue="type",
    hue_order=["indel", "snp"],
    corner=True,
    diag_kind="kde", 
    plot_kws={"s": 3, "alpha": 0.5}
)