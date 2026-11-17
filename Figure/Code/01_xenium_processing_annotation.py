#!/usr/bin/env python3

import numpy as np
import scanpy as sc
import stlearn as st

# ---- 1. QC ------------------------------------------------------------------
st.pp.filter_genes(adata, min_counts=10)   # genes in >=10 cells
st.pp.filter_cells(adata, min_counts=10)   # cells with >=10 transcripts

# ---- 2. Normalization + Freeman-Tukey variance stabilization ----------------
st.pp.normalize_total(adata)
adata.X = np.sqrt(adata.X.toarray()) + np.sqrt(adata.X.toarray() + 1)

# ---- 3. Embedding and clustering ------------------------
st.em.run_pca(adata, n_comps=50, random_state=0)
st.pp.neighbors(adata, n_neighbors=25, use_rep="X_pca", random_state=0)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
sc.tl.umap(adata)

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
adata.obs["anno"] = adata.obs["leiden"].replace({...})