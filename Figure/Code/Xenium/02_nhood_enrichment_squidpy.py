#!/usr/bin/env python3

import pandas as pd
import squidpy as sq

rows = []
for sample in adata.obs["sample"].unique():
    a = adata[adata.obs["sample"] == sample].copy()

    # Delaunay spatial graph + permutation-based neighborhood enrichment
    sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)
    sq.gr.nhood_enrichment(a, cluster_key="anno")

    z = a.uns["anno_nhood_enrichment"]["zscore"]
    cats = list(a.obs["anno"].cat.categories)
    for i, c1 in enumerate(cats):
        for j, c2 in enumerate(cats):
            rows.append({"sample": sample, "cluster_1": c1,
                         "cluster_2": c2, "zscore": z[i, j]})

pd.DataFrame(rows).to_csv("nhood_enrichment_all_samples.csv", index=False)
