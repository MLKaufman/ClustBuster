"""Generate the small redistributable H5AD used for manual ClustBuster testing."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main() -> None:
    rng = np.random.default_rng(20260830)
    genes = ["CD3D", "IL7R", "LYZ", "S100A8", "MS4A1", "CD79A", "NKG7", "GNLY", "PPBP", "FCGR3A"]
    cluster_names = ["T cells", "Myeloid", "B cells", "NK cells", "Platelets"]
    cells_per_cluster = 18
    clusters = np.repeat(cluster_names, cells_per_cluster)
    cell_count = len(clusters)

    counts = rng.poisson(0.18, size=(cell_count, len(genes))).astype(np.int32)
    marker_pairs = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
    for cluster_index, marker_indices in enumerate(marker_pairs):
        start = cluster_index * cells_per_cluster
        stop = start + cells_per_cluster
        counts[start:stop, marker_indices] += rng.poisson(
            4.5, size=(cells_per_cluster, len(marker_indices))
        )

    angles = np.linspace(0, 2 * np.pi, len(cluster_names), endpoint=False)
    centers = np.column_stack([4 * np.cos(angles), 4 * np.sin(angles)])
    embedding = np.vstack(
        [
            centers[index] + rng.normal(0, 0.42, size=(cells_per_cluster, 2))
            for index in range(len(cluster_names))
        ]
    )
    pca = rng.normal(size=(cell_count, 6))
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical(np.repeat(["0", "1", "2", "3", "4"], cells_per_cluster)),
            "cell_type_hint": pd.Categorical(clusters),
            "sample": pd.Categorical(np.resize(["donor-A", "donor-B", "donor-C"], cell_count)),
            "condition": pd.Categorical(np.resize(["control", "treated"], cell_count)),
        },
        index=[f"demo-cell-{index:03d}" for index in range(cell_count)],
    )
    var = pd.DataFrame(
        {"gene_symbol": genes, "feature_type": pd.Categorical(["Gene Expression"] * len(genes))},
        index=genes,
    )
    count_matrix = sparse.csr_matrix(counts)
    normalized = sparse.csr_matrix(np.log1p(counts).astype(np.float32))
    adata = ad.AnnData(X=normalized, obs=obs, var=var)
    adata.layers["counts"] = count_matrix
    adata.layers["log1p"] = normalized.copy()
    adata.raw = ad.AnnData(X=count_matrix.copy(), obs=obs.copy(), var=var.copy())
    adata.obsm["X_umap"] = embedding
    adata.obsm["X_pca"] = pca
    adata.uns["fixture"] = {
        "generator": "scripts/create_test_fixtures/generate_h5ad.py",
        "seed": 20260830,
        "purpose": "ClustBuster manual and integration testing",
    }

    destination = Path(__file__).resolve().parents[2] / "testdata" / "clustbuster-demo.h5ad"
    destination.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(destination, compression="gzip")
    print(f"Wrote {destination} ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()

