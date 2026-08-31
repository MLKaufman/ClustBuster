"""Generate a deterministic H5Seurat compatibility fixture."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from readseurat.convert.seurat.convert import convert_to_h5seurat
from scipy import sparse


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    output = root / "testdata" / "so.h5seurat"
    cells = [f"h5-cell-{index:03d}" for index in range(12)]
    genes = ["CD3D", "IL7R", "LYZ", "MS4A1", "NKG7"]
    counts = sparse.csr_matrix(
        np.arange(len(cells) * len(genes), dtype=np.float32).reshape(len(cells), len(genes))
        % 7
    )
    adata = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(
            {"seurat_clusters": pd.Categorical(np.repeat(["0", "1", "2"], 4))},
            index=cells,
        ),
        var=pd.DataFrame(index=genes),
    )
    adata.layers["counts"] = counts.copy()
    coordinates = np.arange(len(cells), dtype=float)
    adata.obsm["X_umap"] = np.column_stack([coordinates, -coordinates])

    output.parent.mkdir(parents=True, exist_ok=True)
    convert_to_h5seurat(
        adata,
        str(output),
        assay_name="RNA",
        project="ClustBusterH5Fixture",
        version="5.0.0",
    )
    print(f"Wrote {output.relative_to(root)}")


if __name__ == "__main__":
    main()
