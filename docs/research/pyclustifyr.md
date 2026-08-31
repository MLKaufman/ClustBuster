# pyclustifyr integration research

Checked on 2026-08-30 against
[`MLKaufman/pyclustifyr`](https://github.com/MLKaufman/pyclustifyr), revision
`db8761a87072b814f95ce7e0767540b4d10689bd` (package version `0.1.0`).

## Contract used by ClustBuster

- `clustify(input, ref_mat, ..., per_cell=True, compute_method=..., if_log=False)`
  accepts gene-by-item pandas matrices and returns a query-by-reference score matrix.
- ClustBuster performs sparse-safe pseudobulk averaging first, passes the resulting
  cluster profiles as items, and uses `per_cell=True` to avoid densifying the full
  cell-level matrix inside the package.
- Query and reference matrices are aligned to the same validated gene names before
  the call. The adapter defaults to Spearman similarity and permits Pearson or cosine.
- ClustBuster derives a single preview row per source cluster from the returned score
  matrix. It preserves pyclustifyr-style `unassigned` calls below the explicit cutoff
  and marks tied best calls with a `-CLASH!` suffix.

The contract test runs the bundled AnnData/reference pair through the real pinned
package and asserts all five known cluster labels, the score-matrix shape, package
revision, gene overlap, and that `adata.obs` remains unchanged.

## Current limitations

- The upstream project has no tagged release or PyPI distribution, so ClustBuster
  pins a Git revision instead of a release version.
- The package does not expose `__version__`; the adapter reads installed package
  metadata and separately records the pinned revision.
- ClustBuster currently exposes cluster-level mean pseudobulk scoring only. Per-cell,
  permutation, marker-list, Kendall, and KL-divergence modes are not exposed.
- Reference matching is exact first and then case-insensitive only when unambiguous.
  Missing and ambiguous genes are reported; insufficient overlap stops the run.
- Whether the selected query expression is log-normalized must currently be stated by
  the user. When enabled, sparse values are back-transformed before averaging and the
  cluster means are transformed with `log1p`, matching pyclustifyr's mean-pseudobulk
  behavior.
