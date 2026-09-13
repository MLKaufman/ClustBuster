# Reference metadata for ClustBuster

ClustBuster reads reference catalogs without modifying them. Existing catalogs remain
compatible: new metadata columns are optional. Search covers the reference title,
species, tissue, condition, assay, cell-type column names, and optional metadata,
including descriptions, submitters, and per-cell-type aliases. Searches ignore case
and normalize punctuation, so `Luminal.progenitor` matches `luminal progenitor`.
Dropdown choices show title, species, tissue, condition, assay, cell-type count, and
an ID prefix. The full stable ID is available in the details panel.

## Populate existing fields first

Supply `description`, `species_scientific_name`, `species_common_name`, `tissue`,
`condition`, `developmental_stage`, `assay`, `platform`, `feature_id_type`,
`normalization`, `value_type`, `source_title`, `citation`, `doi`, `pmid`, `source_url`,
`publication_year`, `data_license`, `submitter`, and `notes` where known.
Keep missing values NULL, rather than inventing values. For unpublished work,
identify the source study and record that it is unpublished in notes.

For the current mammary reference, verify the scientific name (`Mus musculus`),
normalization, gene identifier type, platform, provenance, and usage permissions.
The application cannot infer these reliably from an expression matrix alone.

Keep `column_names` synchronized with the matrix header, including its `gene`
column. ClustBuster excludes `gene` when listing cell types. Maintain the actual
row count, column count, checksum, and file metadata when uploading a new version.
Search uses catalog column names; selection validates the actual matrix/checksum.

## Optional additional database fields

These fields are now understood by the ClustBuster provider:

| Field | Suggested type | Meaning |
|---|---|---|
| taxonomy_id | INTEGER | NCBI taxonomy identifier, e.g. 10090 or 9606 |
| reference_version | VARCHAR | Version of this reference, separate from catalog release |
| construction_method | VARCHAR | Averaging/pseudobulk method, filtering and replicate handling |
| intended_use | VARCHAR | Appropriate tissues, conditions and annotation resolution |
| limitations | VARCHAR | Missing populations, confounders and known limitations |
| cell_count | BIGINT | Total source cells |
| sample_count | BIGINT | Independent source samples |
| donor_count | BIGINT | Independent donors/animals |
| cell_type_metadata | JSON or VARCHAR | Metadata keyed by exact matrix column name |

Example `cell_type_metadata` value (illustrative, not real reference data):

```json
{
  "Fibroblasts": {
    "cell_count": 1200,
    "sample_count": 6,
    "donor_count": 3,
    "description": "Description of the reference population",
    "ontology_id": "CL:0000057",
    "aliases": ["fibroblast", "stromal fibroblast"]
  }
}
```

The cell-type coverage table shows these fields when present; absent optional
columns are omitted. Aliases and descriptions are searchable. Use exact matrix
column names as keys; an ontology label does not replace the matrix column identity.

For a normalized database, maintain a separate `reference_cell_types` table keyed
by `(matrix_id, cell_type)` and expose the JSON object through the catalog view or
export step. This gives ClustBuster a stable read-only contract without coupling it
to the curator's internal schema.

## Local JSON references

The same optional fields are supported in local JSON sidecars. Include `cell_types`
as an array of matrix column names. Older CSV/TSV/TXT references can discover these
from the header; Parquet sidecars should provide `cell_types` explicitly for search.

## Compatibility reporting

Selecting a reference validates its file in a background task. With a loaded
dataset, the report shows shared, missing and ambiguous genes using exactly the
same matching function as pyclustifyr annotation: exact matches first, otherwise a
unique case-insensitive match. Multiple possible matches are omitted. The searchable
gene report lists each reference gene and its status.

Human/mouse dataset species detection uses locally downloaded enrichment libraries.
If unavailable or ambiguous, species is reported as uncertain. Species mismatch is
reported when recognized. Gene-symbol overlap does not establish biological
compatibility and does not perform ortholog conversion. The selected expression
source's genes, including raw when selected, determine the overlap report.

Changing reference or workspace configuration clears the old prediction preview.
Background results from an older selection cannot replace the current selection.


## Upload your own reference in ClustBuster

On Refmats, select **Upload your own**, choose a CSV or tab-delimited TSV/TXT file,
optionally enter a name, species, tissue and normalization description, then click
**Validate uploaded reference**. The first column must contain unique gene identifiers;
the remaining headers must be unique cell-type names. An unnamed first column, as
in an R export with row names, is accepted. All expression values must be numeric,
finite, and non-missing. Matrices with cell types in rows must be transposed first.

The validated upload appears in the reference selector and uses the same cell-type
coverage, compatibility report, annotation preview, correlation plot and Apply
workflow as a catalog reference. Run annotation after validation; Apply predictions
updates your workspace only when explicitly clicked. Set the minimum shared-gene
threshold as appropriate for your reference.

Uploads are private to the session and never added to the shared Refmat database.
They are deleted with the session workspace. Changing the file clears the old
selection and predictions; validate the replacement before annotating. Metadata is
captured when Validate is clicked. Switch back to **Reference catalog** to use the
installed references. Upload mode also works when the catalog is unavailable.
