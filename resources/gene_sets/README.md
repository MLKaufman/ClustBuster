# Offline gene-set cache

Use Settings → Download for offline use to populate this directory with GMT libraries
from Enrichr. Set CLUSTBUSTER_GENE_SET_CACHE_ROOT to relocate the cache. Each library
has a JSON download record with its source URL, timestamp and SHA-256. GMT/JSON files
are ignored by Git; see the main README for offline methods and staging instructions.

Mouse GO and Reactome GMTs are stored separately in `mouse/` (MSigDB 2025.1.Mm).
Download both species from Settings for automatic species detection. No homolog
conversion is performed. Copy the entire directory, including `mouse/`, to offline hosts.
