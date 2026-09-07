"""Exercise bundled resources and offline analysis inside a release container."""
from pathlib import Path
from unittest.mock import patch

from clustbuster.app import _default_enrichment_mode
from clustbuster.config import AppConfig
from clustbuster.integrations.enrichr import DEFAULT_LIBRARY
from clustbuster.integrations.offline_enrichment import GeneSetCache, OfflineEnrichmentClient
from clustbuster.models import ReferenceFilters
from clustbuster.resources.providers import (
    marker_provider_from_config,
    reference_provider_from_config,
)

config = AppConfig.from_env()
cache = GeneSetCache(config.gene_set_cache_root)
assert _default_enrichment_mode(cache.root) == 'offline'
assert len(list(cache.root.rglob('*.gmt'))) == 10
with patch('clustbuster.integrations.enrichr.urlopen', side_effect=AssertionError('Network used')):
    for species, genes in [
        ('human', ('APOE', 'LYZ', 'CD3D', 'CD3E', 'MS4A1', 'CD79A', 'NKG7', 'ACTB')),
        ('mouse', ('Apoe', 'Lyz2', 'Cd3d', 'Cd3e', 'Ms4a1', 'Cd79a', 'Nkg7', 'Actb')),
    ]:
        client = OfflineEnrichmentClient(cache, DEFAULT_LIBRARY, genes, 'auto')
        assert client.species == species
        assert not client.enrich(genes[:4], description='container smoke').values.empty
        print(f'{species} offline enrichment passed')
marker = marker_provider_from_config(config)
assert marker.status().available
matches = marker.search_cell_types('', limit=5)
assert matches
reference = reference_provider_from_config(config)
assert reference.status().available
refs = reference.list_references(ReferenceFilters())
assert refs
for ref in refs:
    reference.load_reference(ref.reference_id)
print(f'Atlas marker search and loading all {len(refs)} Refmats passed')
(cache.root / 'persistence-smoke.txt').write_text('persisted')
assert Path('/workspace').is_dir()
