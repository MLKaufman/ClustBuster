"""Seed persistent offline resources, then start the requested container command."""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile


def main() -> None:
    bundled = Path('/app/resources/gene_sets')
    cache = Path(os.environ.get('CLUSTBUSTER_GENE_SET_CACHE_ROOT', '/data/gene_sets'))
    cache.mkdir(parents=True, exist_ok=True)
    for source in bundled.rglob('*'):
        if not source.is_file() or source.suffix not in {'.gmt', '.json'}:
            continue
        destination = cache / source.relative_to(bundled)
        if destination.exists():
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        with NamedTemporaryFile(dir=destination.parent, delete=False) as temporary:
            temporary_path = Path(temporary.name)
        try:
            shutil.copyfile(source, temporary_path)
            temporary_path.replace(destination)
        finally:
            temporary_path.unlink(missing_ok=True)

    atlas = Path('/app/resources/atlas')
    if (atlas / 'manifest.json').is_file():
        import json
        manifest = json.loads((atlas / 'manifest.json').read_text())
        defaults = {
            'CLUSTBUSTER_MARKER_DB_PATH': str(atlas / 'markercodex.duckdb'),
            'CLUSTBUSTER_REFERENCE_CATALOG_PATH': str(atlas / 'reference_matrices.duckdb'),
            'CLUSTBUSTER_REFERENCE_FILES_ROOT': str(atlas / 'reference_matrices/files'),
            'CLUSTBUSTER_ATLAS_VERSION': manifest['atlas_revision'],
        }
        for key, value in defaults.items():
            os.environ.setdefault(key, value)
    os.execvp(sys.argv[1], sys.argv[1:])


if __name__ == '__main__':
    main()
