"""Copy a local Atlas snapshot into the container build context with checksums."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import UTC, datetime
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('atlas', type=Path, help='Sovereign Atlas checkout')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
destination = root / 'resources/atlas'
if destination.exists():
    raise SystemExit('resources/atlas already exists; move it aside before staging a new snapshot.')
data = args.atlas / 'data'
required = ['markercodex.duckdb', 'reference_matrices.duckdb', 'reference_matrices/files']
for name in required:
    if not (data / name).exists():
        raise SystemExit(f'Missing Atlas resource: {data / name}')
revision = subprocess.check_output(['git', '-C', str(args.atlas), 'rev-parse', 'HEAD'], text=True).strip()
for name in required:
    source, target = data / name, destination / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, target)
    else:
        shutil.copy2(source, target)
checksums = {
    str(path.relative_to(destination)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in sorted(destination.rglob('*')) if path.is_file()
}
(destination / 'manifest.json').write_text(json.dumps({
    'atlas_revision': revision, 'staged_at': datetime.now(UTC).isoformat(), 'sha256': checksums,
}, indent=2) + '\n')
print(f'Staged {len(checksums)} Atlas files from {revision}')
