#!/usr/bin/env bash
set -euo pipefail

clustbuster_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
atlas_repo="${CLUSTBUSTER_ATLAS_REPO:-${clustbuster_root}/../sovereign-atlas}"
atlas_data="${atlas_repo}/data"

cd "${clustbuster_root}"

if [[ -f "${atlas_data}/markercodex.duckdb" \
   && -f "${atlas_data}/reference_matrices.duckdb" \
   && -d "${atlas_data}/reference_matrices/files" ]]; then
    export CLUSTBUSTER_MARKER_DB_PATH="${atlas_data}/markercodex.duckdb"
    export CLUSTBUSTER_REFERENCE_CATALOG_PATH="${atlas_data}/reference_matrices.duckdb"
    export CLUSTBUSTER_REFERENCE_FILES_ROOT="${atlas_data}/reference_matrices/files"
    if [[ -z "${CLUSTBUSTER_ATLAS_VERSION:-}" ]]; then
        CLUSTBUSTER_ATLAS_VERSION="$(git -C "${atlas_repo}" rev-parse --short HEAD)"
        export CLUSTBUSTER_ATLAS_VERSION
    fi
    echo "Using Sovereign Atlas ${CLUSTBUSTER_ATLAS_VERSION} from ${atlas_repo}"
else
    echo "Sovereign Atlas not found at ${atlas_repo}; using bundled demo resources"
fi

uv sync --extra dev
CLUSTBUSTER_ENABLE_SEURAT_IMPORT=1 \
CLUSTBUSTER_ENABLE_SCE_IMPORT=1 \
uv run clustbuster --reload --launch-browser
