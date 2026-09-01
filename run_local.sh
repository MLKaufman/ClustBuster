#!/usr/bin/env bash
set -euo pipefail

uv sync --extra dev
CLUSTBUSTER_ENABLE_SEURAT_IMPORT=1 \
CLUSTBUSTER_ENABLE_SCE_IMPORT=1 \
uv run clustbuster --reload --launch-browser
