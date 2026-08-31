"""Centralized application configuration."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


class ConfigurationError(ValueError):
    """Raised when startup configuration is invalid."""


def _positive_int(value: str, name: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ConfigurationError(f"{name} must be an integer") from exc
    if parsed <= 0:
        raise ConfigurationError(f"{name} must be greater than zero")
    return parsed


@dataclass(frozen=True, slots=True)
class AppConfig:
    """Validated configuration populated once at application startup."""

    host: str = "127.0.0.1"
    port: int = 8000
    workspace_root: Path = Path("/tmp/clustbuster")
    max_upload_mb: int = 2048
    cache_limit_mb: int = 1024
    enrichment_timeout_seconds: int = 20
    log_level: str = "INFO"
    enable_seurat_import: bool = False
    enable_sce_import: bool = False

    @classmethod
    def from_env(cls) -> AppConfig:
        prefix = "CLUSTBUSTER_"
        log_level = os.getenv(f"{prefix}LOG_LEVEL", "INFO").upper()
        if log_level not in {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}:
            raise ConfigurationError("CLUSTBUSTER_LOG_LEVEL is not a valid log level")
        return cls(
            host=os.getenv(f"{prefix}HOST", "127.0.0.1"),
            port=_positive_int(os.getenv(f"{prefix}PORT", "8000"), f"{prefix}PORT"),
            workspace_root=Path(
                os.getenv(f"{prefix}WORKSPACE_ROOT", "/tmp/clustbuster")
            ).expanduser(),
            max_upload_mb=_positive_int(
                os.getenv(f"{prefix}MAX_UPLOAD_MB", "2048"), f"{prefix}MAX_UPLOAD_MB"
            ),
            cache_limit_mb=_positive_int(
                os.getenv(f"{prefix}CACHE_LIMIT_MB", "1024"), f"{prefix}CACHE_LIMIT_MB"
            ),
            enrichment_timeout_seconds=_positive_int(
                os.getenv(f"{prefix}ENRICHMENT_TIMEOUT_SECONDS", "20"),
                f"{prefix}ENRICHMENT_TIMEOUT_SECONDS",
            ),
            log_level=log_level,
            enable_seurat_import=os.getenv(f"{prefix}ENABLE_SEURAT_IMPORT", "0") == "1",
            enable_sce_import=os.getenv(f"{prefix}ENABLE_SCE_IMPORT", "0") == "1",
        )
