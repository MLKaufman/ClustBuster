"""Construct biological-resource providers from centralized configuration."""

from clustbuster.config import AppConfig
from clustbuster.resources.base import MarkerProvider, ReferenceProvider
from clustbuster.resources.markers.csv import CsvMarkerProvider
from clustbuster.resources.markers.duckdb import DuckDbMarkerProvider
from clustbuster.resources.references.duckdb import DuckDbReferenceProvider
from clustbuster.resources.references.local import LocalReferenceProvider


def marker_provider_from_config(config: AppConfig) -> MarkerProvider:
    if config.marker_db_path is not None:
        return DuckDbMarkerProvider(
            config.marker_db_path,
            resource_version=config.atlas_version,
        )
    return CsvMarkerProvider(config.marker_catalog_path)


def reference_provider_from_config(config: AppConfig) -> ReferenceProvider:
    if config.reference_catalog_path is not None:
        files_root = config.reference_files_root
        if files_root is None:
            files_root = config.reference_catalog_path.parent / "reference_matrices" / "files"
        return DuckDbReferenceProvider(
            config.reference_catalog_path,
            files_root,
            resource_version=config.atlas_version,
        )
    return LocalReferenceProvider(config.reference_root)
