import pytest

from clustbuster.config import AppConfig, ConfigurationError


def test_config_reads_namespaced_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CLUSTBUSTER_PORT", "9000")
    monkeypatch.setenv("CLUSTBUSTER_ENRICHMENT_TIMEOUT_SECONDS", "30")
    monkeypatch.setenv("CLUSTBUSTER_ENABLE_SEURAT_IMPORT", "1")
    monkeypatch.setenv("CLUSTBUSTER_ENABLE_SCE_IMPORT", "1")
    monkeypatch.setenv("CLUSTBUSTER_MARKER_DB_PATH", "/atlas/markers.duckdb")
    monkeypatch.setenv("CLUSTBUSTER_REFERENCE_CATALOG_PATH", "/atlas/references.duckdb")
    monkeypatch.setenv("CLUSTBUSTER_REFERENCE_FILES_ROOT", "/atlas/files")
    monkeypatch.setenv("CLUSTBUSTER_ATLAS_VERSION", "preview-7")
    config = AppConfig.from_env()
    assert config.port == 9000
    assert config.enrichment_timeout_seconds == 30
    assert config.enable_seurat_import is True
    assert config.enable_sce_import is True
    assert str(config.marker_db_path) == "/atlas/markers.duckdb"
    assert str(config.reference_catalog_path) == "/atlas/references.duckdb"
    assert str(config.reference_files_root) == "/atlas/files"
    assert config.atlas_version == "preview-7"


def test_config_rejects_invalid_positive_integers(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CLUSTBUSTER_MAX_UPLOAD_MB", "0")
    with pytest.raises(ConfigurationError, match="greater than zero"):
        AppConfig.from_env()
