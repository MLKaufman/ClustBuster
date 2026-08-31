import pytest

from clustbuster.config import AppConfig, ConfigurationError


def test_config_reads_namespaced_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CLUSTBUSTER_PORT", "9000")
    monkeypatch.setenv("CLUSTBUSTER_ENRICHMENT_TIMEOUT_SECONDS", "30")
    monkeypatch.setenv("CLUSTBUSTER_ENABLE_SEURAT_IMPORT", "1")
    monkeypatch.setenv("CLUSTBUSTER_ENABLE_SCE_IMPORT", "1")
    config = AppConfig.from_env()
    assert config.port == 9000
    assert config.enrichment_timeout_seconds == 30
    assert config.enable_seurat_import is True
    assert config.enable_sce_import is True


def test_config_rejects_invalid_positive_integers(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CLUSTBUSTER_MAX_UPLOAD_MB", "0")
    with pytest.raises(ConfigurationError, match="greater than zero"):
        AppConfig.from_env()
