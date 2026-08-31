from pathlib import Path

import yaml

ROOT = Path(__file__).parents[2]
DEPLOYMENT = ROOT / "deployment"


def _yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text())


def test_production_compose_pins_control_plane_and_isolates_networks() -> None:
    compose = _yaml(DEPLOYMENT / "compose.production.yml")
    services = compose["services"]

    assert services["shinyproxy"]["image"] == "openanalytics/shinyproxy:3.2.4"
    assert services["traefik"]["image"] == "traefik:v3.6"
    assert services["shinyproxy"]["read_only"] is True
    assert services["traefik"]["read_only"] is True
    assert services["shinyproxy"]["cap_drop"] == ["ALL"]
    assert services["traefik"]["cap_drop"] == ["ALL"]
    assert any(
        mount.startswith("/tmp:rw,exec,nosuid,nodev,")
        for mount in services["shinyproxy"]["tmpfs"]
    )
    assert "/var/run/docker.sock:/var/run/docker.sock:ro" in services["shinyproxy"][
        "volumes"
    ]
    assert all("docker.sock" not in volume for volume in services["traefik"]["volumes"])
    assert compose["networks"]["apps"]["internal"] is True
    assert set(services["traefik"]["networks"]) == {"edge"}
    assert set(services["shinyproxy"]["networks"]) == {"edge", "apps", "egress"}


def test_shinyproxy_enforces_authentication_lifecycle_and_resource_limits() -> None:
    config = _yaml(DEPLOYMENT / "shinyproxy" / "application.yml")
    proxy = config["proxy"]
    spec = proxy["specs"][0]

    assert proxy["authentication"] == "simple"
    assert proxy["users"][0]["password"] == "${CLUSTBUSTER_ADMIN_PASSWORD}"
    assert proxy["users"][1]["password"] == "${CLUSTBUSTER_USER_PASSWORD}"
    assert proxy["users"][1]["groups"] == ["scientists"]
    assert proxy["default-stop-proxy-on-logout"] is True
    assert proxy["stop-proxies-on-shutdown"] is True
    assert proxy["recover-running-proxies"] is False
    assert proxy["docker"]["internal-networking"] is True
    assert proxy["docker"]["privileged"] is False
    assert spec["container-network"] == "clustbuster-apps"
    assert spec["container-network-connections"] == ["clustbuster-egress"]
    assert spec["max-instances"] == 1
    assert spec["stop-on-logout"] is True
    assert "CLUSTBUSTER_SESSION_CPUS" in spec["container-cpu-limit"]
    assert "CLUSTBUSTER_SESSION_MEMORY" in spec["container-memory-limit"]
    assert spec["container-env"]["CLUSTBUSTER_ENABLE_SEURAT_IMPORT"] == "1"
    assert spec["container-env"]["CLUSTBUSTER_ENABLE_SCE_IMPORT"] == "1"


def test_traefik_is_tls_only_and_uses_static_service_discovery() -> None:
    static = _yaml(DEPLOYMENT / "traefik" / "traefik.yml")
    dynamic = _yaml(DEPLOYMENT / "traefik" / "dynamic.yml")

    assert "docker" not in static["providers"]
    assert static["entryPoints"]["web"]["http"]["redirections"]["entryPoint"][
        "scheme"
    ] == "https"
    router = dynamic["http"]["routers"]["clustbuster"]
    assert router["entryPoints"] == ["websecure"]
    assert "tls" in router
    assert dynamic["tls"]["options"]["modern"]["minVersion"] == "VersionTLS12"


def test_example_environment_contains_no_password() -> None:
    env = (DEPLOYMENT / ".env.production.example").read_text().splitlines()
    password_lines = [line for line in env if "PASSWORD=" in line]
    assert password_lines == [
        "CLUSTBUSTER_ADMIN_PASSWORD=",
        "CLUSTBUSTER_USER_PASSWORD=",
    ]
