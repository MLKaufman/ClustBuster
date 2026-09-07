# ClustBuster operations guide

The production stack places Traefik in front of ShinyProxy. Traefik terminates TLS
and forwards HTTP/WebSocket traffic; ShinyProxy authenticates users and creates one
ephemeral ClustBuster container per active user session.

## Prerequisites

- Docker Engine with Compose v2
- a Linux Docker socket group ID (`stat -c '%g' /var/run/docker.sock`)
- a trusted TLS certificate and private key
- enough host capacity for the configured session count and resource limits

The default capacity is four sessions at four CPUs and 24 GB RAM each. These are
starting limits, not measured guarantees for every dataset. Reduce them for a
smaller host and complete staged load testing before opening access.

The 61,844-cell `so-large.rds` stress fixture used approximately 18.3 GiB at peak
during import on Docker Desktop. Keep meaningful headroom above that measurement;
different objects with denser matrices or more assays can require substantially more.

On Docker Desktop for macOS, the mounted socket is normally owned by group `0`
inside Linux containers, so use `DOCKER_GID=0`. On Linux, use the socket group ID
reported by the command above; do not assume the values are interchangeable.

## Configure and start

1. Copy `deployment/.env.production.example` to
   `deployment/.env.production`, set `DOCKER_GID`, and provide a strong unique
   `CLUSTBUSTER_ADMIN_PASSWORD` and `CLUSTBUSTER_USER_PASSWORD`. The environment
   file is ignored by Git.
2. Place the certificate chain at `deployment/certs/fullchain.pem` and its private
   key at `deployment/certs/privkey.pem`. For local staging only, generate a
   short-lived self-signed certificate with:

   ```bash
   deployment/generate-local-certificate.sh deployment/certs clustbuster.localhost
   ```

3. Build the session image and validate the resolved configuration:

   ```bash
   docker compose --env-file deployment/.env.production \
     -f deployment/compose.production.yml --profile build build clustbuster-image
   docker compose --env-file deployment/.env.production \
     -f deployment/compose.production.yml config --quiet
   ```

4. Start the control plane:

   ```bash
   docker compose --env-file deployment/.env.production \
     -f deployment/compose.production.yml up -d shinyproxy traefik
   ```

HTTP redirects to HTTPS. Configure local DNS or a reverse-DNS entry for the host;
the dedicated endpoint accepts any hostname because it has only one router.
When using the example nonstandard staging ports, open
`https://127.0.0.1:8443` directly because an HTTP scheme redirect cannot infer the
separately mapped HTTPS port.

## Authentication

The included administrator and scientist simple-auth accounts are suitable for
isolated home-network staging and concurrent-session validation.
ShinyProxy stores simple-auth credentials in its runtime configuration, so use OIDC,
LDAP, or SAML before exposing ClustBuster beyond a trusted network. Keep identity
provider secrets outside the repository and restrict access to the environment file.

## Session lifecycle and privacy

- ShinyProxy creates one ClustBuster container per authenticated user and permits
  one instance per user.
- A missing browser heartbeat stops an app after 15 minutes by default.
- A session is stopped after eight hours and immediately on logout.
- The global default is four simultaneous sessions.
- Uploads, caches, state, and exports live only inside the session container. Removing
  the container removes this data. Users must download exports before logout or expiry.
- Session traffic to ShinyProxy stays on the internal `clustbuster-apps` network. A
  separate `clustbuster-egress` network permits outbound services such as Enrichr;
  Traefik is not attached to either app network.
- ShinyProxy receives the Docker Unix socket because its lifecycle role requires
  container create/start/stop/remove operations. Traefik uses only a static file
  provider and receives no Docker socket.

## Sovereign Atlas resources

For production, stage one immutable Sovereign Atlas release directory on the host,
mount it read-only at `/atlas` in every ClustBuster session container, and set:

```text
CLUSTBUSTER_MARKER_DB_PATH=/atlas/markercodex.duckdb
CLUSTBUSTER_REFERENCE_CATALOG_PATH=/atlas/reference_matrices.duckdb
CLUSTBUSTER_REFERENCE_FILES_ROOT=/atlas/reference_matrices/files
CLUSTBUSTER_ATLAS_VERSION=<release tag or commit>
```

Do not point running containers at a checkout being modified by the curator. Publish
or stage a complete new directory and restart sessions onto it; existing sessions
can finish against the previous immutable directory. This avoids mixing database
metadata with matrix files from different revisions and makes rollback a path or
image change.

## Upgrade and rollback

1. Back up `deployment/.env.production`, certificates, and the deployment YAML files.
   No uploaded datasets belong in backups.
2. Build a new immutable ClustBuster image tag and update `CLUSTBUSTER_IMAGE`.
3. Run the tests and `docker compose ... config --quiet`.
4. Stop active sessions during a maintenance window, then recreate ShinyProxy and
   Traefik. Roll back by restoring the previous image tag and configuration.

## Health and troubleshooting

```bash
docker compose --env-file deployment/.env.production \
  -f deployment/compose.production.yml ps
docker compose --env-file deployment/.env.production \
  -f deployment/compose.production.yml logs --tail=200 shinyproxy traefik
docker ps --filter label=app.kubernetes.io/name=clustbuster
```

- `permission denied` on `/var/run/docker.sock`: correct `DOCKER_GID` and recreate
  ShinyProxy; do not make the socket world-writable.
- `could not load FFI provider`: keep ShinyProxy's size-limited `/tmp` tmpfs
  executable. Its Docker Unix-socket client extracts a native jffi stub there.
- session container cannot be reached: verify both ShinyProxy and the session container
  joined `clustbuster-apps`.
- websocket disconnects: verify TLS terminates at Traefik and no upstream proxy imposes
  a shorter idle timeout.
- `OOMKilled`: inspect the dataset size and raise the session memory limit only when the
  host has corresponding capacity.
- certificate errors: verify the full chain/key paths and file permissions.

## Cleanup

Normal logout, heartbeat expiry, and maximum lifetime remove session containers.
After a maintenance stop, list any remaining ShinyProxy-managed containers before
removing them. Never prune volumes or images blindly on a host containing unrelated
workloads.

## Portable release images

`deployment/compose.portable.yml` runs a prebuilt image on localhost:8000 with a
read-only root filesystem and a writable named volume at `/data/gene_sets`.
ARM64 and AMD64 release images are packaged separately. Copy the matching image
archive and portable Compose file to the destination and use `docker load -i`.
Set `CLUSTBUSTER_IMAGE` and `CLUSTBUSTER_PLATFORM` for ARM64; defaults target AMD64.
No host Python or R installation is required.

The container entrypoint copies bundled enrichment libraries into the persistent
cache on first use. Existing files are retained across restarts and image upgrades;
Settings downloads can replace them explicitly. Both development Compose and
ShinyProxy session containers mount a named enrichment-cache volume. Do not remove
this volume unless you intend to discard downloaded libraries. Annotation state,
uploaded data, and exports remain session-scoped; download exports before stopping.

Before building a release, download the six human and four supported mouse libraries
in Settings, then stage Atlas with:

```sh
uv run python scripts/stage_release_resources.py /path/to/sovereign-atlas
```

This copies a snapshot and checksum manifest into ignored `resources/atlas/`.
The container activates this snapshot automatically; explicit `CLUSTBUSTER_*` Atlas
settings still override it. A clean checkout needs resources staged before packaging.
To replace an existing snapshot, move `resources/atlas/` aside before staging again.

Build and verify both platforms, for example:

```sh
docker buildx build --platform linux/arm64 --load -t clustbuster:2026.09.06-arm64 .
docker buildx build --platform linux/amd64 --load -t clustbuster:2026.09.06-amd64 .
scripts/container_smoke.sh clustbuster:2026.09.06-arm64
scripts/container_smoke.sh clustbuster:2026.09.06-amd64
```

The release smoke test requires the full staged resources. It checks non-root
startup on a read-only filesystem, HTTP health, native species offline enrichment,
MarkerCodex search, loading all bundled Refmats, and cache persistence in a second
container with networking disabled. Testing AMD64 under ARM emulation does not
replace validation on a physical destination host. Preserve the reference metadata
and source attribution when transferring the bundle.
