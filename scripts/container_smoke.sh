#!/usr/bin/env bash
set -euo pipefail

image="${1:-clustbuster:test}"
container_name="clustbuster-smoke-${RANDOM}-$$"

cleanup() {
  docker rm -fv "${container_name}" >/dev/null 2>&1 || true
  docker volume rm "${container_name}-cache" >/dev/null 2>&1 || true
}
trap cleanup EXIT

platform="linux/$(docker image inspect --format '{{.Architecture}}' "${image}")"
image_uid="$(docker run --rm --platform "${platform}" --entrypoint id "${image}" -u)"
if [[ "${image_uid}" == "0" ]]; then
  echo "Container image runs as root" >&2
  exit 1
fi

docker run -d --platform "${platform}" \
  --name "${container_name}" \
  --read-only \
  -v "${container_name}-cache:/data/gene_sets" \
  --tmpfs /workspace:size=256M,mode=1777 \
  --tmpfs /tmp:size=128M,mode=1777 \
  --security-opt no-new-privileges:true \
  --cap-drop ALL \
  -p 127.0.0.1::8000 \
  "${image}" >/dev/null

for _ in {1..60}; do
  health="$(docker inspect --format '{{.State.Health.Status}}' "${container_name}")"
  if [[ "${health}" == "healthy" ]]; then
    break
  fi
  sleep 1
done

if [[ "${health:-}" != "healthy" ]]; then
  docker logs "${container_name}" >&2
  echo "Container did not become healthy" >&2
  exit 1
fi

runtime_uid="$(docker exec "${container_name}" id -u)"
if [[ "${runtime_uid}" == "0" ]]; then
  echo "Running container process is root" >&2
  exit 1
fi

host_port="$(docker port "${container_name}" 8000/tcp | sed 's/.*://')"
curl --fail --silent --show-error "http://127.0.0.1:${host_port}/" | grep --quiet ClustBuster
docker exec "${container_name}" test -w /workspace

docker exec -i "${container_name}" python /usr/local/bin/clustbuster-entrypoint.py python < "$(dirname "$0")/container_resource_smoke.py"
docker run --rm --platform "${platform}" --network none \
  -v "${container_name}-cache:/data/gene_sets" "${image}" \
  python -c 'from pathlib import Path; assert Path("/data/gene_sets/persistence-smoke.txt").read_text() == "persisted"'
echo "ClustBuster container and persistent resource smoke tests passed (${platform})"
