#!/usr/bin/env bash
set -euo pipefail

image="${1:-clustbuster:test}"
container_name="clustbuster-smoke-${RANDOM}-$$"

cleanup() {
  docker rm -f "${container_name}" >/dev/null 2>&1 || true
}
trap cleanup EXIT

image_uid="$(docker run --rm --entrypoint id "${image}" -u)"
if [[ "${image_uid}" == "0" ]]; then
  echo "Container image runs as root" >&2
  exit 1
fi

docker run -d \
  --name "${container_name}" \
  --read-only \
  --tmpfs /workspace:size=256M,mode=1777 \
  --tmpfs /tmp:size=128M,mode=1777 \
  --security-opt no-new-privileges:true \
  --cap-drop ALL \
  -p 127.0.0.1::8000 \
  "${image}" >/dev/null

for _ in {1..30}; do
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

echo "ClustBuster container smoke test passed"
