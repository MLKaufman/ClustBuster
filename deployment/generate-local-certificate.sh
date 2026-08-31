#!/bin/sh
set -eu

output_dir="${1:-deployment/certs}"
hostname="${2:-clustbuster.localhost}"
mkdir -p "$output_dir"
openssl req -x509 -newkey rsa:3072 -sha256 -nodes -days 30 \
  -keyout "$output_dir/privkey.pem" \
  -out "$output_dir/fullchain.pem" \
  -subj "/CN=$hostname" \
  -addext "subjectAltName=DNS:$hostname,DNS:localhost,IP:127.0.0.1"
chmod 600 "$output_dir/privkey.pem"
chmod 644 "$output_dir/fullchain.pem"
printf 'Wrote a 30-day local certificate for %s to %s\n' "$hostname" "$output_dir"
