FROM python:3.13-slim AS builder

COPY --from=ghcr.io/astral-sh/uv:0.12.7 /uv /usr/local/bin/uv
RUN apt-get update \
    && apt-get install --yes --no-install-recommends git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_NO_CACHE=1
WORKDIR /build
COPY pyproject.toml uv.lock README.md LICENSE ./
COPY src ./src
RUN uv sync --frozen --no-dev --no-editable

FROM python:3.13-slim
ENV PATH=/opt/venv/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CLUSTBUSTER_HOST=0.0.0.0 \
    CLUSTBUSTER_PORT=8000 \
    CLUSTBUSTER_WORKSPACE_ROOT=/workspace \
    MPLCONFIGDIR=/tmp/matplotlib
RUN groupadd --system clustbuster \
    && useradd --system --gid clustbuster --home-dir /app clustbuster \
    && mkdir -p /app /workspace \
    && chown -R clustbuster:clustbuster /app /workspace
COPY --from=builder /opt/venv /opt/venv
COPY --chown=clustbuster:clustbuster resources /app/resources
WORKDIR /app
USER clustbuster
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/', timeout=3)" || exit 1
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "8000", "clustbuster.app:app"]
