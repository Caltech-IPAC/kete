#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

if [[ ! -d "$HOME/.kete" ]]; then
  echo "Kete cache directory not found. Please run the pytests first." >&2
  echo "You can run the tests with:" >&2
  echo "    uv run pytest" >&2
  exit 1
fi

if ! command -v podman &>/dev/null; then
  echo "Podman is not installed. Please install podman to build the container." >&2
  exit 1
fi

if [[ "${1:-}" == "--push" || "${1:-}" == "-p" ]]; then
  PUSH=true
else
  PUSH=false
fi

CONTAINER_REPO="ghcr.io/caltech-ipac"
CONTAINER_NAME="kete-cache"

# Container image metadata
TITLE="Kete cache"
DESCRIPTION="Kete cache as an OCI container image"
VENDOR="Caltech/IPAC"
REPO="https://github.com/Caltech-IPAC/kete"
URL="${REPO}/pkgs/container/kete-cache"
DOCS="https://caltech-ipac.github.io/kete"
CALVER="$(date '+%Y.%m.%d')"
CREATED_AT="$(date --rfc-3339=date)"
AUTHORS="$("${SCRIPT_DIR}"/get-authors.py)"
VERSION="$("${SCRIPT_DIR}"/get-version.py)"

declare -a images=(
  "${CONTAINER_REPO}/${CONTAINER_NAME}:latest"
  "${CONTAINER_REPO}/${CONTAINER_NAME}:${CALVER}"
  "${CONTAINER_REPO}/${CONTAINER_NAME}:${VERSION}"
)

declare -a labels=(
  "org.opencontainers.image.created=${CREATED_AT}"
  "org.opencontainers.image.authors=${AUTHORS}"
  "org.opencontainers.image.url=${URL}"
  "org.opencontainers.image.documentation=${DOCS}"
  "org.opencontainers.image.source=${REPO}"
  "org.opencontainers.image.version=${VERSION}"
  "org.opencontainers.image.vendor=${VENDOR}"
  "org.opencontainers.image.licenses=BSD-3-Clause"
  "org.opencontainers.image.title=${TITLE}"
  "org.opencontainers.image.description=${DESCRIPTION}"
)

declare -a build_args
declare -a manifest_args

for label in "${labels[@]}"; do
  build_args+=(
    --label="$label"
    --annotation="$label"
  )
  manifest_args+=(
    --annotation="$label"
  )
done

for image in "${images[@]}"; do
  if podman manifest exists "$image" &>/dev/null; then
    podman manifest rm "$image"
  fi
  podman manifest create "${manifest_args[@]}" "$image"
done

for arch in "amd64" "arm64" "ppc64le" "s390x" "riscv64"; do
  podman build \
    --platform="linux/$arch" \
    -t "${CONTAINER_REPO}/${CONTAINER_NAME}:${CALVER}-${arch}" \
    "${build_args[@]}" \
    -f "${SCRIPT_DIR}/Containerfile" "$HOME/.kete"
  for image in "${images[@]}"; do
    podman manifest add "${manifest_args[@]}" "$image" "${CONTAINER_REPO}/${CONTAINER_NAME}:${CALVER}-${arch}"
  done
done

if [[ "$PUSH" == true ]]; then
  for image in "${images[@]}"; do
    podman manifest push --format=oci --all "$image" "docker://${image}"
  done
fi
