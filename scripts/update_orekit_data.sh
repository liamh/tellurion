#!/usr/bin/env bash
# Refresh the vendored Orekit data snapshot from upstream.
# Usage: scripts/update_orekit_data.sh [git-ref]   (default: main)
set -euo pipefail
REF="${1:-main}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$ROOT/tellurion/_orekit_data"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
git clone --depth 1 --branch "$REF" https://gitlab.orekit.org/orekit/orekit-data.git "$TMP/src"
COMMIT="$(git -C "$TMP/src" rev-parse HEAD)"
README="$(cat "$DEST/README.md")"
rm -rf "$DEST"
mkdir -p "$DEST"
cp -r "$TMP/src/." "$DEST/"
rm -rf "$DEST/.git"
printf '%s\n' "$README" > "$DEST/README.md"
echo "$COMMIT" > "$DEST/SNAPSHOT_COMMIT"
echo "Orekit data updated to $COMMIT"
