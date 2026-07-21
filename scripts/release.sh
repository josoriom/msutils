#!/usr/bin/env bash
# check, build every platform, write manifest.json into artifacts/<version>/
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

version="$(sed -n 's/^version = "\(.*\)"/\1/p' core/Cargo.toml | head -1)"
[ -n "$version" ] || { echo "no version in core/Cargo.toml"; exit 1; }
[ -f core/Cargo.lock ] || { echo "core/Cargo.lock is missing — it must be committed"; exit 1; }
git diff --quiet -- core/Cargo.lock || { echo "core/Cargo.lock changed — commit it before releasing"; exit 1; }
out="artifacts/$version"

make --no-print-directory doctor
make --no-print-directory check-version
make --no-print-directory test
make --no-print-directory header-check

staging="artifacts/.building-$version"
rm -rf "$staging" "$out"
trap 'rm -rf "$staging"' EXIT
OUT="$staging" make --no-print-directory all
QUANTION_ARTIFACTS_ROOT="$root/$staging" make --no-print-directory check-api
[ -d "$staging" ] || { echo "nothing was built into $staging"; exit 1; }

if [ -f artifacts/include/quantion.h ]; then
  mkdir -p "$staging/include"
  cp artifacts/include/quantion.h "$staging/include/"
fi

sums="$(mktemp)"
trap 'rm -rf "$staging"; rm -f "$sums"' EXIT
( cd "$staging" && find . -type f ! -name manifest.json | sed 's#^\./##' | sort | while read -r f; do
    echo "$(shasum -a 256 "$f" | cut -d' ' -f1)  $f"
  done ) > "$sums"

git_revision="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
ionic_revision="$(sed -n 's/.*ionic[^#]*#\([0-9a-f]\{7,40\}\).*/\1/p' core/Cargo.lock 2>/dev/null | head -1)"
[ -n "$ionic_revision" ] || { echo "could not read the ionic revision from core/Cargo.lock"; exit 1; }

VERSION="$version" GIT="$git_revision" IONIC="$ionic_revision" \
OUT="$staging/manifest.json" SUMS="$sums" python3 - <<'PY'
import json, os
config = json.load(open("config.json"))[0]
files = {
    line.split("  ", 1)[1].rstrip(): "sha256:" + line.split("  ", 1)[0]
    for line in open(os.environ["SUMS"])
}
json.dump({
    "version": os.environ["VERSION"],
    "format": config["format"],
    "description": config["description"],
    "git": os.environ["GIT"],
    "ionic": os.environ["IONIC"],
    "files": files,
}, open(os.environ["OUT"], "w"), indent=2)
PY

mv "$staging" "$out"
trap 'rm -f "$sums"' EXIT

echo
echo "built $out (manifest.json lists a sha256 for every file)"
echo "next: review 'git diff', commit, then run:  make publish VERSION=$version"
