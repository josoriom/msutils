#!/usr/bin/env bash
# fail if a fresh build changes any copy of quantion.h
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

copies=(
  artifacts/include/quantion.h
  wrappers/js/src/include/quantion.h
  wrappers/r/src/include/quantion.h
)

saved="$(mktemp -d)"
trap 'rm -rf "$saved"' EXIT
for copy in "${copies[@]}"; do
  [ -f "$copy" ] || { echo "missing header: $copy"; exit 1; }
  cp "$copy" "$saved/$(echo "$copy" | tr / _)"
done

cargo build --locked --manifest-path core/Cargo.toml >/dev/null

failed=0
for copy in "${copies[@]}"; do
  if ! cmp -s "$copy" "$saved/$(echo "$copy" | tr / _)"; then
    echo "stale header: $copy differs from a fresh build"
    failed=1
  fi
done

reference="${copies[0]}"
for copy in "${copies[@]:1}"; do
  cmp -s "$reference" "$copy" || { echo "header copies differ: $reference vs $copy"; failed=1; }
done

[ "$failed" -eq 0 ] || exit 1
echo "every quantion.h matches a fresh build"
