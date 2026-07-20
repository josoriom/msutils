#!/usr/bin/env bash
# set the version everywhere, release, commit, tag, push. takes VERSION=X.Y.Z
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

version="${1:-}"
[ -n "$version" ] || { echo "usage: make publish VERSION=X.Y.Z   (e.g. make publish VERSION=0.2.0)"; exit 1; }
version="${version#v}"
major="${version%%.*}"
rest="${version#*.}"
minor="${rest%%.*}"
patch="${rest#*.}"

for part in "$major" "$minor" "$patch"; do
  case "$part" in
    ''|*[!0-9]*) echo "version must be three numbers, got '$version'"; exit 1 ;;
  esac
done

branch="$(git rev-parse --abbrev-ref HEAD)"
[ "$branch" = "main" ] || { echo "publish only from main, you are on '$branch'"; exit 1; }

[ -z "$(git status --porcelain)" ] || {
  echo "the working tree is not clean — commit or stash first:"
  git status --short
  exit 1
}

git rev-parse "v$version" >/dev/null 2>&1 && { echo "tag v$version already exists"; exit 1; }
git fetch --quiet origin
[ "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)" ] || {
  echo "main is not in sync with origin/main — pull or push first"; exit 1
}

message="feat: release v$version"

echo "publishing v$version"
echo

echo "setting version $version everywhere..."
sed -i.bak "s/^version = \".*\"/version = \"$version\"/" core/Cargo.toml && rm -f core/Cargo.toml.bak
sed -i.bak "s/^version = \".*\"/version = \"$version\"/" wrappers/python/pyproject.toml && rm -f wrappers/python/pyproject.toml.bak
sed -i.bak "s/^__version__ = \".*\"/__version__ = \"$version\"/" wrappers/python/quantion/__init__.py \
  && rm -f wrappers/python/quantion/__init__.py.bak
sed -i.bak "s/^Version: .*/Version: $version/" wrappers/r/DESCRIPTION && rm -f wrappers/r/DESCRIPTION.bak

VERSION="$version" python3 - <<'PY'
import json, os
path = "wrappers/js/package.json"
package = json.load(open(path))
package["version"] = os.environ["VERSION"]
json.dump(package, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY

VERSION="$version" DESCRIPTION="$message" python3 - <<'PY'
import json, os
path = "config.json"
entries = json.load(open(path)) if os.path.exists(path) else []
fmt = entries[0]["format"] if entries else {"current": 1, "min_supported": 1, "max_supported": 1}
entries = [e for e in entries if e["package"] != os.environ["VERSION"]]
entries.insert(0, {
    "package": os.environ["VERSION"],
    "format": fmt,
    "description": os.environ["DESCRIPTION"],
})
json.dump(entries, open(path, "w"), indent=2)
open(path, "a").write("\n")
PY

echo "version set in Cargo.toml, pyproject.toml, DESCRIPTION, package.json, config.json"
echo

make --no-print-directory release

echo
echo "verifying artifact checksums against manifest.json..."
OUT="artifacts/$version" python3 - <<'PY'
import json, hashlib, os, sys
directory = os.environ["OUT"]
manifest = json.load(open(os.path.join(directory, "manifest.json")))
bad = [
    name for name, want in manifest["files"].items()
    if "sha256:" + hashlib.sha256(open(os.path.join(directory, name), "rb").read()).hexdigest() != want
]
sys.exit("checksum mismatch: " + ", ".join(bad)) if bad else print("checksums match")
PY

echo
git add core/Cargo.toml core/Cargo.lock wrappers/python/pyproject.toml \
        wrappers/python/quantion/__init__.py wrappers/r/DESCRIPTION \
        wrappers/js/package.json config.json "artifacts/$version"
git commit -m "$message"
git tag -a "v$version" -m "$message"
git push origin HEAD
git push origin "v$version"

echo
echo "published v$version"
