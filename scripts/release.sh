#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

crate="quantion"

version="$(sed -n 's/^version = "\(.*\)"/\1/p' core/Cargo.toml | head -1)"
[ -n "$version" ] || { echo "no version in core/Cargo.toml"; exit 1; }
[ -f core/Cargo.lock ] || { echo "core/Cargo.lock is missing — it must be committed"; exit 1; }
git diff --quiet -- core/Cargo.lock || { echo "core/Cargo.lock changed — commit it before releasing"; exit 1; }

out="artifacts/$version"

# a run that dies partway leaves a directory holding some fresh binaries and
# some stale ones. that is not a release, and the manifest it still carried
# would describe binaries that are no longer there, so the build phase drops
# the manifest and only the final step writes it back.
sums=""
cleanup() {
  code=$?
  [ -n "$sums" ] && rm -f "$sums"
  if [ "$code" -ne 0 ] && [ ! -f "$out/manifest.json" ]; then
    echo
    echo "release did not finish — $out has no manifest.json, so no wrapper will"
    echo "load it and 'make publish' will refuse it."
    echo "re-run 'make release' to finish: every target that already built is"
    echo "recorded and will be skipped."
  fi
}
trap cleanup EXIT

make --no-print-directory doctor
make --no-print-directory check-version
make --no-print-directory test
make --no-print-directory header
make --no-print-directory header-check

builds="macos-arm64 macos-arm64/lib$crate.dylib
macos-x86_64 macos-x86_64/lib$crate.dylib
linux-amd64 linux-x86_64/lib$crate.so
linux-arm64 linux-arm64/lib$crate.so
windows-amd64 windows-x86_64/lib$crate.dll
wasm wasm/$crate.wasm"

# everything the cdylib is built from. a change to any of it makes every
# binary in artifacts/<version> stale, whatever its bytes still look like.
sources="core/Cargo.toml
core/Cargo.lock
core/build.rs
core/cbindgen.toml
core/.cargo/config.toml
core/src"

fingerprint_sources() {
  SOURCES="$sources" python3 - <<'PY'
import hashlib, os


def walk(root):
    if os.path.isfile(root):
        yield root
        return
    for base, dirs, names in os.walk(root):
        dirs.sort()
        for name in sorted(names):
            yield os.path.join(base, name)


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


paths = sorted({p for root in os.environ["SOURCES"].splitlines() for p in walk(root)})
listing = "\n".join(f"{file_digest(path)}  {path}" for path in paths)
print("sha256:" + hashlib.sha256(listing.encode()).hexdigest())
PY
}

# which sources each binary was built from, one entry per target. it lives
# outside artifacts/<version> so it never ships inside a release, and it
# outlives a failed run so the retry rebuilds only what did not finish.
state="artifacts/.build-state.json"

record_built() {
  STATE="$state" VERSION="$version" FINGERPRINT="$fingerprint" TARGET="$1" python3 - <<'PY'
import json, os

path = os.environ["STATE"]
try:
    with open(path) as handle:
        state = json.load(handle)
except (OSError, ValueError):
    state = {}
if not isinstance(state, dict):
    state = {}

state.setdefault(os.environ["VERSION"], {})[os.environ["TARGET"]] = os.environ["FINGERPRINT"]

temporary = path + ".tmp"
with open(temporary, "w") as handle:
    json.dump(state, handle, indent=2, sort_keys=True)
os.replace(temporary, path)
PY
}

# targets that still need building: the ones built from other sources than the
# ones on disk now, plus any whose binary went missing, was truncated, is not a
# library for its platform, or drifted from the manifest that describes it.
stale_targets() {
  OUT="$out" BUILDS="$builds" STATE="$state" VERSION="$version" FINGERPRINT="$fingerprint" \
  python3 - <<'PY'
import hashlib, json, os

out = os.environ["OUT"]
fingerprint = os.environ["FINGERPRINT"]


def load(path, default):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return default


built_from = load(os.environ["STATE"], {}).get(os.environ["VERSION"], {})
recorded = load(os.path.join(out, "manifest.json"), {}).get("files", {})

first_bytes = {
    ".dylib": [b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xca\xfe\xba\xbe"],
    ".so": [b"\x7fELF"],
    ".dll": [b"MZ"],
    ".wasm": [b"\x00asm"],
}


def looks_built(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    allowed = first_bytes.get(os.path.splitext(path)[1])
    if allowed is None:
        return True
    with open(path, "rb") as handle:
        head = handle.read(4)
    return any(head.startswith(start) for start in allowed)


def matches_manifest(path, name):
    wanted = recorded.get(name)
    if wanted is None:
        return True
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return wanted == "sha256:" + digest.hexdigest()


for line in os.environ["BUILDS"].splitlines():
    target, name = line.split(" ", 1)
    path = os.path.join(out, name)
    if built_from.get(target) != fingerprint:
        print(target)
    elif not looks_built(path) or not matches_manifest(path, name):
        print(target)
PY
}

fingerprint="$(fingerprint_sources)"
stale="$(stale_targets)"

if [ -z "$stale" ]; then
  echo "every binary in $out was built from the current core sources"
else
  echo "building from the current core sources: $(echo $stale)"
  # from here until the manifest is written the directory is a mix of fresh and
  # stale binaries. every wrapper finds a version by its manifest.json, so
  # dropping it hides the half-built directory instead of letting anything load
  # it, and publish.sh refuses to ship it.
  rm -f "$out/manifest.json"
  for target in $stale; do
    echo "building $target"
    OUT="$out" make --no-print-directory "$target"
    record_built "$target"
  done
fi

left="$(stale_targets)"
[ -z "$left" ] || { echo "still missing after building: $(echo "$left" | tr '\n' ' ')"; exit 1; }

QUANTION_ARTIFACTS_ROOT="$root/$out" make --no-print-directory check-api

if [ -f artifacts/include/quantion.h ]; then
  mkdir -p "$out/include"
  cp artifacts/include/quantion.h "$out/include/"
fi

sums="$(mktemp)"
( cd "$out" && find . -type f ! -name manifest.json | sed 's#^\./##' | sort | while read -r f; do
    echo "$(shasum -a 256 "$f" | cut -d' ' -f1)  $f"
  done ) > "$sums"

git_revision="$(git rev-parse HEAD 2>/dev/null || echo unknown)"
ionic_revision="$(sed -n 's/.*ionic[^#]*#\([0-9a-f]\{7,40\}\).*/\1/p' core/Cargo.lock 2>/dev/null | head -1)"
[ -n "$ionic_revision" ] || { echo "could not read the ionic revision from core/Cargo.lock"; exit 1; }

VERSION="$version" GIT="$git_revision" IONIC="$ionic_revision" SOURCE="$fingerprint" \
OUT="$out/manifest.json" SUMS="$sums" python3 - <<'PY'
import json, os
config = json.load(open("config.json"))[0]
files = {
    line.split("  ", 1)[1].rstrip(): "sha256:" + line.split("  ", 1)[0]
    for line in open(os.environ["SUMS"])
}
manifest = {
    "version": os.environ["VERSION"],
    "format": config["format"],
    "description": config["description"],
    "git": os.environ["GIT"],
    "ionic": os.environ["IONIC"],
    "source": os.environ["SOURCE"],
    "files": files,
}

# writing the manifest is what turns the directory back into a release, so it
# lands whole or not at all — a half-written one would be worse than none.
path = os.environ["OUT"]
temporary = path + ".tmp"
with open(temporary, "w") as handle:
    json.dump(manifest, handle, indent=2)
os.replace(temporary, path)
PY

echo
echo "built $out (manifest.json lists a sha256 for every file)"
echo "next: review 'git diff', commit, then run:  make publish VERSION=$version"
