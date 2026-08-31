#!/usr/bin/env bash
# fail if the 5 version numbers disagree, or format.current != the ABI version
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

wanted="$(sed -n 's/^version = "\(.*\)"/\1/p' core/Cargo.toml | head -1)"
[ -n "$wanted" ] || { echo "no version in core/Cargo.toml"; exit 1; }

[ -f config.json ] || {
  echo "config.json missing - create it with: make publish VERSION=$wanted"
  exit 1
}

js="$(python3 -c "import json;print(json.load(open('wrappers/js/package.json'))['version'])")"
py="$(sed -n 's/^version = "\(.*\)"/\1/p' wrappers/python/pyproject.toml | head -1)"
r="$(sed -n 's/^Version: *\(.*\)/\1/p' wrappers/r/DESCRIPTION | head -1)"
pyinit="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' wrappers/python/quantion/__init__.py | head -1)"
config="$(python3 -c "import json;print(json.load(open('config.json'))[0]['package'])")"
abi="$(sed -n 's/^pub const QUANTION_ABI_VERSION[^=]*=[^0-9]*\([0-9][0-9]*\).*/\1/p' core/src/ffi.rs | head -1)"

failed=0
check() {
  [ "$2" = "$wanted" ] || { echo "version mismatch: $1 = '$2' (expected '$wanted' from Cargo.toml)"; failed=1; }
}
check "wrappers/js/package.json" "$js"
check "wrappers/python/pyproject.toml" "$py"
check "wrappers/python/quantion/__init__.py" "$pyinit"
check "wrappers/r/DESCRIPTION" "$r"
check "config.json[0].package" "$config"

[ "$failed" -eq 0 ] || exit 1
echo "all versions agree at $wanted (abi $abi)"
