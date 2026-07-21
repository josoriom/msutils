#!/usr/bin/env bash
# fail if a wrapper reads a remote ion file wrongly, or reads more than it needs
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
here="$root/scripts/check-remote"
cd "$root"

export QUANTION_ARTIFACTS_ROOT="${QUANTION_ARTIFACTS_ROOT:-$root/artifacts}"

fixture="core/tests/fixtures/api/api.ion"
[ -f "$fixture" ] || { echo "fixture missing: $fixture"; exit 1; }

echo "check-remote — one file over HTTP, only the bytes each query needs"
echo

work_dir="$(mktemp -d)"
r_work=""
languages=()
missing=()
failed=0

run() {
  local name="$1" script="$2"
  shift 2
  if [ ! -f "$here/$script" ]; then
    missing+=("$name")
    return
  fi
  echo "$name:"
  local measured="$work_dir/$name.txt"
  if "$@" "$here/$script" > "$measured"; then
    if python3 "$here/judge.py" "$measured" "$name"; then
      languages+=("$name")
    else
      failed=1
    fi
  else
    echo "  the $name runner failed"
    failed=1
  fi
  echo
}

server_log="$(mktemp)"
python3 - "$fixture" > "$server_log" <<'PY' &
import sys, signal
sys.path.insert(0, "scripts/check-remote")
from server import RangeServer
server = RangeServer(sys.argv[1])
print(server.url, flush=True)
signal.pause()
PY
server_pid=$!

clean_up() {
  kill "$server_pid" 2>/dev/null || true
  wait "$server_pid" 2>/dev/null || true
  rm -rf "$work_dir" "$r_work" "$server_log"
}
trap clean_up EXIT

for _ in $(seq 1 50); do
  [ -s "$server_log" ] && break
  sleep 0.1
done
export QUANTION_REMOTE_URL="$(head -1 "$server_log")"
[ -n "$QUANTION_REMOTE_URL" ] || { echo "the range server did not start"; exit 1; }

run python run_python.py python3

if command -v Rscript >/dev/null 2>&1 && [ -f "$here/run_r.R" ]; then
  r_work="$(mktemp -d)"
  cp -R "$root/wrappers/r" "$r_work/rpkg"
  if R CMD INSTALL --no-docs --no-byte-compile -l "$r_work" "$r_work/rpkg" > "$r_work/install.log" 2>&1; then
    r_user_lib="$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"
    r_path_sep="$(Rscript -e 'cat(.Platform$path.sep)')"
    export R_LIBS_USER="$r_work${r_user_lib:+$r_path_sep$r_user_lib}"
    run r run_r.R Rscript
  else
    echo "the R package failed to install:"
    tail -20 "$r_work/install.log"
    failed=1
  fi
else
  missing+=("r")
fi

run js run_js.mjs node

python3 - "$work_dir" ${languages[@]+"${languages[@]}"} <<'PY' || failed=1
import sys, os
work, names = sys.argv[1], sys.argv[2:]
keys = ("opening_bytes", "opening_requests", "query_bytes", "query_requests")
measured = {}
for name in names:
    values = {}
    for line in open(os.path.join(work, name + ".txt")):
        key, sep, value = line.strip().partition(" = ")
        if sep:
            values[key.strip()] = value.strip()
    measured[name] = values
problems = []
for key in keys:
    seen = {name: measured[name].get(key) for name in names}
    if len(set(seen.values())) > 1:
        problems.append(f"  FAIL  {key} differs across languages: {seen}")
print("\n".join(problems) if problems else "  ok  every language fetched the same bytes in the same requests")
sys.exit(1 if problems else 0)
PY

if [ ${#missing[@]} -gt 0 ]; then
  echo "NOT CHECKED: ${missing[*]} — no runner in scripts/check-remote/"
  echo "   those languages read remote files unchecked"
fi

if [ "$failed" -ne 0 ]; then
  echo "check-remote failed"
  exit 1
fi

echo "check-remote passed for: ${languages[*]-}"
