#!/usr/bin/env bash
# build for the machine you are on, into artifacts/<version>/<platform>/
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

version="$(sed -n 's/^version = "\(.*\)"/\1/p' core/Cargo.toml | head -1)"
[ -n "$version" ] || { echo "no version in core/Cargo.toml"; exit 1; }

system="$(uname -s)"
machine="$(uname -m)"

case "$system" in
  Darwin)
    case "$machine" in
      arm64|aarch64) platform="macos-arm64" ;;
      *)             platform="macos-x86_64" ;;
    esac
    library="libquantion.dylib"
    ;;
  Linux)
    case "$machine" in
      aarch64|arm64) platform="linux-arm64" ;;
      *)             platform="linux-x86_64" ;;
    esac
    library="libquantion.so"
    ;;
  MINGW*|MSYS*|CYGWIN*|Windows_NT)
    case "$machine" in
      aarch64|arm64) platform="windows-arm64" ;;
      *)             platform="windows-x86_64" ;;
    esac
    library="quantion.dll"
    ;;
  *)
    echo "unsupported system: $system $machine"
    exit 1
    ;;
esac

cargo build --locked --manifest-path core/Cargo.toml --release

out="artifacts/$version/$platform"
mkdir -p "$out"

built="core/target/release/$library"
[ -f "$built" ] || { echo "cargo did not produce $built"; exit 1; }

# every wrapper looks for libquantion.dll on Windows, but cargo emits quantion.dll
case "$platform" in
  windows-*) cp "$built" "$out/libquantion.dll" ;;
  *)         cp "$built" "$out/$library" ;;
esac

rm -f "artifacts/$version/manifest.json"

echo "built $out"
