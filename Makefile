CRATE          := quantion
CRATE_MANIFEST := core/Cargo.toml
ARTIFACTS      := artifacts
LLVM_PREFIX    ?= /opt/homebrew/opt/llvm

DOCKER_BULLSEYE_IMAGE := rust:1-bullseye
CARGO_CACHE           := $(HOME)/.cache/quantion-cargo
WINDOWS_ARM64_IMAGE   := dockcross/windows-arm64:latest
WINDOWS_ARM64_TARGET  := aarch64-pc-windows-gnullvm

CRATE_VERSION := $(shell sed -n 's/^version = "\(.*\)"/\1/p' $(CRATE_MANIFEST) | head -1)
OUT           ?= $(ARTIFACTS)/$(CRATE_VERSION)

.DEFAULT_GOAL := help

.PHONY: help all mac linux windows wasm \
        macos-arm64 macos-x86_64 linux-amd64 linux-arm64 windows-amd64 windows-arm64 \
        doctor check-version test header-check native quick-check check-api check-api-update \
        release publish clean

# build every platform into artifacts/<version>/
all: mac linux windows wasm

# build both macOS targets
mac:     macos-arm64 macos-x86_64
# build both Linux targets
linux:   linux-amd64 linux-arm64
# build Windows x86_64
windows: windows-amd64

# Apple Silicon -> libquantion.dylib
macos-arm64:
	rustup target add aarch64-apple-darwin
	cargo build --locked --manifest-path $(CRATE_MANIFEST) --release --target aarch64-apple-darwin
	mkdir -p $(OUT)/macos-arm64
	cp core/target/aarch64-apple-darwin/release/lib$(CRATE).dylib $(OUT)/macos-arm64/

# Intel Mac -> libquantion.dylib
macos-x86_64:
	rustup target add x86_64-apple-darwin
	cargo build --locked --manifest-path $(CRATE_MANIFEST) --release --target x86_64-apple-darwin
	mkdir -p $(OUT)/macos-x86_64
	cp core/target/x86_64-apple-darwin/release/lib$(CRATE).dylib $(OUT)/macos-x86_64/

# Linux x86_64 in Docker -> libquantion.so
linux-amd64:
	mkdir -p $(CARGO_CACHE)
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/core/target-linux-amd64 \
	  -e CARGO_HOME=/cargo \
	  -u $$(id -u):$$(id -g) \
	  -v $(CARGO_CACHE):/cargo \
	  -v $$PWD:/work -w /work \
	  --entrypoint /bin/bash $(DOCKER_BULLSEYE_IMAGE) -lc '\
	    set -e; \
	    export PATH=/usr/local/cargo/bin:$$PATH; \
	    cargo build --locked --manifest-path $(CRATE_MANIFEST) --release; \
	    strip /work/core/target-linux-amd64/release/lib$(CRATE).so'
	mkdir -p $(OUT)/linux-x86_64
	cp core/target-linux-amd64/release/lib$(CRATE).so $(OUT)/linux-x86_64/

# Linux arm64 in Docker -> libquantion.so
linux-arm64:
	mkdir -p $(CARGO_CACHE)
	docker run --rm --platform=linux/arm64 \
	  -e CARGO_TARGET_DIR=/work/core/target-linux-arm64 \
	  -e CARGO_HOME=/cargo \
	  -u $$(id -u):$$(id -g) \
	  -v $(CARGO_CACHE):/cargo \
	  -v $$PWD:/work -w /work \
	  --entrypoint /bin/bash $(DOCKER_BULLSEYE_IMAGE) -lc '\
	    set -e; \
	    export PATH=/usr/local/cargo/bin:$$PATH; \
	    cargo build --locked --manifest-path $(CRATE_MANIFEST) --release; \
	    strip /work/core/target-linux-arm64/release/lib$(CRATE).so'
	mkdir -p $(OUT)/linux-arm64
	cp core/target-linux-arm64/release/lib$(CRATE).so $(OUT)/linux-arm64/

# browser target -> quantion.wasm (needs LLVM)
wasm:
	rustup target add wasm32-unknown-unknown
	CC_wasm32_unknown_unknown="$(LLVM_PREFIX)/bin/clang" \
	AR_wasm32_unknown_unknown="$(LLVM_PREFIX)/bin/llvm-ar" \
	CFLAGS_wasm32_unknown_unknown="--target=wasm32-unknown-unknown" \
	cargo build --locked --manifest-path $(CRATE_MANIFEST) --release --target wasm32-unknown-unknown
	mkdir -p $(OUT)/wasm
	cp core/target/wasm32-unknown-unknown/release/$(CRATE).wasm $(OUT)/wasm/

# Windows x86_64 in Docker -> libquantion.dll
windows-amd64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/core/target-windows-amd64 \
	  -v $$PWD:/work -w /work $(DOCKER_BULLSEYE_IMAGE) bash -lc '\
	    set -euo pipefail; \
	    export PATH=/usr/local/cargo/bin:$$PATH; \
	    apt-get update && apt-get install -y --no-install-recommends \
	      curl ca-certificates gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 mingw-w64; \
	    if ! command -v rustup >/dev/null 2>&1; then \
	      curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; \
	    fi; \
	    rustup target add x86_64-pc-windows-gnu; \
	    RUSTFLAGS="-C linker=x86_64-w64-mingw32-gcc" \
	      cargo build --locked --manifest-path $(CRATE_MANIFEST) --release --target x86_64-pc-windows-gnu'
	mkdir -p $(OUT)/windows-x86_64
	cp core/target-windows-amd64/x86_64-pc-windows-gnu/release/$(CRATE).dll $(OUT)/windows-x86_64/lib$(CRATE).dll

# Windows arm64. opt-in: no wrapper can load it, so `all` skips it
windows-arm64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/core/target-windows-arm64 \
	  -v $$PWD:/work -w /work $(WINDOWS_ARM64_IMAGE) \
	  bash -lc 'set -euo pipefail; \
	    export PATH="$$HOME/.cargo/bin:/usr/local/cargo/bin:$$PATH"; \
	    command -v rustup >/dev/null 2>&1 || { curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; . "$$HOME/.cargo/env"; }; \
	    rustup target add $(WINDOWS_ARM64_TARGET); \
	    cargo build --locked --manifest-path $(CRATE_MANIFEST) --release --target $(WINDOWS_ARM64_TARGET)'
	mkdir -p $(OUT)/windows-arm64
	cp core/target-windows-arm64/$(WINDOWS_ARM64_TARGET)/release/$(CRATE).dll $(OUT)/windows-arm64/lib$(CRATE).dll

# fail if a build tool is missing
doctor:
	@missing=""; \
	for t in cargo rustup docker python3 node; do \
	  command -v $$t >/dev/null 2>&1 || missing="$$missing $$t"; \
	done; \
	command -v Rscript >/dev/null 2>&1 || missing="$$missing Rscript"; \
	[ -x "$(LLVM_PREFIX)/bin/clang" ] || missing="$$missing llvm@$(LLVM_PREFIX)"; \
	if [ -n "$$missing" ]; then echo "❌ missing host tools:$$missing"; exit 1; fi; \
	echo "✓ all host tools present"

# fail if the 5 version numbers disagree, or format.current != the ABI version
check-version:
	@scripts/check-version.sh

# run the Rust tests
test:
	cargo test --all-features --manifest-path $(CRATE_MANIFEST)

# fail if a fresh build changes any committed quantion.h
header-check:
	@scripts/header-check.sh

# build for the machine you are on, no Docker
native:
	@scripts/native.sh

# build the host binary, then check all four languages agree. the fast loop
quick-check:
	@$(MAKE) --no-print-directory native
	@$(MAKE) --no-print-directory check-api

# fail if Rust, R, Python and JS return different numbers
check-api:
	@scripts/check-api.sh

# rewrite check-api/expected.txt from the Rust run
check-api-update:
	@UPDATE=1 scripts/check-api.sh

# check, build every platform, write manifest.json. writes artifacts/<version>/
release:
	@scripts/release.sh

# set the version everywhere, release, commit, tag, push. takes VERSION=X.Y.Z
publish:
	@scripts/publish.sh "$(VERSION)"

# delete cargo build folders. keeps artifacts/
clean:
	cargo clean --manifest-path $(CRATE_MANIFEST)
	rm -rf core/target-linux-amd64 core/target-linux-arm64 \
	       core/target-windows-amd64 core/target-windows-arm64

# list every target with its one-line description
help:
	@awk '/^# / { text = substr($$0, 3); next } \
	      /^[a-z][a-zA-Z0-9_-]*:/ && text { printf "  %-16s %s\n", substr($$1, 1, length($$1) - 1), text; text = "" } \
	      { if ($$0 !~ /^# /) text = "" }' $(MAKEFILE_LIST)
