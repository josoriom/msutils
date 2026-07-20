import os
import sys
import shutil
import platform
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py


def detect_platform() -> str:
    system  = sys.platform
    machine = platform.machine().lower()

    os_name = {"darwin": "macos", "linux": "linux", "win32": "windows"}.get(system)
    if not os_name:
        raise RuntimeError(f"Unsupported OS: {system}")

    if machine in ("arm64", "aarch64"):
        arch = "arm64"
    elif machine in ("amd64", "x86_64"):
        arch = "x86_64"
    else:
        raise RuntimeError(f"Unsupported architecture: {machine}")

    return f"{os_name}-{arch}"


def copy_native_binary(platform_dir: str, dst_root: Path) -> None:
    here = Path(__file__).parent.resolve()
    artifacts_root = (here / ".." / ".." / "artifacts").resolve()
    src = artifacts_root / platform_dir

    if not src.is_dir() and artifacts_root.is_dir():
        versions = sorted(
            (p for p in artifacts_root.iterdir() if (p / platform_dir).is_dir()),
            key=lambda p: tuple(int(part) if part.isdigit() else 0 for part in p.name.split(".")),
            reverse=True,
        )
        if versions:
            src = versions[0] / platform_dir

    if not src.is_dir():
        raise FileNotFoundError(
            f"\n"
            f"  Native binary not found at: {src}\n"
            f"\n"
            f"  Run `make {platform_dir}` at the repo root first, then reinstall:\n"
            f"\n"
            f"    pip install --force-reinstall .\n"
            f"  or\n"
            f"    pip install --force-reinstall "
            f"git+https://github.com/josoriom/quantion#subdirectory=wrappers/python\n"
        )

    dst = dst_root / "quantion" / "native" / platform_dir
    dst.mkdir(parents=True, exist_ok=True)

    copied = []
    for f in src.iterdir():
        if f.is_file():
            shutil.copy2(f, dst / f.name)
            copied.append(f.name)

    if not copied:
        raise FileNotFoundError(f"No files found in {src} — did `make` finish?")

    print(f"[quantion] copied {platform_dir}: {', '.join(copied)}")


class BuildWithNative(build_py):
    def run(self):
        platform_dir = detect_platform()
        build_root = Path(self.build_lib).resolve()
        copy_native_binary(platform_dir, build_root)
        super().run()


setup(cmdclass={"build_py": BuildWithNative})