import { copyFileSync, existsSync, mkdirSync, readFileSync, rmSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const wrapper = join(here, "..");
const version = JSON.parse(readFileSync(join(wrapper, "package.json"), "utf8")).version;
const artifacts = join(wrapper, "..", "..", "artifacts", version);

const platforms = [
  ["macos-arm64", "libquantion.dylib"],
  ["macos-x86_64", "libquantion.dylib"],
  ["linux-arm64", "libquantion.so"],
  ["linux-x86_64", "libquantion.so"],
  ["windows-x86_64", "libquantion.dll"],
];

const clean = process.argv.includes("--clean");

if (clean) {
  for (const [platform] of platforms) {
    rmSync(join(wrapper, "native", platform), { recursive: true, force: true });
  }
} else {
  for (const [platform, file] of platforms) {
    const source = join(artifacts, platform, file);
    if (!existsSync(source)) {
      throw new Error(`quantion: ${source} is missing. Run 'make all' at the repo root first.`);
    }
  }

  for (const [platform, file] of platforms) {
    const target = join(wrapper, "native", platform);
    rmSync(target, { recursive: true, force: true });
    mkdirSync(target, { recursive: true });
    copyFileSync(join(artifacts, platform, file), join(target, file));
  }
}
