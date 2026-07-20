import { copyFileSync, existsSync, mkdirSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const wrapper = join(here, "..");
const version = JSON.parse(readFileSync(join(wrapper, "package.json"), "utf8")).version;
const source = join(wrapper, "..", "..", "artifacts", version, "wasm", "quantion.wasm");

if (!existsSync(source)) {
  throw new Error(`quantion: ${source} is missing. Run 'make wasm' at the repo root first.`);
}
mkdirSync(join(wrapper, "dist"), { recursive: true });
copyFileSync(source, join(wrapper, "dist", "quantion.wasm"));
