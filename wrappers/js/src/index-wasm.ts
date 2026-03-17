import { WasmBackend } from "./utilities/backendWasm";
import { setBackend } from "./utilities/api";

const backend = new WasmBackend();

export async function init(): Promise<void> {
  await backend.init();
  setBackend(backend);
}

export * from "./utilities/api";
export { MzMlFile } from "./utilities/mzmlFile";
