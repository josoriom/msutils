import { WasmBackend } from "./utilities/backendWasm";
import { setBackend } from "./utilities/api";

const backend = new WasmBackend();
const _ready = backend.init().then(() => setBackend(backend));

export async function init(): Promise<void> {
  await _ready;
}

export * from "./utilities/api";
export { SampleFile } from "./utilities/sampleFile";
