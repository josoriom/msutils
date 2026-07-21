import { WasmBackend } from "./utilities/backendWasm";
import { setBackend, setInitPromise } from "./utilities/api";

const backend = new WasmBackend();
const _ready = backend.init().then(() => setBackend(backend));
setInitPromise(_ready);

export async function init(): Promise<void> {
  await _ready;
}

export * from "./utilities/api";
export { SampleFile } from "./utilities/sampleFile";
