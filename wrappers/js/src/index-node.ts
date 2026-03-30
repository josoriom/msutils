import { NodeBackend } from "./utilities/backendNode";
import { setBackend, setInitPromise } from "./utilities/api";

const backend = new NodeBackend(process);
const _ready = Promise.resolve().then(() => setBackend(backend));
setInitPromise(_ready);

export async function init(): Promise<void> {
  await _ready;
}

export * from "./utilities/api";
export { MzMlFile } from "./utilities/mzmlFile";
export { packPeakOptions } from "./utilities/pack";
