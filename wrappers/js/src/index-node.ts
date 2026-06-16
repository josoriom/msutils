import { NodeBackend } from "./utilities/backendNode";
import { setBackend, setInitPromise } from "./utilities/api";

const backend = new NodeBackend(process);
setBackend(backend);

const _ready = Promise.resolve();
setInitPromise(_ready);

export async function init(): Promise<void> {
  await _ready;
}

export * from "./utilities/api";
export { SampleFile } from "./utilities/sampleFile";
export { encodeTargetIds, unpackTargets } from "./utilities/pack";
export { ContainmentCache, type ByteRange } from "./utilities/containmentCache";
