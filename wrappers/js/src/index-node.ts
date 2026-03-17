import { NodeBackend } from "./utilities/backendNode";
import { setBackend } from "./utilities/api";

const backend = new NodeBackend(process);
setBackend(backend);

export * from "./utilities/api";
export { MzMlFile } from "./utilities/mzmlFile";
export { packPeakOptions } from "./utilities/pack";
