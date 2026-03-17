import type { Backend, FileHandle } from "./backend";
import { camelizeKeys } from "./shared";

export class MzMlFile {
  /** @internal */
  _handle: FileHandle | null;
  /** @internal */
  _backend: Backend;

  constructor(handle: FileHandle, backend: Backend) {
    this._handle = handle;
    this._backend = backend;
  }

  dispose(): void {
    if (this._handle) {
      this._backend.freeFile(this._handle);
      this._handle = null;
    }
  }

  toJson(): any {
    this.assertValid("toJson");
    return camelizeKeys(this._backend.fileToJson(this._handle!));
  }

  toMzml(): string {
    this.assertValid("toMzml");
    return this._backend.fileToMzml(this._handle!);
  }

  toBin(level = 12, f32Compress = false): Uint8Array {
    this.assertValid("toBin");
    return this._backend.fileToBin(this._handle!, level, f32Compress);
  }

  private assertValid(caller: string): void {
    if (!this._handle) {
      throw new Error(`${caller}: file has been disposed`);
    }
  }
}
