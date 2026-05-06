import type { Backend, FileHandle } from "./backend";
import { camelizeKeys } from "./shared";

/**
 * A loaded mass spectrometry sample. Holds an internal handle to the native
 * file object. Call {@link dispose} or use a try/finally block to release
 * memory when done.
 */
export class SampleFile {
  /** @internal */
  _handle: FileHandle | null;
  /** @internal */
  _backend: Backend;

  constructor(handle: FileHandle, backend: Backend) {
    this._handle = handle;
    this._backend = backend;
  }

  /**
   * Release the native memory held by this sample. Safe to call more than once.
   * After disposal all other methods throw.
   */
  dispose(): void {
    if (this._handle) {
      this._backend.freeFile(this._handle);
      this._handle = null;
    }
  }

  /**
   * Return the sample data as a plain JavaScript object.
   *
   * @returns Parsed sample data with camelCase keys.
   */
  toJson(): any {
    this.assertValid("toJson");
    return camelizeKeys(this._backend.fileToJson(this._handle!));
  }

  /**
   * Serialize the sample back to mzML format.
   *
   * @returns Full mzML file content as a string.
   */
  toMzml(): string {
    this.assertValid("toMzml");
    return this._backend.fileToMzml(this._handle!);
  }

  /**
   * Encode the sample as compressed ion binary bytes.
   *
   * @param level - Compression level, 0 (none) to 22 (max). Default 12.
   * @param f32Compress - Compress intensity values to 32-bit float. Default false.
   * @returns Raw ion binary bytes.
   */
  toIon(level = 12, f32Compress = false): Uint8Array {
    this.assertValid("toIon");
    return this._backend.fileToBin(this._handle!, level, f32Compress);
  }

  private assertValid(caller: string): void {
    if (!this._handle) {
      throw new Error(`${caller}: file has been disposed`);
    }
  }
}
