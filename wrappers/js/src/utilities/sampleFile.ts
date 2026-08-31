import type { Backend, FileHandle } from "./backend";
import { fileHandleRegistry } from "./handleRegistry";

type Cell = ReturnType<typeof fileHandleRegistry.register>;

/**
 * A loaded sample file. Holds an internal handle to the native
 * file object. Cleanup is automatic, but call {@link dispose} when
 * you're done with large samples or tight processing loops to free memory sooner.
 */
export class SampleFile {
  /** @internal */
  _backend: Backend;

  private _cell: Cell | null;

  constructor(handle: FileHandle, backend: Backend) {
    this._backend = backend;

    if (backend.handlesAreGcFinalized) {
      this._cell = { handle, free: () => {} };
    } else {
      this._cell = fileHandleRegistry.register(this, handle, (h) =>
        backend.freeFile(h),
      );
    }
  }

  get _handle(): FileHandle | null {
    return this._cell?.handle ?? null;
  }

  /**
   * Release the native memory held by this sample.
   */
  dispose(): void {
    const cell = this._cell;
    if (!cell) return;

    if (this._backend.handlesAreGcFinalized) {
      const h = cell.handle;
      if (h != null) {
        this._backend.freeFile(h);
        cell.handle = null;
      }
    } else {
      fileHandleRegistry.release(this, cell);
    }
  }

  /**
   * Serialize the sample back to mzML format.
   * @returns Full mzML file content as a string.
   */
  toMzml(): string {
    this.assertValid("toMzml");
    return this._backend.fileToMzml(this._handle!);
  }

  /**
   * Encode the sample as compressed ion binary bytes.
   *
   * @param options - Encoding options.
   * @param options.level - Compression level, 0 (none) to 22 (max). Default 12.
   * @param options.f32Compress - Compress intensity values to 32-bit float. Default false.
   * @returns Raw ion binary bytes.
   */
  toIon(options: { level?: number; f32Compress?: boolean } = {}): Uint8Array {
    const { level = 12, f32Compress = false } = options;
    this.assertValid("toIon");
    return this._backend.fileToBin(this._handle!, level, f32Compress);
  }

  private assertValid(caller: string): void {
    if (!this._handle) {
      throw new Error(`${caller}: file has been disposed`);
    }
  }
}
