export interface ByteRange {
  offset: bigint;
  length: bigint;
}

interface CachedRange {
  start: bigint;
  end: bigint;
  bytes: Uint8Array;
  pinned: boolean;
}

const COALESCE_GAP_THRESHOLD = 65536n;

export class ContainmentCache {
  private cached: CachedRange[];

  constructor() {
    this.cached = [];
  }

  add(range: ByteRange, bytes: Uint8Array): void {
    if (BigInt(bytes.length) !== range.length) {
      throw new Error(
        `cache add: byte length ${bytes.length} does not match range.length ${range.length}`
      );
    }

    const new_start = range.offset;
    const new_end = range.offset + range.length;

    this.cached.push({ start: new_start, end: new_end, bytes, pinned: false });
    this.cached.sort((a, b) => (a.start < b.start ? -1 : 1));
  }

  read(offset: bigint, length: bigint): Uint8Array | null {
    const target_end = offset + length;

    for (const cached of this.cached) {
      if (cached.start <= offset && target_end <= cached.end) {
        const slice_start = Number(offset - cached.start);
        const slice_end = slice_start + Number(length);
        return cached.bytes.slice(slice_start, slice_end);
      }
    }
    return null;
  }

  has(range: ByteRange): boolean {
    return this.read(range.offset, range.length) !== null;
  }

  missing(ranges: ByteRange[]): ByteRange[] {
    return ranges.filter((range) => !this.has(range));
  }

  clear(): void {
    this.cached = [];
  }

  pin(): void {
    for (const cached of this.cached) cached.pinned = true;
  }

  dropUnpinned(): void {
    this.cached = this.cached.filter((cached) => cached.pinned);
  }

  bytes(): number {
    let total = 0;
    for (const cached of this.cached) total += cached.bytes.length;
    return total;
  }
}
