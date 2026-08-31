import { ContainmentCache } from "./containmentCache";

export type ByteRangeResult = { offset: bigint; length: bigint };

export type RemoteSource = {
  url: string;
  cache: ContainmentCache;
  total: bigint;
  multipart: boolean | null;
};

export const HEADER_BYTES = 1024n;

const LARGEST_GAP = 131072n;

function gapFor(total: bigint): bigint {
  const eighth = total / 8n;
  return eighth < LARGEST_GAP ? eighth : LARGEST_GAP;
}

const MAX_RANGES_PER_REQUEST = 100;

const MULTIPART_MIN_RANGES = 8;

const PART_SIZE = 16n * 1024n * 1024n;

const MAX_PARTS_IN_FLIGHT = 16;

export function newRemoteSource(url: string): RemoteSource {
  return { url, cache: new ContainmentCache(), total: 0n, multipart: null };
}

export async function fetchRange(
  url: string,
  range: ByteRangeResult,
): Promise<Uint8Array> {
  const start = range.offset;
  const end = range.offset + range.length - 1n;

  const response = await fetch(url, {
    headers: {
      Range: `bytes=${start}-${end}`,
      "Accept-Encoding": "identity",
    },
  });

  if (response.status !== 206) {
    throw new Error(
      `range request failed: expected 206 Partial Content, got ${response.status}`,
    );
  }

  const contentRange = response.headers.get("Content-Range");
  const expectedPrefix = `bytes ${start}-${end}/`;
  if (contentRange !== null && !contentRange.startsWith(expectedPrefix)) {
    throw new Error(
      `range response has wrong Content-Range: got ${contentRange}, expected to start with ${expectedPrefix}`,
    );
  }

  const bytes = new Uint8Array(await response.arrayBuffer());
  if (BigInt(bytes.length) !== range.length) {
    throw new Error(
      `range response length mismatch: got ${bytes.length}, expected ${range.length}`,
    );
  }

  return bytes;
}

export function coalesceRanges(
  ranges: ByteRangeResult[],
  gap: bigint = LARGEST_GAP,
): ByteRangeResult[] {
  if (ranges.length === 0) return [];

  const sorted = [...ranges].sort((a, b) =>
    a.offset < b.offset ? -1 : a.offset > b.offset ? 1 : 0,
  );

  const result: ByteRangeResult[] = [];
  let current = sorted[0];

  for (let index = 1; index < sorted.length; index++) {
    const next = sorted[index];
    if (next.offset - (current.offset + current.length) <= gap) {
      current = {
        offset: current.offset,
        length: next.offset + next.length - current.offset,
      };
    } else {
      result.push(current);
      current = next;
    }
  }
  result.push(current);
  return result;
}

function indexOfBytes(
  haystack: Uint8Array,
  needle: Uint8Array,
  from: number,
): number {
  const limit = haystack.length - needle.length;
  outer: for (let start = from; start <= limit; start++) {
    for (let index = 0; index < needle.length; index++) {
      if (haystack[start + index] !== needle[index]) continue outer;
    }
    return start;
  }
  return -1;
}

function boundaryOf(contentType: string): string | null {
  const match = /;\s*boundary=("?)([^";]+)\1/i.exec(contentType);
  return match ? match[2] : null;
}

function offsetOfPart(headers: string): bigint | null {
  const match = /^content-range:\s*bytes\s+(\d+)-(\d+)\//im.exec(headers);
  return match ? BigInt(match[1]) : null;
}

/**
 * Split a `multipart/byteranges` body into the parts it carries. Each part is
 * located by its own `Content-Range`, so the parts need not arrive in the order
 * they were asked for. Returns null when the body is not parseable as multipart.
 */
export function splitByteRanges(
  body: Uint8Array,
  boundary: string,
): { offset: bigint; bytes: Uint8Array }[] | null {
  const ascii = new TextDecoder("latin1");
  const delimiter = new TextEncoder().encode(`--${boundary}`);
  const terminator = new TextEncoder().encode("\r\n\r\n");

  const parts: { offset: bigint; bytes: Uint8Array }[] = [];
  let cursor = indexOfBytes(body, delimiter, 0);
  if (cursor < 0) return null;

  while (cursor >= 0) {
    cursor += delimiter.length;
    if (body[cursor] === 0x2d && body[cursor + 1] === 0x2d) break;

    const headerEnd = indexOfBytes(body, terminator, cursor);
    if (headerEnd < 0) return null;

    const headers = ascii.decode(body.subarray(cursor, headerEnd));
    const offset = offsetOfPart(headers);
    if (offset === null) return null;

    const bodyStart = headerEnd + terminator.length;
    const next = indexOfBytes(body, delimiter, bodyStart);
    if (next < 0) return null;

    parts.push({ offset, bytes: body.slice(bodyStart, next - 2) });
    cursor = next;
  }

  return parts.length > 0 ? parts : null;
}

/**
 * Ask for every range in a single request. Resolves to null when the server
 * does not honour multi-range requests, which object stores such as S3, GCS and
 * R2 signal by answering 200 with the whole object; the body is discarded
 * unread in that case so nothing extra is downloaded.
 */
export async function fetchRanges(
  url: string,
  ranges: ByteRangeResult[],
): Promise<Uint8Array[] | null> {
  const spec = ranges
    .map((range) => `${range.offset}-${range.offset + range.length - 1n}`)
    .join(", ");

  const response = await fetch(url, {
    headers: {
      Range: `bytes=${spec}`,
      "Accept-Encoding": "identity",
    },
  });

  if (response.status !== 206) {
    await response.body?.cancel();
    return null;
  }

  const boundary = boundaryOf(response.headers.get("Content-Type") ?? "");
  if (boundary === null) {
    await response.body?.cancel();
    return null;
  }

  const parts = splitByteRanges(
    new Uint8Array(await response.arrayBuffer()),
    boundary,
  );
  if (parts === null) return null;

  const ordered: Uint8Array[] = [];
  for (const range of ranges) {
    const end = range.offset + range.length;
    const holder = parts.find(
      (part) =>
        part.offset <= range.offset &&
        part.offset + BigInt(part.bytes.length) >= end,
    );
    if (!holder) return null;
    const start = Number(range.offset - holder.offset);
    ordered.push(holder.bytes.subarray(start, start + Number(range.length)));
  }
  return ordered;
}

export async function fetchRangeParts(
  url: string,
  range: ByteRangeResult,
): Promise<Uint8Array> {
  const count = Number((range.length + PART_SIZE - 1n) / PART_SIZE);
  const bytes = new Uint8Array(Number(range.length));
  let next = 0;

  while (next < count) {
    const batch = [];
    for (
      let index = next;
      index < Math.min(next + MAX_PARTS_IN_FLIGHT, count);
      index++
    ) {
      const offset = range.offset + BigInt(index) * PART_SIZE;
      const length =
        offset + PART_SIZE > range.offset + range.length
          ? range.offset + range.length - offset
          : PART_SIZE;
      const at = Number(offset - range.offset);
      batch.push(
        fetchRange(url, { offset, length }).then((part) => bytes.set(part, at)),
      );
    }
    await Promise.all(batch);
    next += batch.length;
  }

  return bytes;
}

function chunked<T>(items: T[], size: number): T[][] {
  const chunks: T[][] = [];
  for (let start = 0; start < items.length; start += size) {
    chunks.push(items.slice(start, start + size));
  }
  return chunks;
}

async function fetchEachRange(
  source: RemoteSource,
  ranges: ByteRangeResult[],
): Promise<void> {
  const wanted = coalesceRanges(ranges, gapFor(source.total));
  await Promise.all(
    wanted.map(async (range) => {
      const bytes =
        range.length > PART_SIZE
          ? await fetchRangeParts(source.url, range)
          : await fetchRange(source.url, range);
      source.cache.add(range, bytes);
    }),
  );
}

export async function prefetchRanges(
  source: RemoteSource,
  ranges: ByteRangeResult[],
): Promise<void> {
  const missing = source.cache.missing(ranges);
  if (missing.length === 0) return;

  const coalesced = coalesceRanges(missing, gapFor(source.total));

  const large = coalesced.filter((range) => range.length > PART_SIZE);
  await Promise.all(
    large.map(async (range) => {
      const bytes = await fetchRangeParts(source.url, range);
      source.cache.add(range, bytes);
    }),
  );

  const exact = coalesced.filter((range) => range.length <= PART_SIZE);
  if (exact.length === 0) return;
  if (exact.length === 1) {
    const bytes = await fetchRange(source.url, exact[0]);
    source.cache.add(exact[0], bytes);
    return;
  }

  if (source.multipart !== false && exact.length >= MULTIPART_MIN_RANGES) {
    const batches = chunked(exact, MAX_RANGES_PER_REQUEST);
    const results = await Promise.all(
      batches.map((batch) => fetchRanges(source.url, batch)),
    );

    if (results.every((result) => result !== null)) {
      source.multipart = true;
      batches.forEach((batch, index) => {
        const bytes = results[index] as Uint8Array[];
        batch.forEach((range, part) => source.cache.add(range, bytes[part]));
      });
      return;
    }

    source.multipart = false;
    const done = new Set<number>();
    results.forEach((result, index) => {
      if (result === null) return;
      done.add(index);
      batches[index].forEach((range, part) =>
        source.cache.add(range, result[part]),
      );
    });

    const rest = batches.filter((_, index) => !done.has(index)).flat();
    await fetchEachRange(source, rest);
    return;
  }

  await fetchEachRange(source, exact);
}

function readTotalFromContentRange(response: Response): bigint | null {
  const total =
    (response.headers.get("Content-Range") ?? "").split("/")[1] ?? "";
  return /^\d+$/.test(total) ? BigInt(total) : null;
}

async function fetchTotalByHead(url: string): Promise<bigint> {
  const response = await fetch(url, { method: "HEAD" });
  if (!response.ok) {
    throw new Error(`size request failed for ${url}: got ${response.status}`);
  }

  const encoding = response.headers.get("Content-Encoding");
  if (encoding !== null && encoding !== "identity") {
    throw new Error(
      `size request for ${url} was ${encoding}-encoded, so Content-Length is not the file size`,
    );
  }

  const length = response.headers.get("Content-Length") ?? "";
  if (!/^\d+$/.test(length)) {
    throw new Error(
      `size request for ${url} returned no usable Content-Length`,
    );
  }

  return BigInt(length);
}

const SIGNATURE = [0x49, 0x4f, 0x4e, 0x49, 0x43];
const FILE_SIZE_OFFSET = 400;

export async function fetchHeader(source: RemoteSource): Promise<Uint8Array> {
  const response = await fetch(source.url, {
    headers: {
      Range: `bytes=0-${HEADER_BYTES - 1n}`,
      "Accept-Encoding": "identity",
    },
  });
  if (response.status !== 206) {
    throw new Error(
      `range request failed for ${source.url}: expected 206 Partial Content, got ${response.status}`,
    );
  }

  const bytes = new Uint8Array(await response.arrayBuffer());

  let totalFromHeader: bigint | null = null;
  if (bytes.length >= FILE_SIZE_OFFSET + 8) {
    let signed = true;
    for (let index = 0; index < SIGNATURE.length; index++) {
      if (bytes[index] !== SIGNATURE[index]) signed = false;
    }
    if (signed) {
      const view = new DataView(
        bytes.buffer,
        bytes.byteOffset,
        bytes.byteLength,
      );
      const stored = view.getBigUint64(FILE_SIZE_OFFSET, true);
      if (stored > 0n) totalFromHeader = stored;
    }
  }

  source.total =
    totalFromHeader ??
    readTotalFromContentRange(response) ??
    (await fetchTotalByHead(source.url));

  const range: ByteRangeResult = { offset: 0n, length: BigInt(bytes.length) };
  source.cache.add(range, bytes);
  return bytes;
}
