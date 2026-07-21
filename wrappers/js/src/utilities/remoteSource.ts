import { ContainmentCache } from "./containmentCache";

export type ByteRangeResult = { offset: bigint; length: bigint };

export type RemoteSource = {
  url: string;
  cache: ContainmentCache;
  total: bigint;
};

export const HEADER_BYTES = 1024n;

const LARGEST_GAP = 65536n;

function gapFor(total: bigint): bigint {
  if (total <= 0n) return LARGEST_GAP;
  const eighth = total / 8n;
  return eighth < LARGEST_GAP ? eighth : LARGEST_GAP;
}

export function newRemoteSource(url: string): RemoteSource {
  return { url, cache: new ContainmentCache(), total: 0n };
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
  if (!contentRange?.startsWith(expectedPrefix)) {
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

export async function prefetchRanges(
  source: RemoteSource,
  ranges: ByteRangeResult[],
): Promise<void> {
  const missing = source.cache.missing(ranges);
  const wanted = coalesceRanges(missing, gapFor(source.total));

  await Promise.all(
    wanted.map(async (range) => {
      const bytes = await fetchRange(source.url, range);
      source.cache.add(range, bytes);
    }),
  );
}

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

  const contentRange = response.headers.get("Content-Range") ?? "";
  const total = contentRange.split("/")[1] ?? "";
  if (!/^\d+$/.test(total)) {
    throw new Error(
      `range response for ${source.url} has no total size: Content-Range is ${contentRange}`,
    );
  }
  source.total = BigInt(total);

  const bytes = new Uint8Array(await response.arrayBuffer());
  const range: ByteRangeResult = { offset: 0n, length: BigInt(bytes.length) };
  source.cache.add(range, bytes);
  return bytes;
}
