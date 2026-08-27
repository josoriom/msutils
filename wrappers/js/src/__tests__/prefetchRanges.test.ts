import { describe, test, expect, beforeEach } from "@jest/globals";
import {
  coalesceRanges,
  fetchRangeParts,
  newRemoteSource,
  prefetchRanges,
  splitByteRanges,
  type ByteRangeResult,
} from "../utilities/remoteSource";

const BOUNDARY = "QUANTION";

type Request = { start: number; end: number };

function body(source: Uint8Array, spans: Request[]): Uint8Array {
  const encoder = new TextEncoder();
  const chunks: Uint8Array[] = [];
  for (const span of spans) {
    chunks.push(
      encoder.encode(
        `--${BOUNDARY}\r\nContent-Type: application/octet-stream\r\n` +
          `Content-Range: bytes ${span.start}-${span.end}/${source.length}\r\n\r\n`,
      ),
    );
    chunks.push(source.subarray(span.start, span.end + 1));
    chunks.push(encoder.encode("\r\n"));
  }
  chunks.push(encoder.encode(`--${BOUNDARY}--\r\n`));

  const total = chunks.reduce((sum, chunk) => sum + chunk.length, 0);
  const out = new Uint8Array(total);
  let at = 0;
  for (const chunk of chunks) {
    out.set(chunk, at);
    at += chunk.length;
  }
  return out;
}

function parseSpec(header: string): Request[] {
  return header
    .replace("bytes=", "")
    .split(",")
    .map((part) => part.trim())
    .map((part) => {
      const [start, end] = part.split("-");
      return { start: Number(start), end: Number(end) };
    });
}

function serve(
  data: Uint8Array,
  mode: "multipart" | "merged" | "missing" | "whole",
) {
  const calls: string[] = [];
  const fetcher = jest.fn(async (_url: string, init: any) => {
    const header = init.headers.Range as string;
    calls.push(header);
    const asked = parseSpec(header);

    if (asked.length === 1) {
      const span = asked[0];
      return {
        status: 206,
        headers: {
          get: (name: string) =>
            name === "Content-Range"
              ? `bytes ${span.start}-${span.end}/${data.length}`
              : null,
        },
        arrayBuffer: async () =>
          data.subarray(span.start, span.end + 1).slice().buffer,
        body: null,
      };
    }

    if (mode === "whole") {
      return {
        status: 200,
        headers: { get: () => null },
        arrayBuffer: async () => data.slice().buffer,
        body: { cancel: async () => undefined },
      };
    }

    const served =
      mode === "merged"
        ? [{ start: asked[0].start, end: asked[asked.length - 1].end }]
        : mode === "missing"
          ? [asked[0]]
          : asked;

    const payload = body(data, served);
    return {
      status: 206,
      headers: {
        get: (name: string) =>
          name === "Content-Type"
            ? `multipart/byteranges; boundary=${BOUNDARY}`
            : null,
      },
      arrayBuffer: async () => payload.slice().buffer,
      body: { cancel: async () => undefined },
    };
  });

  return { fetcher, calls };
}

function sample(length: number): Uint8Array {
  const data = new Uint8Array(length);
  for (let index = 0; index < length; index++) data[index] = index % 251;
  return data;
}

describe("prefetchRanges", () => {
  const original = globalThis.fetch;

  beforeEach(() => {
    globalThis.fetch = original;
  });

  test("coalesces at the 128 KB gap", () => {
    const ranges: ByteRangeResult[] = [
      { offset: 0n, length: 1024n },
      { offset: 100000n, length: 1024n },
      { offset: 1000000n, length: 1024n },
    ];
    expect(coalesceRanges(ranges).length).toBe(2);
    expect(coalesceRanges(ranges, 0n).length).toBe(3);
  });

  test("a single-range plan never touches multipart", async () => {
    const data = sample(4096);
    const { fetcher, calls } = serve(data, "multipart");
    globalThis.fetch = fetcher as any;

    const source = newRemoteSource("http://x/f.ion");
    await prefetchRanges(source, [
      { offset: 0n, length: 512n },
      { offset: 512n, length: 512n },
    ]);

    expect(calls).toEqual(["bytes=0-1023"]);
    expect(source.multipart).toBeNull();
    expect(source.cache.read(0n, 1024n)).toEqual(data.subarray(0, 1024));
  });

  test("a server that merges two close ranges satisfies both with no refetch", async () => {
    const data = sample(8192);
    const { fetcher, calls } = serve(data, "merged");
    globalThis.fetch = fetcher as any;

    const source = newRemoteSource("http://x/f.ion");
    await prefetchRanges(source, [
      { offset: 0n, length: 512n },
      { offset: 4096n, length: 512n },
    ]);

    expect(calls.length).toBe(1);
    expect(source.multipart).toBe(true);
    expect(source.cache.read(0n, 512n)).toEqual(data.subarray(0, 512));
    expect(source.cache.read(4096n, 512n)).toEqual(data.subarray(4096, 4608));
  });

  test("a part genuinely missing a range falls back", async () => {
    const data = sample(8192);
    const { fetcher, calls } = serve(data, "missing");
    globalThis.fetch = fetcher as any;

    const source = newRemoteSource("http://x/f.ion");
    await prefetchRanges(source, [
      { offset: 0n, length: 512n },
      { offset: 4096n, length: 512n },
    ]);

    expect(source.multipart).toBe(false);
    expect(calls.length).toBeGreaterThan(1);
    expect(source.cache.read(0n, 512n)).toEqual(data.subarray(0, 512));
    expect(source.cache.read(4096n, 512n)).toEqual(data.subarray(4096, 4608));
  });

  test("a 200 whole-object answer cancels the body and falls back", async () => {
    const data = sample(8192);
    const { fetcher } = serve(data, "whole");
    globalThis.fetch = fetcher as any;

    const source = newRemoteSource("http://x/f.ion");
    await prefetchRanges(source, [
      { offset: 0n, length: 512n },
      { offset: 4096n, length: 512n },
    ]);

    expect(source.multipart).toBe(false);
    expect(source.cache.read(0n, 512n)).toEqual(data.subarray(0, 512));
    expect(source.cache.read(4096n, 512n)).toEqual(data.subarray(4096, 4608));
  });

  test("splitByteRanges reads parts by their own Content-Range", () => {
    const data = sample(4096);
    const parts = splitByteRanges(
      body(data, [
        { start: 2048, end: 2559 },
        { start: 0, end: 511 },
      ]),
      BOUNDARY,
    );
    expect(parts).not.toBeNull();
    expect(parts!.map((part) => part.offset)).toEqual([2048n, 0n]);
    expect(parts![1].bytes).toEqual(data.subarray(0, 512));
  });
});

describe("fetchRangeParts", () => {
  const original = globalThis.fetch;

  test("splits a large run and reassembles it byte-exact", async () => {
    const size = 40 * 1024 * 1024;
    const data = sample(size);
    const { fetcher, calls } = serve(data, "multipart");
    globalThis.fetch = fetcher as any;

    const bytes = await fetchRangeParts("http://x/f.ion", {
      offset: 0n,
      length: BigInt(size),
    });

    expect(calls.length).toBe(3);
    expect(bytes.length).toBe(size);
    expect(bytes).toEqual(data);

    globalThis.fetch = original;
  });

  test("prefetchRanges part-splits above 16 MiB", async () => {
    const size = 20 * 1024 * 1024;
    const data = sample(size);
    const { fetcher, calls } = serve(data, "multipart");
    globalThis.fetch = fetcher as any;

    const source = newRemoteSource("http://x/f.ion");
    await prefetchRanges(source, [{ offset: 0n, length: BigInt(size) }]);

    expect(calls.length).toBe(2);
    expect(source.cache.read(0n, BigInt(size))).toEqual(data);

    globalThis.fetch = original;
  });
});
