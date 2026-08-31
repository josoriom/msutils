import { describe, test, expect, beforeAll, afterAll, beforeEach } from "@jest/globals";
import * as http from "http";
import type { AddressInfo } from "net";

jest.mock(
  "../../build/Release/quantion.node",
  () => ({
    bind: jest.fn(),
    parseIonSource: jest.fn(),
    planOpen: jest.fn(),
    planEic: jest.fn(),
    planScans: jest.fn(),
    planImage: jest.fn(),
    calculateEic: jest.fn(),
    getScans: jest.fn(),
    getIonImage: jest.fn(),
    dispose: jest.fn(),
  }),
  { virtual: true },
);

const modulePaths: Array<[string, string]> = [
  ["source", "../utilities/backendNode"],
  ["built", "../../lib/utilities/backendNode.js"],
];

const EMPTY_EIC = { x: new Float64Array(), y: new Float64Array() };

function emptyBridge(
  payloadKind: number,
  sections: Array<[id: number, elementType: number, elementCount: number]>,
): ArrayBuffer {
  const headerBytes = 32;
  const entryBytes = 24;
  const sizes: Record<number, number> = { 1: 8, 2: 4, 3: 8, 4: 1 };
  let at = headerBytes + sections.length * entryBytes;
  const placed = sections.map(([id, elementType, elementCount]) => {
    const length = elementCount * sizes[elementType];
    const offset = at;
    at = (offset + length + 7) & ~7;
    return { id, elementType, offset, length };
  });
  const total = at;
  const buffer = new ArrayBuffer(total);
  const view = new DataView(buffer);
  view.setUint32(0, 0x42544e51, true);
  view.setUint16(4, 1, true);
  view.setUint16(6, payloadKind, true);
  view.setUint32(8, placed.length, true);
  view.setUint32(12, headerBytes, true);
  view.setBigUint64(16, BigInt(total), true);
  view.setBigUint64(24, 0n, true);
  placed.forEach((section, index) => {
    const entry = headerBytes + index * entryBytes;
    view.setUint32(entry, section.id, true);
    view.setUint32(entry + 4, section.elementType, true);
    view.setBigUint64(entry + 8, BigInt(section.offset), true);
    view.setBigUint64(entry + 16, BigInt(section.length), true);
  });
  return buffer;
}

const EMPTY_SCANS = emptyBridge(1, [
  [1, 3, 1],
  [2, 1, 0],
  [3, 1, 0],
  ...[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14].map(
    (id) => [id, 1, 0] as [number, number, number],
  ),
]);
const EMPTY_IMAGE = emptyBridge(2, [
  [15, 2, 6],
  [16, 1, 0],
  [17, 2, 0],
]);


const SAVED_LIB = process.env.QUANTION_LIB;
process.env.QUANTION_LIB = SAVED_LIB ?? "/quantion/mocked/native/library";

describe.each(modulePaths)("remote handle lifetime (%s)", (_name, modulePath) => {
  const { NodeBackend } = require(modulePath);

  let server: http.Server;
  let origin: string;
  let backend: any;
  let native: any;
  let disposedWhileStillTracked: boolean | null;

  beforeAll((done) => {
    server = http.createServer((_request, response) => {
      const bytes = Buffer.alloc(1024, 7);
      response.writeHead(206, {
        "content-length": bytes.length,
        "content-range": `bytes 0-1023/${bytes.length}`,
      });
      response.end(bytes);
    });
    server.listen(0, () => {
      const port = (server.address() as AddressInfo).port;
      origin = `http://127.0.0.1:${port}`;
      done();
    });
  });

  afterAll((done) => {
    server.close(() => done());
  });

  beforeEach(() => {
    disposedWhileStillTracked = null;
    backend = new NodeBackend(process);
    native = {
      parseIonSource: jest.fn(() => ({ handle: "native file" })),
      planOpen: jest.fn(() => []),
      planEic: jest.fn(() => []),
      planScans: jest.fn(() => []),
      planImage: jest.fn(() => []),
      calculateEic: jest.fn(() => EMPTY_EIC),
      getScans: jest.fn(() => EMPTY_SCANS),
      getIonImage: jest.fn(() => EMPTY_IMAGE),
      dispose: jest.fn((handle: object) => {
        disposedWhileStillTracked = backend.remote_by_handle.has(handle);
      }),
    };
    backend.native = native;
  });

  async function openRemote(): Promise<object> {
    return backend.parseIonRemote(new URL(`${origin}/file.ion`), 0);
  }

  test("the remote store is a weak map, so nothing pins a handle", () => {
    expect(backend.remote_by_handle).toBeInstanceOf(WeakMap);
    expect(backend.remote_by_handle).not.toBeInstanceOf(Map);
  });

  test("a numeric handle is rejected by the weak store", async () => {
    native.parseIonSource = jest.fn(() => 42);
    await expect(openRemote()).rejects.toThrow(/weak/i);
  });

  test("an open remote file is tracked and prefetches on eic", async () => {
    const handle = await openRemote();

    expect(backend.remote_by_handle.has(handle)).toBe(true);

    await backend.calculateEic(handle, 100, 0, 1, 5, 0.01);
    expect(native.planEic).toHaveBeenCalledTimes(1);
  });

  test("an open remote file prefetches on scans", async () => {
    const handle = await openRemote();

    await backend.getScans(handle, 0, 0, 10, 1);

    expect(native.planScans).toHaveBeenCalledTimes(1);
    expect(native.getScans).toHaveBeenCalledTimes(1);
  });

  test("an open remote file prefetches on ion image", async () => {
    const handle = await openRemote();

    await backend.getIonImage(handle, 100, 0.1, 1);

    expect(native.planImage).toHaveBeenCalledTimes(1);
    expect(native.getIonImage).toHaveBeenCalledTimes(1);
  });

  test("a local file asks for no scan plan at all", async () => {
    const handle = { handle: "local file" };

    await backend.getScans(handle, 0, 0, 10, 1);
    await backend.getIonImage(handle, 100, 0.1, 1);

    expect(native.planScans).not.toHaveBeenCalled();
    expect(native.planImage).not.toHaveBeenCalled();
    expect(native.getScans).toHaveBeenCalledTimes(1);
    expect(native.getIonImage).toHaveBeenCalledTimes(1);
  });

  test("freeFile drops the entry before it disposes the handle", async () => {
    const handle = await openRemote();

    backend.freeFile(handle);

    expect(disposedWhileStillTracked).toBe(false);
    expect(native.dispose).toHaveBeenCalledWith(handle);
    expect(backend.remote_by_handle.has(handle)).toBe(false);
  });

  test("freeFile releases the byte cache, so later reads stop prefetching", async () => {
    const handle = await openRemote();
    backend.freeFile(handle);

    await backend.calculateEic(handle, 100, 0, 1, 5, 0.01);

    expect(native.planEic).not.toHaveBeenCalled();
    expect(native.calculateEic).toHaveBeenCalledTimes(1);
  });

  test("freeing one remote file leaves another one tracked", async () => {
    const first = await openRemote();
    native.parseIonSource = jest.fn(() => ({ handle: "second file" }));
    const second = await openRemote();

    backend.freeFile(first);

    expect(backend.remote_by_handle.has(first)).toBe(false);
    expect(backend.remote_by_handle.has(second)).toBe(true);
  });
});
