import { describe, test, expect, beforeAll, afterAll, beforeEach } from "@jest/globals";
import * as http from "http";
import type { AddressInfo } from "net";

jest.mock(
  "../../build/Release/quantion.node",
  () => ({
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
      getScans: jest.fn(() => "[]"),
      getIonImage: jest.fn(() => "{}"),
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
