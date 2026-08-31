import { describe, test, expect, beforeEach } from "@jest/globals";
import { parseIon, setBackend } from "../utilities/api";
import type { Backend, FileHandle } from "../utilities/backend";

describe("parseIon dispatch to backend primitives", () => {
  let mockBackend: Backend;

  beforeEach(() => {
    const parseIonPathMock = jest.fn(() => 42 as FileHandle);
    const parseIonRemoteMock = jest.fn(() => 43 as FileHandle);
    const parseIonBufferMock = jest.fn(() => 44 as FileHandle);

    mockBackend = {
      ready: true,
      handlesAreGcFinalized: true,
      init: jest.fn(),
      parseMzML: jest.fn(),
      parseBin: jest.fn(),
      parseIonPath: parseIonPathMock,
      parseIonRemote: parseIonRemoteMock,
      parseIonBuffer: parseIonBufferMock,
      freeFile: jest.fn(),
      fileToMzml: jest.fn(),
      fileToBin: jest.fn(),
      calculateEic: jest.fn(),
      getScans: jest.fn(),
      findPeaks: jest.fn(),
      getPeak: jest.fn(),
      findNoiseLevel: jest.fn(),
      calculateBaseline: jest.fn(),
      getPeaksFromEic: jest.fn(),
      getPeaksFromChrom: jest.fn(),
      findFeatures: jest.fn(),
      findFeature: jest.fn(),
    } as any;

    setBackend(mockBackend);
  });

  describe("path string dispatch", () => {
    test("parseIon with path string calls parseIonPath only", async () => {
      const path = "/data/file.ion";
      await parseIon(path);

      expect(mockBackend.parseIonPath).toHaveBeenCalledWith(path, 0);
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });

    test("parseIon with path string and cache size calls parseIonPath with cache size", async () => {
      const path = "/data/file.ion";
      const cacheSize = 1000;
      await parseIon(path, { maxCacheSize: cacheSize });

      expect(mockBackend.parseIonPath).toHaveBeenCalledWith(path, cacheSize);
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });
  });

  describe("URL dispatch", () => {
    test("parseIon with https URL calls parseIonRemote only", async () => {
      const url = new URL("https://example.com/file.ion");
      await parseIon(url);

      expect(mockBackend.parseIonRemote).toHaveBeenCalledWith(url, 0);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });

    test("parseIon with http URL calls parseIonRemote only", async () => {
      const url = new URL("http://example.com/file.ion");
      await parseIon(url);

      expect(mockBackend.parseIonRemote).toHaveBeenCalledWith(url, 0);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });

    test("parseIon rejects a file: URL", async () => {
      await expect(parseIon(new URL("file:///data/file.ion"))).rejects.toThrow(
        /only http and https/,
      );
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });

    test("parseIon accepts an https URL string and routes it to parseIonRemote", async () => {
      await parseIon("https://example.com/file.ion");

      const calls = (mockBackend.parseIonRemote as jest.Mock).mock.calls;
      expect(calls.length).toBe(1);
      expect(calls[0][0].href).toBe("https://example.com/file.ion");
      expect(calls[0][1]).toBe(0);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });

    test("parseIon accepts an http URL string and routes it to parseIonRemote", async () => {
      await parseIon("http://example.com/file.ion");

      const calls = (mockBackend.parseIonRemote as jest.Mock).mock.calls;
      expect(calls.length).toBe(1);
      expect(calls[0][0].href).toBe("http://example.com/file.ion");
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
    });

    test("parseIon with URL and cache size calls parseIonRemote with cache size", async () => {
      const url = new URL("https://example.com/file.ion");
      const cacheSize = 2000;
      await parseIon(url, { maxCacheSize: cacheSize });

      expect(mockBackend.parseIonRemote).toHaveBeenCalledWith(url, cacheSize);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonBuffer).not.toHaveBeenCalled();
    });
  });

  describe("buffer dispatch", () => {
    test("parseIon with Uint8Array calls parseIonBuffer only", async () => {
      const bytes = new Uint8Array([1, 2, 3]);
      await parseIon(bytes);

      const callArg = (mockBackend.parseIonBuffer as jest.Mock).mock.calls[0];
      expect(callArg[0]).toEqual(bytes);
      expect(callArg[1]).toBe(0);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
    });

    test("parseIon with ArrayBuffer calls parseIonBuffer (normalized to Uint8Array)", async () => {
      const arrayBuffer = new ArrayBuffer(4);
      await parseIon(arrayBuffer);

      const callArg = (mockBackend.parseIonBuffer as jest.Mock).mock.calls[0];
      expect(callArg[0]).toBeInstanceOf(Uint8Array);
      expect(callArg[0].byteLength).toBe(4);
      expect(callArg[1]).toBe(0);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
    });

    test("parseIon with buffer and cache size calls parseIonBuffer with cache size", async () => {
      const bytes = new Uint8Array([5, 6, 7]);
      const cacheSize = 3000;
      await parseIon(bytes, { maxCacheSize: cacheSize });

      const callArg = (mockBackend.parseIonBuffer as jest.Mock).mock.calls[0];
      expect(callArg[0]).toEqual(bytes);
      expect(callArg[1]).toBe(cacheSize);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
    });

    test("parseIon with ArrayBuffer and cache size normalizes and passes cache size", async () => {
      const arrayBuffer = new ArrayBuffer(8);
      const cacheSize = 5000;
      await parseIon(arrayBuffer, { maxCacheSize: cacheSize });

      const callArg = (mockBackend.parseIonBuffer as jest.Mock).mock.calls[0];
      expect(callArg[0]).toBeInstanceOf(Uint8Array);
      expect(callArg[0].byteLength).toBe(8);
      expect(callArg[1]).toBe(cacheSize);
      expect(mockBackend.parseIonPath).not.toHaveBeenCalled();
      expect(mockBackend.parseIonRemote).not.toHaveBeenCalled();
    });
  });

  describe("dispatch isolation", () => {
    test("each input kind calls exactly one primitive", async () => {
      const path = "/data/file.ion";
      const url = new URL("https://example.com/file.ion");
      const bytes = new Uint8Array([1, 2, 3]);

      await parseIon(path);
      expect(
        (mockBackend.parseIonPath as jest.Mock).mock.calls.length,
      ).toBe(1);
      expect(
        (mockBackend.parseIonRemote as jest.Mock).mock.calls.length,
      ).toBe(0);
      expect(
        (mockBackend.parseIonBuffer as jest.Mock).mock.calls.length,
      ).toBe(0);

      (mockBackend.parseIonPath as jest.Mock).mockClear();
      (mockBackend.parseIonRemote as jest.Mock).mockClear();
      (mockBackend.parseIonBuffer as jest.Mock).mockClear();

      await parseIon(url);
      expect(
        (mockBackend.parseIonPath as jest.Mock).mock.calls.length,
      ).toBe(0);
      expect((mockBackend.parseIonRemote as jest.Mock).mock.calls.length).toBe(1);
      expect(
        (mockBackend.parseIonBuffer as jest.Mock).mock.calls.length,
      ).toBe(0);

      (mockBackend.parseIonPath as jest.Mock).mockClear();
      (mockBackend.parseIonRemote as jest.Mock).mockClear();
      (mockBackend.parseIonBuffer as jest.Mock).mockClear();

      await parseIon(bytes);
      expect(
        (mockBackend.parseIonPath as jest.Mock).mock.calls.length,
      ).toBe(0);
      expect(
        (mockBackend.parseIonRemote as jest.Mock).mock.calls.length,
      ).toBe(0);
      expect(
        (mockBackend.parseIonBuffer as jest.Mock).mock.calls.length,
      ).toBe(1);
    });
  });
});
