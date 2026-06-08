import { describe, test, expect } from "@jest/globals";
import { getIonSource, type IonSource } from "../utilities/ionSource";

describe("ionSource classifier", () => {
  describe("path source", () => {
    test("empty string throws", () => {
      expect(() => getIonSource("")).toThrow(
        "path string must not be empty",
      );
    });

    test("string with http scheme throws", () => {
      expect(() => getIonSource("http://example.com/file.ion")).toThrow(
        "URL strings are not accepted",
      );
    });

    test("string with https scheme throws", () => {
      expect(() => getIonSource("https://example.com/file.ion")).toThrow(
        "URL strings are not accepted",
      );
    });

    test("string with file scheme throws", () => {
      expect(() => getIonSource("file:///data/file.ion")).toThrow(
        "URL strings are not accepted",
      );
    });

    test("absolute path returns path source", () => {
      const source = getIonSource("/data/file.ion");
      expect(source).toEqual({ kind: "path", path: "/data/file.ion" });
    });

    test("relative path returns path source", () => {
      const source = getIonSource("./data/file.ion");
      expect(source).toEqual({ kind: "path", path: "./data/file.ion" });
    });

    test("windows absolute path returns path source", () => {
      const source = getIonSource("C:\\data\\file.ion");
      expect(source).toEqual({ kind: "path", path: "C:\\data\\file.ion" });
    });
  });

  describe("url source", () => {
    test("http URL returns url source", () => {
      const url = new URL("http://example.com/file.ion");
      const source = getIonSource(url);
      expect(source.kind).toBe("url");
      expect((source as { kind: "url"; url: URL }).url).toBe(url);
    });

    test("https URL returns url source", () => {
      const url = new URL("https://example.com/file.ion");
      const source = getIonSource(url);
      expect(source.kind).toBe("url");
      expect((source as { kind: "url"; url: URL }).url).toBe(url);
    });

    test("file URL returns url source", () => {
      const url = new URL("file:///data/file.ion");
      const source = getIonSource(url);
      expect(source.kind).toBe("url");
      expect((source as { kind: "url"; url: URL }).url).toBe(url);
    });
  });

  describe("buffer source", () => {
    test("Uint8Array returns buffer source", () => {
      const bytes = new Uint8Array([1, 2, 3]);
      const source = getIonSource(bytes);
      expect(source.kind).toBe("buffer");
      expect((source as { kind: "buffer"; bytes: Uint8Array }).bytes).toEqual(
        bytes,
      );
    });

    test("ArrayBuffer returns buffer source", () => {
      const buffer = new ArrayBuffer(4);
      const source = getIonSource(buffer);
      expect(source.kind).toBe("buffer");
      const bytes = (source as { kind: "buffer"; bytes: Uint8Array }).bytes;
      expect(bytes).toBeInstanceOf(Uint8Array);
      expect(bytes.byteLength).toBe(4);
    });
  });

  describe("invalid input", () => {
    test("null throws", () => {
      expect(() => getIonSource(null as any)).toThrow(
        "must be a path string, a URL, or a byte buffer",
      );
    });

    test("undefined throws", () => {
      expect(() => getIonSource(undefined as any)).toThrow(
        "must be a path string, a URL, or a byte buffer",
      );
    });

    test("number throws", () => {
      expect(() => getIonSource(42 as any)).toThrow(
        "must be a path string, a URL, or a byte buffer",
      );
    });

    test("plain object throws", () => {
      expect(() => getIonSource({} as any)).toThrow(
        "must be a path string, a URL, or a byte buffer",
      );
    });
  });
});
