import { describe, test, expect, beforeAll, afterAll } from "@jest/globals";
import * as http from "http";
import type { AddressInfo } from "net";

type RemoteSourceModule = {
  newRemoteSource: (url: string) => any;
  fetchHeader: (source: any) => Promise<Uint8Array>;
};

const modulePaths: Array<[string, string]> = [
  ["source", "../utilities/remoteSource"],
  ["built", "../../lib/utilities/remoteSource.js"],
];

const FULL_TOTAL = 34353;
const TINY_LENGTH = 812;

function body(length: number): Buffer {
  return Buffer.alloc(length, 7);
}

describe.each(modulePaths)("fetchHeader (%s)", (_name, modulePath) => {
  const remoteSource: RemoteSourceModule = require(modulePath);

  let server: http.Server;
  let origin: string;

  beforeAll((done) => {
    server = http.createServer((request, response) => {
      if (request.url === "/full") {
        const bytes = body(1024);
        response.writeHead(206, {
          "content-length": bytes.length,
          "content-range": `bytes 0-1023/${FULL_TOTAL}`,
        });
        response.end(bytes);
        return;
      }
      if (request.url === "/tiny") {
        const bytes = body(TINY_LENGTH);
        response.writeHead(206, {
          "content-length": bytes.length,
          "content-range": `bytes 0-${TINY_LENGTH - 1}/${TINY_LENGTH}`,
        });
        response.end(bytes);
        return;
      }
      if (request.url === "/star") {
        const bytes = body(1024);
        response.writeHead(206, {
          "content-length": bytes.length,
          "content-range": "bytes 0-1023/*",
        });
        response.end(bytes);
        return;
      }
      if (request.url === "/nototal") {
        const bytes = body(1024);
        response.writeHead(206, { "content-length": bytes.length });
        response.end(bytes);
        return;
      }
      if (request.url === "/plain") {
        const bytes = body(1024);
        response.writeHead(200, { "content-length": bytes.length });
        response.end(bytes);
        return;
      }
      response.writeHead(404);
      response.end("Not Found");
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

  async function failureFor(route: string): Promise<Error> {
    const source = remoteSource.newRemoteSource(`${origin}${route}`);
    try {
      await remoteSource.fetchHeader(source);
    } catch (error) {
      return error as Error;
    }
    throw new Error(`fetchHeader should have failed for ${route}`);
  }

  test("a valid range response records the total size", async () => {
    const source = remoteSource.newRemoteSource(`${origin}/full`);
    const bytes = await remoteSource.fetchHeader(source);

    expect(bytes.length).toBe(1024);
    expect(source.total).toBe(BigInt(FULL_TOTAL));
    expect(source.cache.read(0n, 1024n)).not.toBeNull();
  });

  test("a file shorter than the header window opens", async () => {
    const source = remoteSource.newRemoteSource(`${origin}/tiny`);
    const bytes = await remoteSource.fetchHeader(source);

    expect(bytes.length).toBe(TINY_LENGTH);
    expect(source.total).toBe(BigInt(TINY_LENGTH));
  });

  test("the cached range matches the received length, not the requested one", async () => {
    const source = remoteSource.newRemoteSource(`${origin}/tiny`);
    await remoteSource.fetchHeader(source);

    const cached = source.cache.read(0n, BigInt(TINY_LENGTH));
    expect(cached).not.toBeNull();
    expect(cached.length).toBe(TINY_LENGTH);
    expect(source.cache.read(0n, 1024n)).toBeNull();
  });

  test("an unknown total size is a named error, not a SyntaxError", async () => {
    const error = await failureFor("/star");

    expect(error).not.toBeInstanceOf(SyntaxError);
    expect(error.message).toContain(`${origin}/star`);
    expect(error.message).toContain("bytes 0-1023/*");
  });

  test("a missing Content-Range is a named error", async () => {
    const error = await failureFor("/nototal");

    expect(error).not.toBeInstanceOf(SyntaxError);
    expect(error.message).toContain(`${origin}/nototal`);
    expect(error.message).toContain("no total size");
  });

  test("a 200 response names the URL", async () => {
    const error = await failureFor("/plain");

    expect(error.message).toContain(`${origin}/plain`);
    expect(error.message).toContain("200");
  });

  test("a 404 response names the URL", async () => {
    const error = await failureFor("/missing");

    expect(error.message).toContain(`${origin}/missing`);
    expect(error.message).toContain("404");
  });
});
