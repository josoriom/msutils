import {
  describe,
  test,
  expect,
  beforeAll,
  afterAll,
  beforeEach,
} from "@jest/globals";
import * as http from "http";
import type { AddressInfo } from "net";

type RemoteSourceModule = {
  newRemoteSource: (url: string) => any;
  fetchHeader: (source: any) => Promise<Uint8Array>;
  fetchRange: (
    url: string,
    range: { offset: bigint; length: bigint },
  ) => Promise<Uint8Array>;
};

const modulePaths: Array<[string, string]> = [
  ["source", "../utilities/remoteSource"],
  ["built", "../../lib/utilities/remoteSource.js"],
];

const FULL_TOTAL = 34353;
const TINY_LENGTH = 812;
const HIDDEN_TOTAL = 56892162;
const BLOCK_OFFSET = 2048n;
const BLOCK_LENGTH = 256;

const totalByRoute: Record<string, number> = {
  "/nototal": HIDDEN_TOTAL,
  "/star": HIDDEN_TOTAL,
};

function body(length: number): Buffer {
  return Buffer.alloc(length, 7);
}

describe.each(modulePaths)("fetchHeader (%s)", (_name, modulePath) => {
  const remoteSource: RemoteSourceModule = require(modulePath);

  let server: http.Server;
  let origin: string;
  let headRoutes: string[] = [];

  beforeAll((done) => {
    server = http.createServer((request, response) => {
      if (request.method === "HEAD") {
        const route = request.url ?? "";
        headRoutes.push(route);

        if (route === "/headfails") {
          response.writeHead(500);
          response.end();
          return;
        }
        if (route === "/gzipped") {
          response.writeHead(200, {
            "content-length": HIDDEN_TOTAL,
            "content-encoding": "gzip",
          });
          response.end();
          return;
        }
        if (route === "/nolength") {
          response.writeHead(200, { "transfer-encoding": "chunked" });
          response.end();
          return;
        }

        const total = totalByRoute[route];
        if (total === undefined) {
          response.writeHead(404);
          response.end();
          return;
        }
        response.writeHead(200, { "content-length": total });
        response.end();
        return;
      }

      if (
        request.url === "/headfails" ||
        request.url === "/gzipped" ||
        request.url === "/nolength"
      ) {
        const bytes = body(1024);
        response.writeHead(206, { "content-length": bytes.length });
        response.end(bytes);
        return;
      }
      if (request.url === "/block-hidden") {
        const bytes = body(BLOCK_LENGTH);
        response.writeHead(206, { "content-length": bytes.length });
        response.end(bytes);
        return;
      }
      if (request.url === "/block-wrong-offset") {
        const bytes = body(BLOCK_LENGTH);
        response.writeHead(206, {
          "content-length": bytes.length,
          "content-range": `bytes 4096-${4096 + BLOCK_LENGTH - 1}/${HIDDEN_TOTAL}`,
        });
        response.end(bytes);
        return;
      }
      if (request.url === "/block-short") {
        const bytes = body(BLOCK_LENGTH - 1);
        response.writeHead(206, { "content-length": bytes.length });
        response.end(bytes);
        return;
      }
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

  beforeEach(() => {
    headRoutes = [];
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

  test("a visible Content-Range needs no size request", async () => {
    const source = remoteSource.newRemoteSource(`${origin}/full`);
    await remoteSource.fetchHeader(source);

    expect(headRoutes).toEqual([]);
  });

  test("a hidden Content-Range takes the total from a size request", async () => {
    const source = remoteSource.newRemoteSource(`${origin}/nototal`);
    const bytes = await remoteSource.fetchHeader(source);

    expect(bytes.length).toBe(1024);
    expect(source.total).toBe(BigInt(HIDDEN_TOTAL));
    expect(headRoutes).toEqual(["/nototal"]);
  });

  test("an unknown total size falls back to a size request", async () => {
    const source = remoteSource.newRemoteSource(`${origin}/star`);
    await remoteSource.fetchHeader(source);

    expect(source.total).toBe(BigInt(HIDDEN_TOTAL));
    expect(headRoutes).toEqual(["/star"]);
  });

  test("a failed size request names the URL", async () => {
    const error = await failureFor("/headfails");

    expect(error).not.toBeInstanceOf(SyntaxError);
    expect(error.message).toContain(`${origin}/headfails`);
    expect(error.message).toContain("500");
  });

  test("a compressed size request is refused, not trusted", async () => {
    const error = await failureFor("/gzipped");

    expect(error.message).toContain(`${origin}/gzipped`);
    expect(error.message).toContain("gzip");
  });

  test("a size request without Content-Length names the URL", async () => {
    const error = await failureFor("/nolength");

    expect(error).not.toBeInstanceOf(SyntaxError);
    expect(error.message).toContain(`${origin}/nolength`);
    expect(error.message).toContain("Content-Length");
  });

  async function blockFailureFor(route: string): Promise<Error> {
    try {
      await remoteSource.fetchRange(`${origin}${route}`, {
        offset: BLOCK_OFFSET,
        length: BigInt(BLOCK_LENGTH),
      });
    } catch (error) {
      return error as Error;
    }
    throw new Error(`fetchRange should have failed for ${route}`);
  }

  test("a block with a hidden Content-Range is accepted", async () => {
    const bytes = await remoteSource.fetchRange(`${origin}/block-hidden`, {
      offset: BLOCK_OFFSET,
      length: BigInt(BLOCK_LENGTH),
    });

    expect(bytes.length).toBe(BLOCK_LENGTH);
  });

  test("a block sent from the wrong offset is still refused", async () => {
    const error = await blockFailureFor("/block-wrong-offset");

    expect(error.message).toContain("wrong Content-Range");
  });

  test("a block with the wrong length is still refused", async () => {
    const error = await blockFailureFor("/block-short");

    expect(error.message).toContain("length mismatch");
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
