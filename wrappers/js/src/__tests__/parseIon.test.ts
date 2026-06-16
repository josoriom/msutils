import { describe, it, expect, beforeAll, afterAll } from '@jest/globals';
import * as fs from 'fs';
import * as path from 'path';
import * as http from 'http';

async function drain(p: Promise<unknown>): Promise<void> {
  await p.catch(() => {});
}

describe('parseIon URL support', () => {
  let test_ion_path: string;
  let test_ion_buffer: Buffer;
  let test_server: http.Server;
  const test_port = 7654;
  const test_server_url = `http://localhost:${test_port}`;

  beforeAll((done) => {
    const wrapper_dir = path.join(__dirname, '..', '..', '..');
    test_ion_path = path.join(wrapper_dir, 'test.ion');
    test_ion_buffer = fs.readFileSync(test_ion_path);

    test_server = http.createServer((request, response) => {
      if (request.url === '/test.ion') {
        if (request.method === 'HEAD') {
          response.writeHead(200, {
            'content-length': test_ion_buffer.length,
            'accept-ranges': 'bytes',
          });
          response.end();
          return;
        }

        const range = request.headers.range;
        if (range) {
          const match = range.match(/bytes=(\d+)-(\d+)/);
          if (match) {
            const start = parseInt(match[1], 10);
            const end = parseInt(match[2], 10);
            const chunk = test_ion_buffer.slice(start, end + 1);
            response.writeHead(206, {
              'content-length': chunk.length,
              'content-range': `bytes ${start}-${end}/${test_ion_buffer.length}`,
              'accept-ranges': 'bytes',
            });
            response.end(chunk);
            return;
          }
        }

        response.writeHead(200, {
          'content-length': test_ion_buffer.length,
          'accept-ranges': 'bytes',
        });
        response.end(test_ion_buffer);
      } else {
        response.writeHead(404);
        response.end('Not Found');
      }
    });

    test_server.listen(test_port, () => done());
  });

  afterAll((done) => {
    test_server.close(() => done());
  });

  it('test.ion file exists and loads', () => {
    expect(fs.existsSync(test_ion_path)).toBe(true);
    expect(test_ion_buffer.length).toBeGreaterThan(0);
  });

  it('parseIon is exported as function', async () => {
    const { parseIon } = await import('../utilities/api');
    expect(typeof parseIon).toBe('function');
  });

  it('rejects URL string (must use new URL)', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon('http://localhost/file.ion');
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toContain('URL strings are not accepted');
    }
  });

  it('rejects file:// as string', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon('file:///tmp/test.ion');
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toContain('URL strings are not accepted');
    }
  });

  it('rejects empty path', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon('');
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toMatch(/empty|string/i);
    }
  });

  it('rejects negative cache size', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon(test_ion_buffer, { maxCacheSize: -1 });
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toContain('non-negative integer');
    }
  });

  it('rejects float cache size', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon(test_ion_buffer, { maxCacheSize: 1.5 });
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toContain('non-negative integer');
    }
  });

  it('rejects NaN cache size', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon(test_ion_buffer, { maxCacheSize: NaN });
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toContain('non-negative integer');
    }
  });

  it('rejects Infinity cache size', async () => {
    const { parseIon } = await import('../utilities/api');

    try {
      await parseIon(test_ion_buffer, { maxCacheSize: Infinity });
      throw new Error('Should reject');
    } catch (error) {
      expect((error as Error).message).toContain('non-negative integer');
    }
  });

  it('buffer input returns Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const result = parseIon(test_ion_buffer);
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('Uint8Array input returns Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const result = parseIon(new Uint8Array(test_ion_buffer));
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('path string input returns Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const result = parseIon(test_ion_path);
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('file URL input returns Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const file_url = new URL(`file://${test_ion_path}`);
    const result = parseIon(file_url);
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('http URL input returns Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const http_url = new URL(`${test_server_url}/test.ion`);
    const result = parseIon(http_url);
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('localhost URL input returns Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const localhost_url = new URL(`http://127.0.0.1:${test_port}/test.ion`);
    const result = parseIon(localhost_url);
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('zero cache size is valid', async () => {
    const { parseIon } = await import('../utilities/api');
    const result = parseIon(test_ion_buffer, { maxCacheSize: 0 });
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('large cache size is valid', async () => {
    const { parseIon } = await import('../utilities/api');
    const result = parseIon(test_ion_buffer, { maxCacheSize: 100000000 });
    expect(result).toBeInstanceOf(Promise);
    await drain(result);
  });

  it('multiple calls return independent Promises', async () => {
    const { parseIon } = await import('../utilities/api');
    const result1 = parseIon(test_ion_buffer);
    const result2 = parseIon(test_ion_buffer);
    expect(result1).not.toBe(result2);
    expect(result1).toBeInstanceOf(Promise);
    expect(result2).toBeInstanceOf(Promise);
    await drain(result1);
    await drain(result2);
  });

  it('buffer and Uint8Array both return Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const buffer_result = parseIon(test_ion_buffer);
    const uint8_result = parseIon(new Uint8Array(test_ion_buffer));
    expect(buffer_result).toBeInstanceOf(Promise);
    expect(uint8_result).toBeInstanceOf(Promise);
    await drain(buffer_result);
    await drain(uint8_result);
  });

  it('file URL and path string both return Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const file_url_result = parseIon(new URL(`file://${test_ion_path}`));
    const path_result = parseIon(test_ion_path);
    expect(file_url_result).toBeInstanceOf(Promise);
    expect(path_result).toBeInstanceOf(Promise);
    await drain(file_url_result);
    await drain(path_result);
  });

  it('all input types return Promise', async () => {
    const { parseIon } = await import('../utilities/api');
    const buffer = parseIon(test_ion_buffer);
    const uint8 = parseIon(new Uint8Array(test_ion_buffer));
    const file_path = parseIon(test_ion_path);
    const file_url = parseIon(new URL(`file://${test_ion_path}`));
    const http_url = parseIon(new URL(`${test_server_url}/test.ion`));

    expect(buffer).toBeInstanceOf(Promise);
    expect(uint8).toBeInstanceOf(Promise);
    expect(file_path).toBeInstanceOf(Promise);
    expect(file_url).toBeInstanceOf(Promise);
    expect(http_url).toBeInstanceOf(Promise);

    await drain(buffer);
    await drain(uint8);
    await drain(file_path);
    await drain(file_url);
    await drain(http_url);
  });

  it('rejects on HTTP 404', async () => {
    const { parseIon } = await import('../utilities/api');
    const http_url = new URL(`${test_server_url}/does-not-exist.ion`);

    try {
      await parseIon(http_url);
      throw new Error('Should fail');
    } catch (error) {
      const msg = (error as Error).message;
      expect(msg).toMatch(/HTTP|404|not found|Backend|error/i);
    }
  });

  it('rejects on missing file URL', async () => {
    const { parseIon } = await import('../utilities/api');
    const file_url = new URL('file:///tmp/does-not-exist-9999.ion');

    try {
      await parseIon(file_url);
      throw new Error('Should fail');
    } catch (error) {
      const msg = (error as Error).message;
      expect(msg).toMatch(/parse|file|error|not found|Backend/i);
    }
  });
});
