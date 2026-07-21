import { describe, test, expect, beforeEach } from "@jest/globals";

const modulePaths: Array<[string, string, string]> = [
  ["source", "../utilities/api", "../utilities/sampleFile"],
  ["built", "../../lib/utilities/api.js", "../../lib/utilities/sampleFile.js"],
];

describe.each(modulePaths)(
  "getPeaksFromChrom index validation (%s)",
  (_name, apiPath, sampleFilePath) => {
    const api = require(apiPath);
    const { SampleFile } = require(sampleFilePath);

    let backend: any;
    let file: any;

    beforeEach(() => {
      backend = {
        ready: true,
        handlesAreGcFinalized: true,
        getPeaksFromChrom: jest.fn(() => []),
        freeFile: jest.fn(),
      };
      api.setBackend(backend);
      file = new SampleFile({}, backend);
    });

    test("a negative index throws instead of returning rows", () => {
      expect(() =>
        api.getPeaksFromChrom(file, [{ id: "a", idx: -1, rt: 1, range: 0.1 }], {}, 1),
      ).toThrow(RangeError);
      expect(backend.getPeaksFromChrom).not.toHaveBeenCalled();
    });

    test("the error names the item and what it needs", () => {
      expect(() =>
        api.getPeaksFromChrom(
          file,
          [
            { id: "a", idx: 0, rt: 1, range: 0.1 },
            { id: "b", idx: -1, rt: 2, range: 0.1 },
          ],
          {},
          1,
        ),
      ).toThrow(/item 1 needs a non-negative integer idx/);
    });

    test("a missing index throws", () => {
      expect(() =>
        api.getPeaksFromChrom(file, [{ id: "a", rt: 1, range: 0.1 }], {}, 1),
      ).toThrow(/non-negative integer idx/);
      expect(backend.getPeaksFromChrom).not.toHaveBeenCalled();
    });

    test("a fractional index throws", () => {
      expect(() =>
        api.getPeaksFromChrom(file, [{ id: "a", idx: 1.5, rt: 1, range: 0.1 }], {}, 1),
      ).toThrow(/non-negative integer idx/);
    });

    test("a NaN index throws", () => {
      expect(() =>
        api.getPeaksFromChrom(file, [{ id: "a", idx: NaN, rt: 1, range: 0.1 }], {}, 1),
      ).toThrow(/non-negative integer idx/);
    });

    test("an infinite index throws", () => {
      expect(() =>
        api.getPeaksFromChrom(
          file,
          [{ id: "a", idx: Infinity, rt: 1, range: 0.1 }],
          {},
          1,
        ),
      ).toThrow(/non-negative integer idx/);
    });

    test("valid indices reach the backend unchanged", () => {
      api.getPeaksFromChrom(
        file,
        [
          { id: "a", idx: 0, rt: 1, range: 0.1 },
          { id: "b", idx: 3, rt: 2, window: 0.2 },
        ],
        {},
        1,
      );

      const call = backend.getPeaksFromChrom.mock.calls[0];
      expect(Array.from(call[1])).toEqual([0, 3]);
      expect(Array.from(call[2])).toEqual([1, 2]);
      expect(Array.from(call[3])).toEqual([0.1, 0.2]);
      expect(call[4]).toBe(2);
    });

    test("the index alias is accepted", () => {
      api.getPeaksFromChrom(file, [{ id: "a", index: 2, rt: 1, range: 0.1 }], {}, 1);

      const call = backend.getPeaksFromChrom.mock.calls[0];
      expect(Array.from(call[1])).toEqual([2]);
    });

    test("an empty item list reaches the backend", () => {
      api.getPeaksFromChrom(file, [], {}, 1);
      expect(backend.getPeaksFromChrom).toHaveBeenCalledTimes(1);
    });
  },
);
