import { describe, test, expect } from "@jest/globals";
import * as fs from "fs";
import * as path from "path";

const removedNames = ["mzmlToIonFile", "fitPeak", "drawPeak"];

const modulePaths: Array<[string, string]> = [
  ["source", "../utilities/api"],
  ["built", "../../lib/utilities/api.js"],
];

describe.each(modulePaths)("removed api stubs (%s)", (_name, modulePath) => {
  const api = require(modulePath);

  test.each(removedNames)("%s is not exported", (name) => {
    expect(name in api).toBe(false);
    expect(api[name]).toBeUndefined();
  });

  test("the surviving conversion and peak calls are still exported", () => {
    expect(typeof api.parseIon).toBe("function");
    expect(typeof api.ionToMzml).toBe("function");
    expect(typeof api.getPeak).toBe("function");
  });
});

describe("removed api stubs (shipped types)", () => {
  const declarations = fs.readFileSync(
    path.join(__dirname, "..", "..", "lib", "utilities", "api.d.ts"),
    "utf8",
  );

  test.each(removedNames)("%s is not declared", (name) => {
    expect(declarations).not.toContain(`declare function ${name}`);
    expect(declarations).not.toContain(`${name}(`);
  });

  test("the shipped types still declare parseIon", () => {
    expect(declarations).toContain("declare function parseIon");
  });
});
