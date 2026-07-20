import { test, expect } from "@jest/globals";
import { execFileSync } from "child_process";
import * as path from "path";

const packageRoot = path.resolve(__dirname, "../..");

function loadWithRoot(root: string): string {
  try {
    execFileSync(process.execPath, ["-e", "require('./lib/index-node.js')"], {
      cwd: packageRoot,
      env: { ...process.env, QUANTION_ARTIFACTS_ROOT: root, QUANTION_LIB: "" },
      encoding: "utf8",
      stdio: "pipe",
    });
    return "";
  } catch (error) {
    return String((error as { stderr?: string }).stderr ?? "");
  }
}

test("a missing artifacts root gives the clear message, not a raw scandir error", () => {
  const message = loadWithRoot("/definitely/not/here");
  expect(message).not.toMatch(/ENOENT/);
  expect(message).toMatch(/QUANTION_ARTIFACTS_ROOT/);
});
