import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";
import path from "node:path";
import fs from "node:fs";

const require = createRequire(import.meta.url);
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../..");
const quantion = require(path.join(root, "wrappers/js/lib/index-node.js"));

const MASS = 90.0555;
const FROM_RT = 0.0;
const TO_RT = 5.0;
const EIC_PPM = 20.0;
const EIC_MZ = 0.005;

function sameBits(left, right) {
  if (left.length !== right.length) return false;
  const view = new DataView(new ArrayBuffer(16));
  for (let index = 0; index < left.length; index += 1) {
    view.setFloat64(0, left[index]);
    view.setFloat64(8, right[index]);
    if (view.getBigUint64(0) !== view.getBigUint64(8)) return false;
  }
  return true;
}

async function stats(url) {
  const base = new URL(url);
  const response = await fetch(`${base.origin}/stats`);
  return response.json();
}

async function reset(url) {
  const base = new URL(url);
  await fetch(`${base.origin}/reset`);
}

async function main() {
  const url = process.env.QUANTION_REMOTE_URL;
  const fixture = path.join(root, "core/tests/fixtures/api/api.ion");
  if (!url) {
    process.stderr.write("QUANTION_REMOTE_URL is not set\n");
    return 1;
  }

  await reset(url);
  const remote = await quantion.parseIon(new URL(url));
  const opening = await stats(url);

  await reset(url);
  const local = await quantion.parseIon(fixture);

  const eicRemote = await quantion.calculateEic(
    remote,
    MASS,
    { from: FROM_RT, to: TO_RT },
    EIC_PPM,
    EIC_MZ,
  );
  const eicLocal = await quantion.calculateEic(
    local,
    MASS,
    { from: FROM_RT, to: TO_RT },
    EIC_PPM,
    EIC_MZ,
  );
  const oneQuery = await stats(url);

  const total = fs.statSync(fixture).size;
  const same =
    sameBits(eicRemote.x, eicLocal.x) && sameBits(eicRemote.y, eicLocal.y);

  console.log(`file_bytes = ${total}`);
  console.log(`opening_bytes = ${opening.bytes_sent}`);
  console.log(`opening_requests = ${opening.requests}`);
  console.log(`query_bytes = ${oneQuery.bytes_sent}`);
  console.log(`query_requests = ${oneQuery.requests}`);
  console.log(`matches_local = ${same ? "yes" : "no"}`);
  return 0;
}

main().then((code) => process.exit(code)).catch((error) => {
  process.stderr.write(String(error?.stack ?? error) + "\n");
  process.exit(1);
});
