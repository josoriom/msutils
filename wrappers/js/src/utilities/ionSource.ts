import type { IonInput } from "../types/types";
import { toUint8 } from "./shared";

export type IonSource =
  | { kind: "path"; path: string }
  | { kind: "url"; url: URL }
  | { kind: "buffer"; bytes: Uint8Array };

function isUrlString(value: string): boolean {
  return (
    /^[A-Za-z][A-Za-z0-9+.-]*:/.test(value) && !/^[A-Za-z]:[\\/]/.test(value)
  );
}

function getPathSource(value: string): IonSource {
  if (value.length === 0) {
    throw new TypeError("parseIon: path string must not be empty");
  }
  if (/^https?:/i.test(value)) {
    return { kind: "url", url: new URL(value) };
  }
  if (isUrlString(value)) {
    throw new TypeError(
      "parseIon: only http and https URLs can be read remotely",
    );
  }
  return { kind: "path", path: value };
}

function isBytes(value: unknown): value is Uint8Array | ArrayBuffer {
  return value instanceof Uint8Array || value instanceof ArrayBuffer;
}

export function getIonSource(input: IonInput): IonSource {
  if (typeof input === "string") return getPathSource(input);
  if (input instanceof URL) {
    if (input.protocol !== "http:" && input.protocol !== "https:") {
      throw new TypeError(
        "parseIon: only http and https URLs can be read remotely",
      );
    }
    return { kind: "url", url: input };
  }
  if (isBytes(input)) return { kind: "buffer", bytes: toUint8(input) };
  throw new TypeError(
    "parseIon: source must be a path string, a URL, or a byte buffer",
  );
}
