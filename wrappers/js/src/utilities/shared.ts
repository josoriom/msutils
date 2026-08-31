import type { BinaryInput } from "../types/types";

export function toCores(cores: unknown): number {
  const v = typeof cores === "number" && Number.isFinite(cores) ? cores | 0 : 1;
  return v > 0 ? v : 1;
}

export function toUint8(input: BinaryInput): Uint8Array {
  if (input instanceof Uint8Array) return input;
  return new Uint8Array(input);
}
