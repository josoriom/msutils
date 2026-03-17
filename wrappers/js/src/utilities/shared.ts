import type { BinaryInput } from "../types/types";

export function snakeToCamel(s: string): string {
  return s.replace(/_([a-z])/g, (_, c) => c.toUpperCase());
}

export function camelizeKeys<T>(value: T): T {
  if (Array.isArray(value)) return value.map(camelizeKeys) as unknown as T;
  if (value && typeof value === "object" && !ArrayBuffer.isView(value)) {
    const result: any = {};
    for (const [key, val] of Object.entries(value))
      result[snakeToCamel(key)] = camelizeKeys(val);
    return result as T;
  }
  return value;
}

export function toCores(cores: unknown): number {
  const v = typeof cores === "number" && Number.isFinite(cores) ? cores | 0 : 1;
  return v > 0 ? v : 1;
}

export function toUint8(input: BinaryInput): Uint8Array {
  if (input instanceof Uint8Array) return input;
  return new Uint8Array(input);
}
