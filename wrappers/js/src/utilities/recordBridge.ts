import {
  PAYLOAD_CHROM_PEAKS,
  PAYLOAD_CONSENSUS_FEATURES,
  PAYLOAD_EIC_PEAKS,
  PAYLOAD_FEATURES,
  PAYLOAD_FIT_RESULT,
  PAYLOAD_FOUND_FEATURES,
  PAYLOAD_PEAKS,
  SECTION_CHROM_FROM,
  SECTION_CHROM_ID_BYTES,
  SECTION_CHROM_ID_STARTS,
  SECTION_CHROM_INDEX,
  SECTION_CHROM_INTEGRAL,
  SECTION_CHROM_INTENSITY,
  SECTION_CHROM_RT,
  SECTION_CHROM_TARGET_RT,
  SECTION_CHROM_TIMESTAMP_BYTES,
  SECTION_CHROM_TIMESTAMP_STARTS,
  SECTION_CHROM_TO,
  SECTION_CHROM_TOTAL_AREA,
  SECTION_CONSENSUS_FREQUENCY,
  SECTION_CONSENSUS_FROM,
  SECTION_CONSENSUS_INTEGRAL,
  SECTION_CONSENSUS_INTENSITY,
  SECTION_CONSENSUS_MZ,
  SECTION_CONSENSUS_RT,
  SECTION_CONSENSUS_TO,
  SECTION_EIC_PEAK_FROM,
  SECTION_EIC_PEAK_ID_BYTES,
  SECTION_EIC_PEAK_ID_STARTS,
  SECTION_EIC_PEAK_INTEGRAL,
  SECTION_EIC_PEAK_INTENSITY,
  SECTION_EIC_PEAK_MZ,
  SECTION_EIC_PEAK_NOISE,
  SECTION_EIC_PEAK_ORT,
  SECTION_EIC_PEAK_RT,
  SECTION_EIC_PEAK_TO,
  SECTION_FEATURE_FROM,
  SECTION_FEATURE_INTEGRAL,
  SECTION_FEATURE_INTENSITY,
  SECTION_FEATURE_MZ,
  SECTION_FEATURE_NOISE,
  SECTION_FEATURE_POINT_COUNT,
  SECTION_FEATURE_RT,
  SECTION_FEATURE_TO,
  SECTION_FIT_CENTER,
  SECTION_FIT_FWHM,
  SECTION_FIT_HEIGHT,
  SECTION_FIT_R2,
  SECTION_FIT_SHAPE,
  SECTION_FIT_TAIL,
  SECTION_FOUND_FROM,
  SECTION_FOUND_ID_BYTES,
  SECTION_FOUND_ID_STARTS,
  SECTION_FOUND_INTEGRAL,
  SECTION_FOUND_INTENSITY,
  SECTION_FOUND_MZ,
  SECTION_FOUND_NOISE,
  SECTION_FOUND_POINT_COUNT,
  SECTION_FOUND_RT,
  SECTION_FOUND_TO,
  SECTION_PEAK_FROM,
  SECTION_PEAK_INTEGRAL,
  SECTION_PEAK_INTENSITY,
  SECTION_PEAK_NOISE,
  SECTION_PEAK_POINT_COUNT,
  SECTION_PEAK_R2,
  SECTION_PEAK_RT,
  SECTION_PEAK_TO,
  readCountSection,
  readF64Section,
  readBridge,
  readU32Section,
  readU8Section,
} from "./bridge";

type Column =
  | { name: string; id: number; kind: "number" }
  | { name: string; id: number; kind: "count" }
  | { name: string; startsId: number; bytesId: number; kind: "text" };

const RECORD_COLUMNS: Record<number, Column[]> = {
  [PAYLOAD_PEAKS]: [
    { name: "from", id: SECTION_PEAK_FROM, kind: "number" },
    { name: "to", id: SECTION_PEAK_TO, kind: "number" },
    { name: "rt", id: SECTION_PEAK_RT, kind: "number" },
    { name: "integral", id: SECTION_PEAK_INTEGRAL, kind: "number" },
    { name: "intensity", id: SECTION_PEAK_INTENSITY, kind: "number" },
    { name: "nPoints", id: SECTION_PEAK_POINT_COUNT, kind: "count" },
    { name: "noise", id: SECTION_PEAK_NOISE, kind: "number" },
    { name: "r2", id: SECTION_PEAK_R2, kind: "number" },
  ],
  [PAYLOAD_FEATURES]: [
    { name: "mz", id: SECTION_FEATURE_MZ, kind: "number" },
    { name: "rt", id: SECTION_FEATURE_RT, kind: "number" },
    { name: "from", id: SECTION_FEATURE_FROM, kind: "number" },
    { name: "to", id: SECTION_FEATURE_TO, kind: "number" },
    { name: "intensity", id: SECTION_FEATURE_INTENSITY, kind: "number" },
    { name: "integral", id: SECTION_FEATURE_INTEGRAL, kind: "number" },
    { name: "nPoints", id: SECTION_FEATURE_POINT_COUNT, kind: "count" },
    { name: "noise", id: SECTION_FEATURE_NOISE, kind: "number" },
  ],
  [PAYLOAD_CHROM_PEAKS]: [
    { name: "index", id: SECTION_CHROM_INDEX, kind: "count" },
    {
      name: "id",
      startsId: SECTION_CHROM_ID_STARTS,
      bytesId: SECTION_CHROM_ID_BYTES,
      kind: "text",
    },
    { name: "ort", id: SECTION_CHROM_TARGET_RT, kind: "number" },
    { name: "rt", id: SECTION_CHROM_RT, kind: "number" },
    { name: "from", id: SECTION_CHROM_FROM, kind: "number" },
    { name: "to", id: SECTION_CHROM_TO, kind: "number" },
    { name: "intensity", id: SECTION_CHROM_INTENSITY, kind: "number" },
    { name: "integral", id: SECTION_CHROM_INTEGRAL, kind: "number" },
    { name: "totalArea", id: SECTION_CHROM_TOTAL_AREA, kind: "number" },
    {
      name: "timestamp",
      startsId: SECTION_CHROM_TIMESTAMP_STARTS,
      bytesId: SECTION_CHROM_TIMESTAMP_BYTES,
      kind: "text",
    },
  ],
  [PAYLOAD_FIT_RESULT]: [
    { name: "shape", id: SECTION_FIT_SHAPE, kind: "number" },
    { name: "height", id: SECTION_FIT_HEIGHT, kind: "number" },
    { name: "center", id: SECTION_FIT_CENTER, kind: "number" },
    { name: "fwhm", id: SECTION_FIT_FWHM, kind: "number" },
    { name: "tail", id: SECTION_FIT_TAIL, kind: "number" },
    { name: "r2", id: SECTION_FIT_R2, kind: "number" },
  ],
  [PAYLOAD_EIC_PEAKS]: [
    {
      name: "id",
      startsId: SECTION_EIC_PEAK_ID_STARTS,
      bytesId: SECTION_EIC_PEAK_ID_BYTES,
      kind: "text",
    },
    { name: "mz", id: SECTION_EIC_PEAK_MZ, kind: "number" },
    { name: "ort", id: SECTION_EIC_PEAK_ORT, kind: "number" },
    { name: "rt", id: SECTION_EIC_PEAK_RT, kind: "number" },
    { name: "from", id: SECTION_EIC_PEAK_FROM, kind: "number" },
    { name: "to", id: SECTION_EIC_PEAK_TO, kind: "number" },
    { name: "intensity", id: SECTION_EIC_PEAK_INTENSITY, kind: "number" },
    { name: "integral", id: SECTION_EIC_PEAK_INTEGRAL, kind: "number" },
    { name: "noise", id: SECTION_EIC_PEAK_NOISE, kind: "number" },
  ],
  [PAYLOAD_CONSENSUS_FEATURES]: [
    { name: "mz", id: SECTION_CONSENSUS_MZ, kind: "number" },
    { name: "rt", id: SECTION_CONSENSUS_RT, kind: "number" },
    { name: "from", id: SECTION_CONSENSUS_FROM, kind: "number" },
    { name: "to", id: SECTION_CONSENSUS_TO, kind: "number" },
    { name: "intensity", id: SECTION_CONSENSUS_INTENSITY, kind: "number" },
    { name: "integral", id: SECTION_CONSENSUS_INTEGRAL, kind: "number" },
    { name: "frequency", id: SECTION_CONSENSUS_FREQUENCY, kind: "number" },
  ],
  [PAYLOAD_FOUND_FEATURES]: [
    {
      name: "id",
      startsId: SECTION_FOUND_ID_STARTS,
      bytesId: SECTION_FOUND_ID_BYTES,
      kind: "text",
    },
    { name: "mz", id: SECTION_FOUND_MZ, kind: "number" },
    { name: "rt", id: SECTION_FOUND_RT, kind: "number" },
    { name: "from", id: SECTION_FOUND_FROM, kind: "number" },
    { name: "to", id: SECTION_FOUND_TO, kind: "number" },
    { name: "intensity", id: SECTION_FOUND_INTENSITY, kind: "number" },
    { name: "integral", id: SECTION_FOUND_INTEGRAL, kind: "number" },
    { name: "nPoints", id: SECTION_FOUND_POINT_COUNT, kind: "count" },
    { name: "noise", id: SECTION_FOUND_NOISE, kind: "number" },
  ],
};

const decoder = new TextDecoder();

export function readRecords(buffer: ArrayBuffer): Record<string, unknown>[] {
  const bridge = readBridge(buffer);
  const columns = RECORD_COLUMNS[bridge.payloadKind];
  if (!columns) {
    throw new Error(
      `quantion bridge: payload kind ${bridge.payloadKind} is not a record table`,
    );
  }

  const count = bridge.recordCount;
  const values = new Map<string, (index: number) => unknown>();

  for (const column of columns) {
    if (column.kind === "number") {
      const numbers = readF64Section(bridge, column.id);
      values.set(column.name, (index) => numbers[index]);
    } else if (column.kind === "count") {
      const counts = readU32Section(bridge, column.id);
      values.set(column.name, (index) => counts[index]);
    } else {
      const starts = readCountSection(bridge, column.startsId);
      const bytes = readU8Section(bridge, column.bytesId);
      values.set(column.name, (index) =>
        decoder.decode(bytes.subarray(starts[index], starts[index + 1])),
      );
    }
  }

  const rows: Record<string, unknown>[] = new Array(count);
  for (let index = 0; index < count; index += 1) {
    const row: Record<string, unknown> = {};
    for (const [name, take] of values) row[name] = take(index);
    rows[index] = row;
  }
  return rows;
}
