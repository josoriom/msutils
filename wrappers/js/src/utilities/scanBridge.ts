import type { CentroidScan } from "../types/types";
import {
  PAYLOAD_SCANS,
  SECTION_BASE_PEAK_INT,
  SECTION_BASE_PEAK_MZ,
  SECTION_INTENSITY,
  SECTION_MS_LEVEL,
  SECTION_MZ,
  SECTION_POINT_STARTS,
  SECTION_POLARITY,
  SECTION_POSITION_X,
  SECTION_POSITION_Y,
  SECTION_POSITION_Z,
  SECTION_RT,
  SECTION_RT_SECONDS,
  SECTION_SELECTED_ION_MZ,
  SECTION_TOTAL_ION_CURRENT,
  expectKind,
  readCountSection,
  readF64Section,
  readBridge,
} from "./bridge";

export function readScans(buffer: ArrayBuffer): CentroidScan[] {
  const bridge = readBridge(buffer);
  expectKind(bridge, PAYLOAD_SCANS);

  const pointStarts = readCountSection(bridge, SECTION_POINT_STARTS);
  const mz = readF64Section(bridge, SECTION_MZ);
  const intensity = readF64Section(bridge, SECTION_INTENSITY);
  const rt = readF64Section(bridge, SECTION_RT);
  const rtSeconds = readF64Section(bridge, SECTION_RT_SECONDS);
  const basePeakMz = readF64Section(bridge, SECTION_BASE_PEAK_MZ);
  const selectedIonMz = readF64Section(bridge, SECTION_SELECTED_ION_MZ);
  const basePeakInt = readF64Section(bridge, SECTION_BASE_PEAK_INT);
  const totalIonCurrent = readF64Section(bridge, SECTION_TOTAL_ION_CURRENT);
  const msLevel = readF64Section(bridge, SECTION_MS_LEVEL);
  const polarity = readF64Section(bridge, SECTION_POLARITY);
  const positionX = readF64Section(bridge, SECTION_POSITION_X);
  const positionY = readF64Section(bridge, SECTION_POSITION_Y);
  const positionZ = readF64Section(bridge, SECTION_POSITION_Z);

  const count = bridge.recordCount;
  if (pointStarts.length !== count + 1) {
    throw new Error("quantion bridge: point starts do not match the scan count");
  }
  if (intensity.length !== mz.length) {
    throw new Error(
      "quantion bridge: intensity length does not match m/z length",
    );
  }
  if (pointStarts[count] !== mz.length) {
    throw new Error("quantion bridge: point starts do not span the m/z section");
  }

  const scans: CentroidScan[] = new Array(count);
  for (let index = 0; index < count; index += 1) {
    const from = pointStarts[index];
    const to = pointStarts[index + 1];
    scans[index] = {
      rt: rt[index],
      mz: mz.subarray(from, to),
      intensity: intensity.subarray(from, to),
      metadata: {
        rtSeconds: rtSeconds[index],
        basePeakMz: basePeakMz[index],
        selectedIonMz: selectedIonMz[index],
        basePeakInt: basePeakInt[index],
        totalIonCurrent: totalIonCurrent[index],
        msLevel: msLevel[index],
        polarity: polarity[index],
        positionX: positionX[index],
        positionY: positionY[index],
        positionZ: positionZ[index],
      },
    };
  }
  return scans;
}
