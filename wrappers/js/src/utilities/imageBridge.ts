import type { IonImage } from "../types/types";
import {
  PAYLOAD_ION_IMAGE,
  SECTION_IMAGE_COUNTS,
  SECTION_IMAGE_DATA,
  SECTION_IMAGE_SHAPE,
  expectKind,
  readF64Section,
  readBridge,
  readU32Section,
} from "./bridge";

export function readIonImage(buffer: ArrayBuffer): IonImage {
  const bridge = readBridge(buffer);
  expectKind(bridge, PAYLOAD_ION_IMAGE);

  const shape = readU32Section(bridge, SECTION_IMAGE_SHAPE);
  if (shape.length !== 6) {
    throw new Error("quantion bridge: image shape must hold six values");
  }

  const data = readF64Section(bridge, SECTION_IMAGE_DATA);
  const counts = readU32Section(bridge, SECTION_IMAGE_COUNTS);
  const [width, height, minX, minY, minZ, maxZ] = shape;

  if (width * height !== bridge.recordCount) {
    throw new Error(
      "quantion bridge: image size does not match the record count",
    );
  }
  if (
    data.length !== bridge.recordCount ||
    counts.length !== bridge.recordCount
  ) {
    throw new Error(
      "quantion bridge: image sections do not match the record count",
    );
  }

  return { width, height, minX, minY, minZ, maxZ, data, counts };
}
