use std::slice;

const _: () = assert!(cfg!(target_endian = "little"));

pub const QUANTION_BRIDGE_MAGIC: u32 = 0x42544E51;
pub const QUANTION_BRIDGE_LAYOUT_VERSION: u16 = 1;
pub const QUANTION_BRIDGE_HEADER_BYTES: u64 = 32;
pub const QUANTION_SECTION_ENTRY_BYTES: u64 = 24;

pub const QUANTION_PAYLOAD_SCANS: u16 = 1;
pub const QUANTION_PAYLOAD_ION_IMAGE: u16 = 2;
pub const QUANTION_PAYLOAD_PEAKS: u16 = 3;
pub const QUANTION_PAYLOAD_FEATURES: u16 = 4;
pub const QUANTION_PAYLOAD_CHROM_PEAKS: u16 = 5;

pub const QUANTION_PAYLOAD_FIT_RESULT: u16 = 6;
pub const QUANTION_PAYLOAD_EIC_PEAKS: u16 = 7;
pub const QUANTION_PAYLOAD_CONSENSUS_FEATURES: u16 = 8;
pub const QUANTION_PAYLOAD_FOUND_FEATURES: u16 = 9;

pub const QUANTION_ELEMENT_F64: u32 = 1;
pub const QUANTION_ELEMENT_U32: u32 = 2;
pub const QUANTION_ELEMENT_U64: u32 = 3;
pub const QUANTION_ELEMENT_U8: u32 = 4;

pub const QUANTION_SECTION_POINT_STARTS: u32 = 1;
pub const QUANTION_SECTION_MZ: u32 = 2;
pub const QUANTION_SECTION_INTENSITY: u32 = 3;
pub const QUANTION_SECTION_RT: u32 = 4;
pub const QUANTION_SECTION_RT_SECONDS: u32 = 5;
pub const QUANTION_SECTION_BASE_PEAK_MZ: u32 = 6;
pub const QUANTION_SECTION_SELECTED_ION_MZ: u32 = 7;
pub const QUANTION_SECTION_BASE_PEAK_INT: u32 = 8;
pub const QUANTION_SECTION_TOTAL_ION_CURRENT: u32 = 9;
pub const QUANTION_SECTION_MS_LEVEL: u32 = 10;
pub const QUANTION_SECTION_POLARITY: u32 = 11;
pub const QUANTION_SECTION_POSITION_X: u32 = 12;
pub const QUANTION_SECTION_POSITION_Y: u32 = 13;
pub const QUANTION_SECTION_POSITION_Z: u32 = 14;

pub const QUANTION_SECTION_IMAGE_SHAPE: u32 = 15;
pub const QUANTION_SECTION_IMAGE_DATA: u32 = 16;
pub const QUANTION_SECTION_IMAGE_COUNTS: u32 = 17;

pub const QUANTION_SECTION_PEAK_FROM: u32 = 18;
pub const QUANTION_SECTION_PEAK_TO: u32 = 19;
pub const QUANTION_SECTION_PEAK_RT: u32 = 20;
pub const QUANTION_SECTION_PEAK_INTEGRAL: u32 = 21;
pub const QUANTION_SECTION_PEAK_INTENSITY: u32 = 22;
pub const QUANTION_SECTION_PEAK_POINT_COUNT: u32 = 23;
pub const QUANTION_SECTION_PEAK_NOISE: u32 = 24;
pub const QUANTION_SECTION_PEAK_R2: u32 = 25;

pub const QUANTION_SECTION_FEATURE_MZ: u32 = 26;
pub const QUANTION_SECTION_FEATURE_RT: u32 = 27;
pub const QUANTION_SECTION_FEATURE_FROM: u32 = 28;
pub const QUANTION_SECTION_FEATURE_TO: u32 = 29;
pub const QUANTION_SECTION_FEATURE_INTENSITY: u32 = 30;
pub const QUANTION_SECTION_FEATURE_INTEGRAL: u32 = 31;
pub const QUANTION_SECTION_FEATURE_POINT_COUNT: u32 = 32;
pub const QUANTION_SECTION_FEATURE_NOISE: u32 = 33;

pub const QUANTION_SECTION_CHROM_INDEX: u32 = 34;
pub const QUANTION_SECTION_CHROM_TARGET_RT: u32 = 35;
pub const QUANTION_SECTION_CHROM_RT: u32 = 36;
pub const QUANTION_SECTION_CHROM_FROM: u32 = 37;
pub const QUANTION_SECTION_CHROM_TO: u32 = 38;
pub const QUANTION_SECTION_CHROM_INTENSITY: u32 = 39;
pub const QUANTION_SECTION_CHROM_INTEGRAL: u32 = 40;
pub const QUANTION_SECTION_CHROM_TOTAL_AREA: u32 = 41;
pub const QUANTION_SECTION_CHROM_ID_STARTS: u32 = 42;
pub const QUANTION_SECTION_CHROM_ID_BYTES: u32 = 43;
pub const QUANTION_SECTION_CHROM_TIMESTAMP_STARTS: u32 = 44;
pub const QUANTION_SECTION_CHROM_TIMESTAMP_BYTES: u32 = 45;

pub const QUANTION_SECTION_FIT_SHAPE: u32 = 46;
pub const QUANTION_SECTION_FIT_HEIGHT: u32 = 47;
pub const QUANTION_SECTION_FIT_CENTER: u32 = 48;
pub const QUANTION_SECTION_FIT_FWHM: u32 = 49;
pub const QUANTION_SECTION_FIT_TAIL: u32 = 50;
pub const QUANTION_SECTION_FIT_R2: u32 = 51;

pub const QUANTION_SECTION_EIC_PEAK_ID_STARTS: u32 = 52;
pub const QUANTION_SECTION_EIC_PEAK_ID_BYTES: u32 = 53;
pub const QUANTION_SECTION_EIC_PEAK_MZ: u32 = 54;
pub const QUANTION_SECTION_EIC_PEAK_ORT: u32 = 55;
pub const QUANTION_SECTION_EIC_PEAK_RT: u32 = 56;
pub const QUANTION_SECTION_EIC_PEAK_FROM: u32 = 57;
pub const QUANTION_SECTION_EIC_PEAK_TO: u32 = 58;
pub const QUANTION_SECTION_EIC_PEAK_INTENSITY: u32 = 59;
pub const QUANTION_SECTION_EIC_PEAK_INTEGRAL: u32 = 60;
pub const QUANTION_SECTION_EIC_PEAK_NOISE: u32 = 61;

pub const QUANTION_SECTION_CONSENSUS_MZ: u32 = 62;
pub const QUANTION_SECTION_CONSENSUS_RT: u32 = 63;
pub const QUANTION_SECTION_CONSENSUS_FROM: u32 = 64;
pub const QUANTION_SECTION_CONSENSUS_TO: u32 = 65;
pub const QUANTION_SECTION_CONSENSUS_INTENSITY: u32 = 66;
pub const QUANTION_SECTION_CONSENSUS_INTEGRAL: u32 = 67;
pub const QUANTION_SECTION_CONSENSUS_FREQUENCY: u32 = 68;

pub const QUANTION_SECTION_FOUND_ID_STARTS: u32 = 69;
pub const QUANTION_SECTION_FOUND_ID_BYTES: u32 = 70;
pub const QUANTION_SECTION_FOUND_MZ: u32 = 71;
pub const QUANTION_SECTION_FOUND_RT: u32 = 72;
pub const QUANTION_SECTION_FOUND_FROM: u32 = 73;
pub const QUANTION_SECTION_FOUND_TO: u32 = 74;
pub const QUANTION_SECTION_FOUND_INTENSITY: u32 = 75;
pub const QUANTION_SECTION_FOUND_INTEGRAL: u32 = 76;
pub const QUANTION_SECTION_FOUND_POINT_COUNT: u32 = 77;
pub const QUANTION_SECTION_FOUND_NOISE: u32 = 78;

const ALIGNMENT: u64 = 8;

fn size_of_element(element_type: u32) -> Option<u64> {
    match element_type {
        QUANTION_ELEMENT_F64 => Some(8),
        QUANTION_ELEMENT_U32 => Some(4),
        QUANTION_ELEMENT_U64 => Some(8),
        QUANTION_ELEMENT_U8 => Some(1),
        _ => None,
    }
}

fn round_up(value: u64) -> Option<u64> {
    value.checked_add(ALIGNMENT - 1).map(|sum| sum & !(ALIGNMENT - 1))
}

#[derive(Clone, Copy)]
struct Section {
    id: u32,
    element_type: u32,
    element_count: u64,
    byte_offset: u64,
    byte_length: u64,
}

pub struct BridgeBuilder {
    payload_kind: u16,
    record_count: u64,
    sections: Vec<Section>,
}

impl BridgeBuilder {
    pub fn new(payload_kind: u16, record_count: u64) -> Self {
        Self {
            payload_kind,
            record_count,
            sections: Vec::new(),
        }
    }

    pub fn add_section(&mut self, id: u32, element_type: u32, element_count: u64) {
        self.sections.push(Section {
            id,
            element_type,
            element_count,
            byte_offset: 0,
            byte_length: 0,
        });
    }

    pub fn build(mut self) -> Option<Bridge> {
        let section_count = u32::try_from(self.sections.len()).ok()?;
        let table_bytes = QUANTION_SECTION_ENTRY_BYTES.checked_mul(u64::from(section_count))?;
        let mut next_offset = QUANTION_BRIDGE_HEADER_BYTES.checked_add(table_bytes)?;

        for section in self.sections.iter_mut() {
            let element_size = size_of_element(section.element_type)?;
            let byte_length = section.element_count.checked_mul(element_size)?;
            section.byte_offset = next_offset;
            section.byte_length = byte_length;
            next_offset = round_up(next_offset.checked_add(byte_length)?)?;
        }

        let total_bytes = next_offset;
        let total = usize::try_from(total_bytes).ok()?;
        let mut bytes = vec![0u8; total];

        bytes[0..4].copy_from_slice(&QUANTION_BRIDGE_MAGIC.to_le_bytes());
        bytes[4..6].copy_from_slice(&QUANTION_BRIDGE_LAYOUT_VERSION.to_le_bytes());
        bytes[6..8].copy_from_slice(&self.payload_kind.to_le_bytes());
        bytes[8..12].copy_from_slice(&section_count.to_le_bytes());
        bytes[12..16]
            .copy_from_slice(&(QUANTION_BRIDGE_HEADER_BYTES as u32).to_le_bytes());
        bytes[16..24].copy_from_slice(&total_bytes.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.record_count.to_le_bytes());

        for (index, section) in self.sections.iter().enumerate() {
            let start = (QUANTION_BRIDGE_HEADER_BYTES as usize)
                + index * (QUANTION_SECTION_ENTRY_BYTES as usize);
            bytes[start..start + 4].copy_from_slice(&section.id.to_le_bytes());
            bytes[start + 4..start + 8].copy_from_slice(&section.element_type.to_le_bytes());
            bytes[start + 8..start + 16].copy_from_slice(&section.byte_offset.to_le_bytes());
            bytes[start + 16..start + 24].copy_from_slice(&section.byte_length.to_le_bytes());
        }

        Some(Bridge {
            bytes,
            sections: self.sections,
        })
    }
}

pub struct Bridge {
    bytes: Vec<u8>,
    sections: Vec<Section>,
}

impl Bridge {
    fn find_section(&self, id: u32, element_type: u32) -> Section {
        let section = self
            .sections
            .iter()
            .find(|section| section.id == id)
            .expect("bridge section was never declared");
        assert_eq!(
            section.element_type, element_type,
            "bridge section written with the wrong element type"
        );
        *section
    }

    fn copy_into(&mut self, byte_offset: u64, source: &[u8]) {
        let start = byte_offset as usize;
        self.bytes[start..start + source.len()].copy_from_slice(source);
    }

    pub fn write_f64_at(&mut self, id: u32, element_offset: u64, values: &[f64]) {
        let section = self.find_section(id, QUANTION_ELEMENT_F64);
        let end = element_offset + values.len() as u64;
        assert!(end <= section.element_count, "bridge f64 write out of range");
        let source = unsafe {
            slice::from_raw_parts(values.as_ptr() as *const u8, std::mem::size_of_val(values))
        };
        self.copy_into(section.byte_offset + element_offset * 8, source);
    }

    pub fn write_f64_section(&mut self, id: u32, values: &[f64]) {
        let section = self.find_section(id, QUANTION_ELEMENT_F64);
        assert_eq!(
            section.element_count,
            values.len() as u64,
            "bridge f64 section written with the wrong element count"
        );
        self.write_f64_at(id, 0, values);
    }

    pub fn write_u32_section(&mut self, id: u32, values: &[u32]) {
        let section = self.find_section(id, QUANTION_ELEMENT_U32);
        assert_eq!(
            section.element_count,
            values.len() as u64,
            "bridge u32 section written with the wrong element count"
        );
        let source = unsafe {
            slice::from_raw_parts(values.as_ptr() as *const u8, std::mem::size_of_val(values))
        };
        self.copy_into(section.byte_offset, source);
    }

    pub fn write_u64_section(&mut self, id: u32, values: &[u64]) {
        let section = self.find_section(id, QUANTION_ELEMENT_U64);
        assert_eq!(
            section.element_count,
            values.len() as u64,
            "bridge u64 section written with the wrong element count"
        );
        let source = unsafe {
            slice::from_raw_parts(values.as_ptr() as *const u8, std::mem::size_of_val(values))
        };
        self.copy_into(section.byte_offset, source);
    }

    pub fn write_text_section(&mut self, id: u32, values: &[&str]) {
        let section = self.find_section(id, QUANTION_ELEMENT_U8);
        let mut at = section.byte_offset;
        for text in values {
            let bytes = text.as_bytes();
            self.copy_into(at, bytes);
            at += bytes.len() as u64;
        }
        assert_eq!(
            at - section.byte_offset,
            section.byte_length,
            "bridge text section written with the wrong byte count"
        );
    }

    pub fn into_bytes(self) -> Box<[u8]> {
        debug_assert_eq!(
            self.bytes.capacity(),
            self.bytes.len(),
            "bridge allocation must not carry spare capacity"
        );
        self.bytes.into_boxed_slice()
    }
}

pub enum Column<'a> {
    Numbers { id: u32, values: &'a [f64] },
    Counts { id: u32, values: &'a [u32] },
    Text { starts_id: u32, bytes_id: u32, values: &'a [&'a str] },
}

pub fn build_record_bridge(
    payload_kind: u16,
    record_count: u64,
    columns: &[Column],
) -> Option<Box<[u8]>> {
    let mut builder = BridgeBuilder::new(payload_kind, record_count);
    for column in columns {
        match column {
            Column::Numbers { id, values } => {
                builder.add_section(*id, QUANTION_ELEMENT_F64, values.len() as u64)
            }
            Column::Counts { id, values } => {
                builder.add_section(*id, QUANTION_ELEMENT_U32, values.len() as u64)
            }
            Column::Text {
                starts_id,
                bytes_id,
                values,
            } => {
                let mut total: u64 = 0;
                for text in values.iter() {
                    total = total.checked_add(text.len() as u64)?;
                }
                builder.add_section(*starts_id, QUANTION_ELEMENT_U64, values.len() as u64 + 1);
                builder.add_section(*bytes_id, QUANTION_ELEMENT_U8, total);
            }
        }
    }
    let mut bridge = builder.build()?;

    for column in columns {
        match column {
            Column::Numbers { id, values } => bridge.write_f64_section(*id, values),
            Column::Counts { id, values } => bridge.write_u32_section(*id, values),
            Column::Text {
                starts_id,
                bytes_id,
                values,
            } => {
                let mut starts = Vec::with_capacity(values.len() + 1);
                let mut next: u64 = 0;
                starts.push(0u64);
                for text in values.iter() {
                    next += text.len() as u64;
                    starts.push(next);
                }
                bridge.write_u64_section(*starts_id, &starts);
                bridge.write_text_section(*bytes_id, values);
            }
        }
    }
    Some(bridge.into_bytes())
}
