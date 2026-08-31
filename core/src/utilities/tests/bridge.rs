#[cfg(test)]
mod tests {
    use crate::utilities::bridge::*;

    fn read_u16(bytes: &[u8], at: usize) -> u16 {
        u16::from_le_bytes(bytes[at..at + 2].try_into().unwrap())
    }

    fn read_u32(bytes: &[u8], at: usize) -> u32 {
        u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap())
    }

    fn read_u64(bytes: &[u8], at: usize) -> u64 {
        u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap())
    }

    fn find_section(bytes: &[u8], wanted: u32) -> (u32, u64, u64) {
        let section_count = read_u32(bytes, 8);
        for index in 0..section_count as usize {
            let start = 32 + index * 24;
            if read_u32(bytes, start) == wanted {
                return (
                    read_u32(bytes, start + 4),
                    read_u64(bytes, start + 8),
                    read_u64(bytes, start + 16),
                );
            }
        }
        panic!("section {wanted} is missing");
    }

    fn read_f64_section(bytes: &[u8], wanted: u32) -> Vec<f64> {
        let (element_type, byte_offset, byte_length) = find_section(bytes, wanted);
        assert_eq!(element_type, QUANTION_ELEMENT_F64);
        let start = byte_offset as usize;
        bytes[start..start + byte_length as usize]
            .chunks_exact(8)
            .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
            .collect()
    }

    fn build_two_scans() -> Box<[u8]> {
        let mut builder = BridgeBuilder::new(QUANTION_PAYLOAD_SCANS, 2);
        builder.add_section(QUANTION_SECTION_POINT_STARTS, QUANTION_ELEMENT_U64, 3);
        builder.add_section(QUANTION_SECTION_MZ, QUANTION_ELEMENT_F64, 5);
        builder.add_section(QUANTION_SECTION_INTENSITY, QUANTION_ELEMENT_F64, 5);
        builder.add_section(QUANTION_SECTION_RT, QUANTION_ELEMENT_F64, 2);
        builder.add_section(QUANTION_SECTION_BASE_PEAK_MZ, QUANTION_ELEMENT_F64, 2);
        let mut bridge = builder.build().expect("bridge builds");

        bridge.write_f64_at(QUANTION_SECTION_MZ, 0, &[100.0, 200.0]);
        bridge.write_f64_at(QUANTION_SECTION_MZ, 2, &[300.0, 400.0, 500.0]);
        bridge.write_f64_at(QUANTION_SECTION_INTENSITY, 0, &[1.0, 2.0]);
        bridge.write_f64_at(QUANTION_SECTION_INTENSITY, 2, &[3.0, 4.0, 5.0]);
        bridge.write_f64_at(QUANTION_SECTION_RT, 0, &[0.5]);
        bridge.write_f64_at(QUANTION_SECTION_RT, 1, &[1.5]);
        bridge.write_f64_at(QUANTION_SECTION_BASE_PEAK_MZ, 0, &[f64::NAN]);
        bridge.write_f64_at(QUANTION_SECTION_BASE_PEAK_MZ, 1, &[-0.0]);
        bridge.write_u64_section(QUANTION_SECTION_POINT_STARTS, &[0, 2, 5]);
        bridge.into_bytes()
    }

    #[test]
    fn test_header_matches_the_wire_format() {
        let bytes = build_two_scans();
        assert_eq!(read_u32(&bytes, 0), QUANTION_BRIDGE_MAGIC);
        assert_eq!(read_u16(&bytes, 4), QUANTION_BRIDGE_LAYOUT_VERSION);
        assert_eq!(read_u16(&bytes, 6), QUANTION_PAYLOAD_SCANS);
        assert_eq!(read_u32(&bytes, 8), 5);
        assert_eq!(read_u32(&bytes, 12), 32);
        assert_eq!(read_u64(&bytes, 16), bytes.len() as u64);
        assert_eq!(read_u64(&bytes, 24), 2);
    }

    #[test]
    fn test_every_section_starts_on_a_multiple_of_eight() {
        let bytes = build_two_scans();
        let section_count = read_u32(&bytes, 8);
        for index in 0..section_count as usize {
            let start = 32 + index * 24;
            let byte_offset = read_u64(&bytes, start + 8);
            assert_eq!(byte_offset % 8, 0, "section {index} is not aligned");
        }
        assert_eq!(bytes.len() % 8, 0, "total is not rounded to eight");
    }

    #[test]
    fn test_sections_stay_inside_the_bridge_and_do_not_overlap() {
        let bytes = build_two_scans();
        let total = bytes.len() as u64;
        let section_count = read_u32(&bytes, 8);
        let mut reach = 32 + u64::from(section_count) * 24;
        for index in 0..section_count as usize {
            let start = 32 + index * 24;
            let byte_offset = read_u64(&bytes, start + 8);
            let byte_length = read_u64(&bytes, start + 16);
            assert!(byte_offset >= reach, "section {index} overlaps the previous");
            assert!(byte_offset <= total);
            assert!(byte_length <= total - byte_offset);
            reach = byte_offset + byte_length;
        }
    }

    #[test]
    fn test_values_survive_the_round_trip() {
        let bytes = build_two_scans();
        assert_eq!(
            read_f64_section(&bytes, QUANTION_SECTION_MZ),
            vec![100.0, 200.0, 300.0, 400.0, 500.0]
        );
        assert_eq!(
            read_f64_section(&bytes, QUANTION_SECTION_INTENSITY),
            vec![1.0, 2.0, 3.0, 4.0, 5.0]
        );
        assert_eq!(read_f64_section(&bytes, QUANTION_SECTION_RT), vec![0.5, 1.5]);
    }

    #[test]
    fn test_not_a_number_and_negative_zero_survive() {
        let bytes = build_two_scans();
        let column = read_f64_section(&bytes, QUANTION_SECTION_BASE_PEAK_MZ);
        assert!(column[0].is_nan(), "NaN must not become zero");
        assert!(
            column[1] == 0.0 && column[1].is_sign_negative(),
            "negative zero must keep its sign"
        );
    }

    #[test]
    fn test_point_starts_describe_every_scan() {
        let bytes = build_two_scans();
        let (element_type, byte_offset, byte_length) =
            find_section(&bytes, QUANTION_SECTION_POINT_STARTS);
        assert_eq!(element_type, QUANTION_ELEMENT_U64);
        assert_eq!(byte_length, 3 * 8);
        let starts: Vec<u64> = (0..3)
            .map(|index| read_u64(&bytes, byte_offset as usize + index * 8))
            .collect();
        assert_eq!(starts, vec![0, 2, 5]);
        let (_, _, mz_bytes) = find_section(&bytes, QUANTION_SECTION_MZ);
        assert_eq!(starts[2], mz_bytes / 8);
    }

    #[test]
    fn test_empty_bridge_is_readable() {
        let mut builder = BridgeBuilder::new(QUANTION_PAYLOAD_SCANS, 0);
        builder.add_section(QUANTION_SECTION_POINT_STARTS, QUANTION_ELEMENT_U64, 1);
        builder.add_section(QUANTION_SECTION_MZ, QUANTION_ELEMENT_F64, 0);
        let mut bridge = builder.build().expect("empty bridge builds");
        bridge.write_u64_section(QUANTION_SECTION_POINT_STARTS, &[0]);
        let bytes = bridge.into_bytes();

        assert_eq!(read_u64(&bytes, 24), 0);
        let (_, byte_offset, byte_length) = find_section(&bytes, QUANTION_SECTION_MZ);
        assert_eq!(byte_length, 0);
        assert!(byte_offset <= bytes.len() as u64);
    }

    #[test]
    fn test_image_shape_round_trips() {
        let mut builder = BridgeBuilder::new(QUANTION_PAYLOAD_ION_IMAGE, 4);
        builder.add_section(QUANTION_SECTION_IMAGE_SHAPE, QUANTION_ELEMENT_U32, 6);
        builder.add_section(QUANTION_SECTION_IMAGE_DATA, QUANTION_ELEMENT_F64, 4);
        builder.add_section(QUANTION_SECTION_IMAGE_COUNTS, QUANTION_ELEMENT_U32, 4);
        let mut bridge = builder.build().expect("image bridge builds");
        bridge.write_u32_section(QUANTION_SECTION_IMAGE_SHAPE, &[2, 2, 0, 0, 0, 1]);
        bridge.write_f64_section(QUANTION_SECTION_IMAGE_DATA, &[1.0, 2.0, 3.0, 4.0]);
        bridge.write_u32_section(QUANTION_SECTION_IMAGE_COUNTS, &[1, 1, 2, 2]);
        let bytes = bridge.into_bytes();

        assert_eq!(read_u16(&bytes, 6), QUANTION_PAYLOAD_ION_IMAGE);
        let (_, shape_offset, _) = find_section(&bytes, QUANTION_SECTION_IMAGE_SHAPE);
        let width = read_u32(&bytes, shape_offset as usize);
        let height = read_u32(&bytes, shape_offset as usize + 4);
        assert_eq!(width * height, read_u64(&bytes, 24) as u32);
        assert_eq!(
            read_f64_section(&bytes, QUANTION_SECTION_IMAGE_DATA),
            vec![1.0, 2.0, 3.0, 4.0]
        );
    }

    #[test]
    fn test_section_ids_never_collide_across_kinds() {
        let scans = [
            QUANTION_SECTION_POINT_STARTS,
            QUANTION_SECTION_MZ,
            QUANTION_SECTION_INTENSITY,
            QUANTION_SECTION_RT,
            QUANTION_SECTION_RT_SECONDS,
            QUANTION_SECTION_BASE_PEAK_MZ,
            QUANTION_SECTION_SELECTED_ION_MZ,
            QUANTION_SECTION_BASE_PEAK_INT,
            QUANTION_SECTION_TOTAL_ION_CURRENT,
            QUANTION_SECTION_MS_LEVEL,
            QUANTION_SECTION_POLARITY,
            QUANTION_SECTION_POSITION_X,
            QUANTION_SECTION_POSITION_Y,
            QUANTION_SECTION_POSITION_Z,
        ];
        let image = [
            QUANTION_SECTION_IMAGE_SHAPE,
            QUANTION_SECTION_IMAGE_DATA,
            QUANTION_SECTION_IMAGE_COUNTS,
        ];
        for scan_id in scans {
            for image_id in image {
                assert_ne!(scan_id, image_id, "section ids must be unique across kinds");
            }
        }
    }
}
