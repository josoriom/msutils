fn main() {
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let root = std::path::PathBuf::from(&crate_dir)
        .parent()
        .unwrap()
        .to_path_buf();

    let artifacts_dir = root.join("artifacts/include");
    std::fs::create_dir_all(&artifacts_dir).unwrap();

    let config =
        cbindgen::Config::from_file(std::path::PathBuf::from(&crate_dir).join("cbindgen.toml"))
            .unwrap();

    let header_path = artifacts_dir.join("msutils.h");

    cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(&header_path);

    let ffi_declarations = "\nextern int32_t parse_mzml(const uint8_t *data, uintptr_t len, void **out);\n\
extern int32_t parse_bin(const uint8_t *data, uintptr_t len, uintptr_t cache_bytes, void **out);\n\
extern int32_t parse_ion_path(const char *path, uintptr_t cache_bytes, void **out);\n\
extern int32_t parse_ion_url(const char *url, uintptr_t cache_bytes, void **out);\n\
extern int32_t plan_open(const uint8_t *header_ptr, uintptr_t header_len, struct Buf *out);\n\
extern int32_t plan_eic(void *h, double target, double from, double to, double ppm, double mz_tol, struct Buf *out);\n\
extern int32_t calculate_eic(void *h, double target, double from, double to, double ppm, double mz_tol, struct Buf *out_x, struct Buf *out_y);\n\
extern int32_t bin_to_json(void *h, struct Buf *out);\n\
extern int32_t bin_to_mzml(void *h, struct Buf *out);\n\
extern int32_t mzml_to_bin(void *h, struct Buf *out, uint8_t level, uint8_t compress);\n\
extern void dispose_mzml(void *h);\n\
extern void free_(uint8_t *ptr, uintptr_t len);\n";

    let header_content = std::fs::read_to_string(&header_path).unwrap();
    let updated_header = header_content + ffi_declarations;
    std::fs::write(&header_path, updated_header).unwrap();

    let header_copy_paths = [
        root.join("wrappers/js/src/include/msutils.h"),
        root.join("wrappers/js/native/include/msutils.h"),
        root.join("wrappers/r/src/include/msutils.h"),
    ];

    for header_copy_path in &header_copy_paths {
        std::fs::create_dir_all(header_copy_path.parent().unwrap()).unwrap();
        std::fs::copy(&header_path, header_copy_path)
            .unwrap_or_else(|_| panic!("failed to copy {}", header_copy_path.display()));
        println!("cargo:rerun-if-changed={}", header_copy_path.display());
    }

    println!("cargo:rerun-if-changed={}", header_path.display());
    println!("cargo:rerun-if-changed=cbindgen.toml");
}
