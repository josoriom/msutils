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

    let header_path = artifacts_dir.join("quantion.h");

    cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(&header_path);

    let header_copy_paths = [
        root.join("wrappers/js/src/include/quantion.h"),
        root.join("wrappers/js/native/include/quantion.h"),
        root.join("wrappers/r/src/include/quantion.h"),
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
