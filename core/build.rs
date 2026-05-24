fn main() {
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let root = std::path::PathBuf::from(&crate_dir).parent().unwrap().to_path_buf();

    let artifacts_dir = root.join("artifacts/include");
    std::fs::create_dir_all(&artifacts_dir).unwrap();

    let config = cbindgen::Config::from_file(
        std::path::PathBuf::from(&crate_dir).join("cbindgen.toml"),
    )
    .unwrap();

    cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(config)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(artifacts_dir.join("msutils.h"));

    // Sync the generated header into every wrapper that needs it.
    let src = artifacts_dir.join("msutils.h");
    let destinations = [
        root.join("wrappers/js/src/include/msutils.h"),
        root.join("wrappers/r/src/include/msutils.h"),
    ];
    for dst in &destinations {
        std::fs::create_dir_all(dst.parent().unwrap()).unwrap();
        std::fs::copy(&src, dst).expect(&format!("Failed to copy msutils.h to {}", dst.display()));
        // Tell cargo to re-run this script if a destination is deleted or
        // edited out-of-band, so the wrapper headers can never be stale.
        println!("cargo:rerun-if-changed={}", dst.display());
    }
    println!("cargo:rerun-if-changed={}", src.display());
    println!("cargo:rerun-if-changed=cbindgen.toml");
}
