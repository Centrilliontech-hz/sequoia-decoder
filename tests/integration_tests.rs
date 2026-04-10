use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

// Configuration for datasets, mirroring the logic in run_tests.py
struct DatasetConfig {
    sep: &'static str,
    rc: i32,
}

fn get_config(dataset_name: &str) -> DatasetConfig {
    match dataset_name {
        "2026-01-28T05-19-17Z" => DatasetConfig {
            sep: "HH,GGG,HH",
            rc: -1,
        },
        "NE-CA-251011-2-6G" => DatasetConfig {
            sep: "HH,GGG,HH",
            rc: -1,
        },
        "CA2-251026-TCR-5G" => DatasetConfig {
            sep: "DD,CCC,DD",
            rc: 1,
        },
        _ => DatasetConfig {
            sep: "HH,GGG,HH",
            rc: -1,
        },
    }
}

fn ensure_release_binary() -> PathBuf {
    let root = env::current_dir().unwrap();
    let release_bin = root.join("target/release/sequoia-decoder");

    println!("Building release binary...");
    let status = Command::new("cargo")
        .arg("build")
        .arg("--release")
        .status()
        .expect("Failed to run cargo build --release");
    assert!(status.success(), "Cargo build --release failed");

    release_bin
}

#[test]
fn test_datasets() {
    let binary = ensure_release_binary();
    let root = env::current_dir().unwrap();
    let assets_dir = root.join("tests/assets");
    let output_dir = root.join("tests/output_rust");

    if output_dir.exists() {
        fs::remove_dir_all(&output_dir).unwrap();
    }
    fs::create_dir_all(&output_dir).unwrap();

    let entries = fs::read_dir(&assets_dir).expect("Failed to read tests/assets");

    for entry in entries {
        let entry = entry.unwrap();
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }

        let dataset_name = path.file_name().unwrap().to_str().unwrap();
        println!("Testing dataset: {}", dataset_name);

        // Find R1 and R2 files
        let mut r1_path = None;
        let mut r2_path = None;

        for file in fs::read_dir(&path).unwrap() {
            let file_path = file.unwrap().path();
            let name = file_path.file_name().unwrap().to_str().unwrap();
            if name.contains("R1") && name.ends_with(".fastq.gz") {
                r1_path = Some(file_path.clone());
            } else if name.contains("R2") && name.ends_with(".fastq.gz") {
                r2_path = Some(file_path);
            }
        }

        if r1_path.is_none() || r2_path.is_none() {
            println!("Skipping {}: R1 or R2 not found", dataset_name);
            continue;
        }

        let r1 = r1_path.unwrap();
        let r2 = r2_path.unwrap();
        let config = get_config(dataset_name);

        let dataset_out_dir = output_dir.join(dataset_name);
        fs::create_dir_all(&dataset_out_dir).unwrap();

        // Construct command
        // Note: Using the new CLI arguments format
        let status = Command::new(&binary)
            .arg("-i")
            .arg(&r1)
            .arg(&r2)
            .arg("--output-dir")
            .arg(&dataset_out_dir)
            .arg("--sep")
            .arg(config.sep)
            .arg("--rc")
            .arg(config.rc.to_string())
            .arg("--threads")
            .arg("4")
            .arg("--head")
            .arg("200000")
            // .arg("--codebook") // Optional now, using built-in default
            .status()
            .expect("Failed to execute process");

        assert!(
            status.success(),
            "Decoder failed for dataset {}",
            dataset_name
        );

        // Verify outputs (names are now fixed in main.rs)
        let output_table = dataset_out_dir.join("decoded_zipcodes.csv.gz");
        let output_r2 = dataset_out_dir.join("R2.fastq.gz");
        let json_summary = dataset_out_dir.join("stats.json");

        assert!(
            output_table.exists(),
            "Output table missing for {}",
            dataset_name
        );
        assert!(output_r2.exists(), "Output R2 missing for {}", dataset_name);
        assert!(
            json_summary.exists(),
            "Summary JSON missing for {}",
            dataset_name
        );

        // Simple content check
        let metadata = fs::metadata(&output_table).unwrap();
        assert!(
            metadata.len() > 0,
            "Output table is empty for {}",
            dataset_name
        );

        println!("Successfully processed {}", dataset_name);
    }
}
