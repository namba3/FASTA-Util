#[path = "../benches/support/mod.rs"]
mod support;

use fasta_util::nucleic_acid::NUCLEIC_ACID_SET;
use std::{
    path::Path,
    process::{Command, Output},
    time::Duration,
};
use support::{
    DEFAULT_INPUT_SIZE, DEFAULT_SAMPLE_MS, INVALID_BASES, mixed_sequence, parse_config,
    parse_config_for, repeated_bytes,
};

fn args(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_owned()).collect()
}

fn run_analysis_benchmark(args: &[&str]) -> Output {
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("scripts/benchmark_analysis.sh");
    Command::new(script)
        .args(args)
        .output()
        .expect("failed to run analysis benchmark script")
}

#[test]
fn analysis_benchmark_help_does_not_require_benchmark_dependencies() {
    let output = run_analysis_benchmark(&["--help"]);

    assert!(output.status.success());
    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("Usage: benchmark_analysis.sh [BASES] [RUNS]"));
    assert!(help.contains("BASES must be at least 64"));
}

#[test]
fn analysis_benchmark_rejects_invalid_size_and_run_count() {
    for args in [["63", "1"], ["1000", "0"], ["abc", "1"]] {
        let output = run_analysis_benchmark(&args);

        assert_eq!(output.status.code(), Some(2), "args: {args:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("BASES must be an integer"),
            "args: {args:?}"
        );
    }
}

#[test]
fn config_uses_documented_defaults_and_ignores_cargo_bench_flag() {
    let config = parse_config(args(&["--bench"]))
        .unwrap()
        .expect("expected benchmark configuration");

    assert_eq!(config.input_size, DEFAULT_INPUT_SIZE);
    assert_eq!(
        config.sample_duration,
        Duration::from_millis(DEFAULT_SAMPLE_MS)
    );
}

#[test]
fn config_accepts_input_size_and_duration_overrides() {
    let config = parse_config(args(&[
        "--input-size",
        "4096",
        "--sample-ms",
        "75",
        "--bench",
    ]))
    .unwrap()
    .expect("expected benchmark configuration");

    assert_eq!(config.input_size, 4096);
    assert_eq!(config.sample_duration, Duration::from_millis(75));
}

#[test]
fn config_help_does_not_start_a_benchmark() {
    assert!(parse_config(args(&["--help"])).unwrap().is_none());
}

#[test]
fn config_usage_names_the_selected_benchmark() {
    let error = parse_config_for(args(&["--input-size", "0"]), "fasta_io")
        .err()
        .expect("zero input size should be rejected");

    assert!(error.contains("cargo bench --bench fasta_io"));
    assert!(error.contains("--input-size SIZE"));
}

#[test]
fn config_rejects_missing_zero_invalid_and_unknown_values() {
    for invalid_args in [
        &["--input-size"][..],
        &["--input-size", "0"][..],
        &["--input-size", "abc"][..],
        &["--sample-ms"][..],
        &["--sample-ms", "0"][..],
        &["--sample-ms", "abc"][..],
        &["--unknown"][..],
    ] {
        assert!(
            parse_config(args(invalid_args)).is_err(),
            "accepted invalid args: {invalid_args:?}"
        );
    }
}

#[test]
fn repeated_bytes_has_requested_length_and_cycles_the_input() {
    assert_eq!(repeated_bytes(b"ACG", 8), b"ACGACGAC");
    assert!(repeated_bytes(b"ACG", 0).is_empty());
}

#[test]
fn mixed_input_is_reproducible_and_contains_only_expected_symbols() {
    let first = mixed_sequence(NUCLEIC_ACID_SET, 50, 0x9e37_79b9, 10_000);
    let second = mixed_sequence(NUCLEIC_ACID_SET, 50, 0x9e37_79b9, 10_000);

    assert_eq!(first, second);
    assert_eq!(first.len(), 10_000);
    assert!(
        first
            .iter()
            .all(|base| { NUCLEIC_ACID_SET.contains(base) || INVALID_BASES.contains(base) })
    );
}

#[test]
fn mixed_input_tracks_requested_valid_percentages() {
    for (valid_percent, seed, min_percent, max_percent) in
        [(50, 0x9e37_79b9, 48, 52), (99, 0x243f_6a88, 98, 100)]
    {
        let sequence = mixed_sequence(NUCLEIC_ACID_SET, valid_percent, seed, 10_000);
        let valid_count = sequence
            .iter()
            .filter(|base| NUCLEIC_ACID_SET.contains(base))
            .count();
        let actual_percent = valid_count * 100 / sequence.len();

        assert!(
            (min_percent..=max_percent).contains(&actual_percent),
            "requested {valid_percent}% valid symbols, got {actual_percent}%"
        );
    }
}
