use std::process::Command;

#[test]
fn header_and_sequence_length_match_requested_size() {
    let requested_size = 123;
    let output = Command::new(env!("CARGO_BIN_EXE_generate_random_data"))
        .arg(requested_size.to_string())
        .output()
        .expect("failed to run random-data generator");

    assert!(
        output.status.success(),
        "generator failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).expect("generator output was not UTF-8");
    let (header, sequence) = stdout
        .split_once('\n')
        .expect("generator output did not contain a FASTA header");
    let sequence_length = sequence.lines().map(str::len).sum::<usize>();

    assert_eq!(header, format!(">TestData {requested_size} random data"));
    assert_eq!(sequence_length, requested_size);
}
