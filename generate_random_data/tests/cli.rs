use std::process::{Command, Output};

fn run_generator(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_generate_random_data"))
        .args(args)
        .output()
        .expect("failed to run random-data generator")
}

#[test]
fn header_and_sequence_length_match_requested_size() {
    let requested_size = 123;
    let output = run_generator(&[&requested_size.to_string()]);

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
    assert!(sequence.lines().all(|line| (1..=50).contains(&line.len())));
    assert!(
        sequence
            .bytes()
            .filter(|byte| !byte.is_ascii_whitespace())
            .all(|base| b"ACGTNUKSYMWRBDHV".contains(&base))
    );
}

#[test]
fn seed_reproduces_the_same_fasta_output() {
    let first = run_generator(&["123", "--seed", "42"]);
    let second = run_generator(&["123", "--seed", "42"]);

    assert!(first.status.success());
    assert!(second.status.success());
    assert_eq!(first.stdout, second.stdout);
}

#[test]
fn custom_line_width_wraps_sequence_at_the_requested_width() {
    let output = run_generator(&["15", "--line-width", "7", "--seed", "42"]);

    assert!(
        output.status.success(),
        "generator failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("generator output was not UTF-8");
    let (header, sequence) = stdout
        .split_once('\n')
        .expect("generator output did not contain a FASTA header");
    let lines = sequence.lines().collect::<Vec<_>>();

    assert_eq!(header, ">TestData 15 random data");
    assert_eq!(lines.iter().map(|line| line.len()).sum::<usize>(), 15);
    assert_eq!(
        lines.iter().map(|line| line.len()).collect::<Vec<_>>(),
        [7, 7, 1]
    );
}

#[test]
fn rejects_a_zero_line_width() {
    let output = run_generator(&["15", "--line-width", "0"]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("line width must be greater than zero"));
}

#[test]
fn rejects_extra_positional_arguments() {
    let output = run_generator(&["123", "extra"]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("unexpected argument 'extra'"));
}
