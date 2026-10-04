use std::{
    io::Write,
    process::{Command, Output, Stdio},
};

fn run_fasta_util(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_fasta-util"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to start fasta-util");

    child
        .stdin
        .as_mut()
        .expect("stdin was not piped")
        .write_all(input)
        .expect("failed to write stdin");
    drop(child.stdin.take());

    child.wait_with_output().expect("failed to collect output")
}

#[test]
fn len_counts_sequence_symbols_from_stdin() {
    let output = run_fasta_util(&["len"], b">record 1\nACGT\n \t\nNU-\n>record 2\n");

    assert!(output.status.success());
    assert_eq!(output.stdout, b"7\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn slice_writes_selected_sequence_from_stdin() {
    let output = run_fasta_util(
        &["slice", "--range", "2..=4", "--chars-per-line", "2"],
        b">record\nACGT\nNU-\n",
    );

    assert!(output.status.success());
    assert_eq!(output.stdout, b">record\nGT\nN\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn invalid_sequence_returns_an_error_without_panicking() {
    let output = run_fasta_util(&["len"], b"ACX\n");
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("invalid nucleic acid"));
    assert!(!stderr.contains("panicked"));
}
