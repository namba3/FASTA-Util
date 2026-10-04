use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::PathBuf,
    process::{Command, Output, Stdio},
    sync::atomic::{AtomicUsize, Ordering},
};

static NEXT_TEMP_FILE_ID: AtomicUsize = AtomicUsize::new(0);

struct TemporaryFile(PathBuf);

impl TemporaryFile {
    fn new(contents: &[u8]) -> Self {
        loop {
            let id = NEXT_TEMP_FILE_ID.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir()
                .join(format!("fasta-util-cli-{}-{id}.tmp", std::process::id()));
            match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(mut file) => {
                    file.write_all(contents).unwrap();
                    return Self(path);
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("failed to create temporary file: {error}"),
            }
        }
    }

    fn path(&self) -> &str {
        self.0.to_str().expect("temporary path is not valid UTF-8")
    }

    fn read(&self) -> Vec<u8> {
        let mut contents = Vec::new();
        File::open(&self.0)
            .unwrap()
            .read_to_end(&mut contents)
            .unwrap();
        contents
    }
}

impl Drop for TemporaryFile {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

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
fn len_counts_sequence_symbols_from_a_multi_record_file() {
    let input = TemporaryFile::new(b">first\r\nACGT\r\n\r\n>second\r\nNU-\r\n");
    let output = run_fasta_util(&["len", "--input", input.path()], b"");

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
fn slice_reads_a_crlf_file_and_writes_normalized_fasta_to_a_file() {
    let input = TemporaryFile::new(b">first\r\nACGT\r\n>second\r\nNU-\r\n");
    let output_file = TemporaryFile::new(b"old contents");
    let output = run_fasta_util(
        &[
            "slice",
            "--input",
            input.path(),
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    assert_eq!(output_file.read(), b">first\nACGT\n>second\nNU-");
}

#[test]
fn slice_keeps_existing_output_when_input_validation_fails_after_partial_output() {
    let input = TemporaryFile::new(b">record\nACGT\nACX\n");
    let output_file = TemporaryFile::new(b"previous output");
    let output = run_fasta_util(
        &[
            "slice",
            "--input",
            input.path(),
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid nucleic acid"));
    assert_eq!(output_file.read(), b"previous output");
}

#[test]
fn slice_does_not_create_output_when_input_validation_fails() {
    let input = TemporaryFile::new(b">record\nACX\n");
    let output_file = TemporaryFile::new(b"remove me");
    fs::remove_file(&output_file.0).unwrap();

    let output = run_fasta_util(
        &[
            "slice",
            "--input",
            input.path(),
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid nucleic acid"));
    assert!(!output_file.0.exists());
}

#[test]
fn slice_rejects_same_input_and_output_file_without_changing_it() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let output = run_fasta_util(
        &["slice", "--input", input.path(), "--output", input.path()],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("same file"));
    assert_eq!(input.read(), b">record\nACGT\n");
}

#[test]
fn invalid_sequence_returns_an_error_without_panicking() {
    let output = run_fasta_util(&["len"], b"ACX\n");
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("invalid nucleic acid"));
    assert!(!stderr.contains("panicked"));
}
