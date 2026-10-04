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

fn indexed_fasta_fixture() -> (TemporaryFile, TemporaryFile) {
    let fasta = b">empty record\r\n>first description\r\nACGT\r\nACGT\r\n>second record\r\nacgt\r\n>third\r\nTT";
    let first_header = b">empty record\r\n";
    let second_header = b">first description\r\n";
    let third_header = b">second record\r\n";
    let fourth_header = b">third\r\n";
    let empty_offset = first_header.len();
    let first_offset = empty_offset + second_header.len();
    let second_offset = first_offset + 2 * 6 + third_header.len();
    let third_offset = second_offset + 6 + fourth_header.len();
    let fai = format!(
        "empty\t0\t{empty_offset}\t0\t0\nfirst\t8\t{first_offset}\t4\t6\nsecond\t4\t{second_offset}\t4\t6\nthird\t2\t{third_offset}\t2\t2\n"
    );
    (
        TemporaryFile::new(fasta),
        TemporaryFile::new(fai.as_bytes()),
    )
}

fn run_fasta_util(args: &[&str], input: &[u8]) -> Output {
    run_fasta_util_with(input, |command| {
        command.args(args);
    })
}

fn run_fasta_util_with(input: &[u8], configure: impl FnOnce(&mut Command)) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fasta-util"));
    configure(&mut command);
    let mut child = command
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

#[cfg(unix)]
fn temporary_non_utf8_file(contents: &[u8]) -> TemporaryFile {
    use std::os::unix::ffi::OsStringExt;

    loop {
        let id = NEXT_TEMP_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let mut name = format!("fasta-util-cli-{}-{id}-", std::process::id()).into_bytes();
        name.push(0xff);
        name.extend_from_slice(b".tmp");
        let path = std::env::temp_dir().join(std::ffi::OsString::from_vec(name));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(mut file) => {
                file.write_all(contents).unwrap();
                return TemporaryFile(path);
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => panic!("failed to create non-UTF-8 temporary file: {error}"),
        }
    }
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
fn len_counts_lowercase_soft_masked_bases_from_a_file() {
    let input = TemporaryFile::new(b">record\r\nacgtn\r\nu-\r\n");
    let output = run_fasta_util(&["len", "--input", input.path()], b"");

    assert!(output.status.success());
    assert_eq!(output.stdout, b"7\n");
    assert!(output.stderr.is_empty());
}

#[cfg(unix)]
#[test]
fn len_accepts_a_non_utf8_input_path() {
    let input = temporary_non_utf8_file(b">record\nACGT\n");
    let output = run_fasta_util_with(b"", |command| {
        command.arg("len").arg("--input").arg(&input.0);
    });

    assert!(output.status.success());
    assert_eq!(output.stdout, b"4\n");
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
fn indexed_slice_matches_streaming_slice_for_wrapped_multi_record_fasta() {
    let (input, index) = indexed_fasta_fixture();
    let cases = [
        ("..", "2"),
        ("0..0", "4"),
        ("3..=9", "4"),
        ("8..=8", "3"),
        ("10..14", "5"),
        ("14..", "4"),
        ("99..", "4"),
    ];

    for (range, chars_per_line) in cases {
        let normal = run_fasta_util(
            &[
                "slice",
                "--input",
                input.path(),
                "--range",
                range,
                "--chars-per-line",
                chars_per_line,
            ],
            b"",
        );
        let indexed = run_fasta_util(
            &[
                "slice",
                "--input",
                input.path(),
                "--fai-index",
                index.path(),
                "--range",
                range,
                "--chars-per-line",
                chars_per_line,
            ],
            b"",
        );

        assert!(normal.status.success(), "normal range {range}: {normal:?}");
        assert!(
            indexed.status.success(),
            "indexed range {range}: {indexed:?}"
        );
        assert_eq!(indexed.stdout, normal.stdout, "range {range}");
        assert!(
            indexed.stderr.is_empty(),
            "range {range}: {:?}",
            indexed.stderr
        );
    }
}

#[test]
fn indexed_slice_matches_streaming_slice_across_read_buffer_chunks() {
    let mut fasta = b">long\n".to_vec();
    let sequence = (0..150_123)
        .map(|index| b"ACGTacgt"[index % 8])
        .collect::<Vec<_>>();
    for line in sequence.chunks(60) {
        fasta.extend_from_slice(line);
        if line.len() == 60 {
            fasta.push(b'\n');
        }
    }
    let input = TemporaryFile::new(&fasta);
    let index = TemporaryFile::new(b"long\t150123\t6\t60\t61\n");
    let cases = [("57..70031", "73"), ("65510..=140007", "61")];

    for (range, chars_per_line) in cases {
        let normal = run_fasta_util(
            &[
                "slice",
                "--input",
                input.path(),
                "--range",
                range,
                "--chars-per-line",
                chars_per_line,
            ],
            b"",
        );
        let indexed = run_fasta_util(
            &[
                "slice",
                "--input",
                input.path(),
                "--fai-index",
                index.path(),
                "--range",
                range,
                "--chars-per-line",
                chars_per_line,
            ],
            b"",
        );

        assert!(normal.status.success(), "normal range {range}: {normal:?}");
        assert!(
            indexed.status.success(),
            "indexed range {range}: {indexed:?}"
        );
        assert_eq!(indexed.stdout, normal.stdout, "range {range}");
    }
}

#[test]
fn indexed_slice_validates_only_selected_sequence_bases() {
    let input = TemporaryFile::new(b">record description\nACXT\n");
    let index = TemporaryFile::new(b"record\t4\t20\t4\t5\n");
    let output = run_fasta_util(
        &[
            "slice",
            "--input",
            input.path(),
            "--fai-index",
            index.path(),
            "--range",
            "0..2",
        ],
        b"",
    );

    assert!(output.status.success(), "{:?}", output.stderr);
    assert_eq!(output.stdout, b">record description\nAC\n");
}

#[test]
fn indexed_slice_rejects_mismatched_index_without_replacing_output() {
    let (input, _) = indexed_fasta_fixture();
    let index = TemporaryFile::new(b"wrong-name\t8\t35\t4\t6\n");
    let output_file = TemporaryFile::new(b"keep this output");
    let output = run_fasta_util(
        &[
            "slice",
            "--input",
            input.path(),
            "--fai-index",
            index.path(),
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("FAI record name"));
    assert_eq!(output_file.read(), b"keep this output");
}

#[test]
fn indexed_slice_rejects_index_as_output_without_changing_it() {
    let (input, index) = indexed_fasta_fixture();
    let original_index = index.read();
    let output = run_fasta_util(
        &[
            "slice",
            "--input",
            input.path(),
            "--fai-index",
            index.path(),
            "--output",
            index.path(),
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("same file"));
    assert_eq!(index.read(), original_index);
}

#[test]
fn indexed_slice_rejects_malformed_and_out_of_bounds_indexes() {
    let (input, _) = indexed_fasta_fixture();
    for contents in [
        &b"record\tlength\t0\t4\t4\n"[..],
        &b"record\t4\t999999\t4\t4\n"[..],
        &b"record\t4\t12\t0\t0\n"[..],
    ] {
        let index = TemporaryFile::new(contents);
        let output = run_fasta_util(
            &[
                "slice",
                "--input",
                input.path(),
                "--fai-index",
                index.path(),
            ],
            b"",
        );
        assert!(!output.status.success(), "index {contents:?} was accepted");
        assert!(!String::from_utf8_lossy(&output.stderr).contains("panicked"));
    }
}

#[test]
fn indexed_slice_requires_a_file_input() {
    let index = TemporaryFile::new(b"record\t4\t7\t4\t4\n");
    let output = run_fasta_util(&["slice", "--fai-index", index.path()], b">record\nACGT\n");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--fai-index requires --input"));
}

#[cfg(unix)]
#[test]
fn slice_accepts_a_non_utf8_output_path() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let output_file = temporary_non_utf8_file(b"old contents");
    let output = run_fasta_util_with(b"", |command| {
        command
            .arg("slice")
            .arg("--input")
            .arg(&input.0)
            .arg("--output")
            .arg(&output_file.0);
    });

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    assert_eq!(output_file.read(), b">record\nACGT");
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
