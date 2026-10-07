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

fn indexed_fasta_with_layout(
    records: &[(&str, &[u8], usize)],
    line_ending: &[u8],
) -> (TemporaryFile, TemporaryFile) {
    let mut fasta = Vec::new();
    let mut fai = String::new();

    for (record_index, (name, sequence, requested_line_bases)) in records.iter().enumerate() {
        fasta.extend_from_slice(b">");
        fasta.extend_from_slice(name.as_bytes());
        fasta.extend_from_slice(b" generated record");
        fasta.extend_from_slice(line_ending);

        let sequence_offset = fasta.len();
        let line_bases = if sequence.is_empty() {
            0
        } else {
            (*requested_line_bases).min(sequence.len())
        };
        let line_width = if sequence.is_empty() {
            0
        } else {
            line_bases + line_ending.len()
        };
        if line_bases > 0 {
            let chunks = sequence.chunks(line_bases).collect::<Vec<_>>();
            for (chunk_index, chunk) in chunks.iter().enumerate() {
                fasta.extend_from_slice(chunk);
                let is_final_chunk_of_final_record =
                    record_index + 1 == records.len() && chunk_index + 1 == chunks.len();
                if !is_final_chunk_of_final_record {
                    fasta.extend_from_slice(line_ending);
                }
            }
        }

        fai.push_str(&format!(
            "{name}\t{}\t{sequence_offset}\t{line_bases}\t{line_width}\n",
            sequence.len()
        ));
    }

    (
        TemporaryFile::new(&fasta),
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

fn run_validate(contents: &[u8], options: &[&str]) -> Output {
    let input = TemporaryFile::new(contents);
    let mut args = vec!["validate", input.path()];
    args.extend_from_slice(options);
    run_fasta_util(&args, b"")
}

#[test]
fn validate_reports_record_count_and_dna_type() {
    let output = run_validate(b">first description\nACGT\n>second\nTTAA\n", &[]);

    assert!(output.status.success());
    assert_eq!(output.stdout, b"OK: 2 records\ntype: DNA\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn validate_uses_the_first_nonwhitespace_header_identifier() {
    let output = run_validate(b">  seq description\nACGT\n", &[]);

    assert!(output.status.success());
    assert_eq!(output.stdout, b"OK: 1 records\ntype: DNA\n");
}

#[test]
fn validate_reports_rna_and_ambiguous_nucleotide_types() {
    let rna = run_validate(b">rna\nACGU\n", &[]);
    let ambiguous = run_validate(b">ambiguous\nACGN-\n", &[]);

    assert!(rna.status.success());
    assert_eq!(rna.stdout, b"OK: 1 records\ntype: RNA\n");
    assert!(ambiguous.status.success());
    assert_eq!(
        ambiguous.stdout,
        b"OK: 1 records\ntype: DNA/RNA ambiguous\n"
    );
}

#[test]
fn validate_rejects_mixed_dna_and_rna_symbols() {
    let output = run_validate(b">dna\nACGT\n>rna\nACGU\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("DNA and RNA symbols (`T` and `U`) are mixed"));
    assert!(stderr.contains(":4:4\n"));
}

#[test]
fn validate_accepts_protein_alphabet_and_gaps_when_selected() {
    let output = run_validate(
        b">protein description\nACDEFGHIKLMNPQRSTVWY\nBJOUXZ*-\n",
        &["--sequence-type", "protein"],
    );

    assert!(output.status.success());
    assert_eq!(output.stdout, b"OK: 1 records\ntype: protein\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn validate_rejects_sequence_before_a_header() {
    let output = run_validate(b"ACGT\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("sequence data appears before the first `>` record"));
    assert!(stderr.contains(":1:1\n"));
    assert!(stderr.contains("error: no FASTA records found"));
}

#[test]
fn validate_reports_empty_record_id_and_empty_sequence() {
    let output = run_validate(b">\n>empty\n>nonempty\nA\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("record identifier is empty"));
    assert!(stderr.contains("record has an empty sequence"));
    assert!(stderr.contains(":1:2\n"));
}

#[test]
fn validate_reports_duplicate_record_ids_with_first_occurrence() {
    let output = run_validate(b">same first\nA\n>same second\nC\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("duplicate record identifier `same`"));
    assert!(stderr.contains("identifier `same` first appeared on line 1"));
    assert!(stderr.contains(":3:2\n"));
}

#[test]
fn validate_reports_invalid_symbol_at_its_line_and_column() {
    let output = run_validate(b">seq\nACGTZACGT\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("error: invalid nucleotide 'Z'"));
    assert!(stderr.contains(":2:5\n"));
    assert!(stderr.contains("2 | ACGTZACGT\n  |     ^"));
}

#[test]
fn validate_reports_multiple_errors_in_one_pass_without_stdout() {
    let output = run_validate(b">duplicate\nACZ\n>duplicate\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(stderr.contains("invalid nucleotide 'Z'"));
    assert!(stderr.contains("duplicate record identifier `duplicate`"));
    assert!(stderr.contains("record has an empty sequence"));
}

#[test]
fn validate_reports_whitespace_inside_sequence_lines() {
    let output = run_validate(b">seq\nAC GT\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("whitespace is not allowed in sequence lines"));
    assert!(stderr.contains(":2:3\n"));
}

#[test]
fn validate_rejects_mixed_lf_and_crlf_endings() {
    let output = run_validate(b">seq\nACGT\nTGCA\r\n", &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("mixed LF and CRLF line endings"));
    assert!(stderr.contains(":3:5\n"));
}

#[test]
fn validate_allows_different_line_endings_between_records() {
    let output = run_validate(b">first\nACGT\n>second\r\nTGCA\r\n", &[]);

    assert!(output.status.success());
    assert_eq!(output.stdout, b"OK: 2 records\ntype: DNA\n");
}

#[test]
fn validate_checks_fai_line_widths_but_allows_a_short_final_line() {
    let valid = run_validate(b">seq\nACGT\nAC\n", &[]);
    let invalid = run_validate(b">seq\nACGT\nAC\nGT\n", &[]);
    let wider = run_validate(b">seq\nAC\nGTA\n", &[]);
    let stderr = String::from_utf8_lossy(&invalid.stderr);
    let wider_stderr = String::from_utf8_lossy(&wider.stderr);

    assert!(valid.status.success());
    assert_eq!(valid.stdout, b"OK: 1 records\ntype: DNA\n");
    assert!(!invalid.status.success());
    assert!(stderr.contains("non-final sequence line is shorter than the `.fai` line width"));
    assert!(stderr.contains(":3:3\n"));
    assert!(!wider.status.success());
    assert!(wider_stderr.contains("sequence line is wider than the first line"));
    assert!(wider_stderr.contains(":3:3\n"));
}

#[test]
fn index_writes_standard_fai_rows_for_wrapped_multirecord_input() {
    let input = TemporaryFile::new(b">chr1 description\nACGT\nTGCA\n>chr2\r\nAACC\r\nGG");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let output = run_fasta_util(&["index", input.path()], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("Indexed 2 records:"));
    assert_eq!(
        fs::read(&index_path).unwrap(),
        b"chr1\t8\t18\t4\t5\nchr2\t6\t35\t4\t6\n"
    );
    let slice = run_fasta_util(
        &[
            "get",
            input.path(),
            "--fai-index",
            index_path.to_str().unwrap(),
            "3-6",
        ],
        b"",
    );
    assert!(
        slice.status.success(),
        "{}",
        String::from_utf8_lossy(&slice.stderr)
    );
    assert_eq!(slice.stdout, b">chr1 description\nGTTG\n");
    fs::remove_file(index_path).unwrap();
}

#[test]
fn index_failure_preserves_an_existing_index_file() {
    let input = TemporaryFile::new(b">record\nACGT\nAC\nGG\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    fs::write(&index_path, b"existing index\n").unwrap();

    let output = run_fasta_util(&["index", input.path()], b"");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("non-final sequence line is shorter"));
    assert_eq!(fs::read(&index_path).unwrap(), b"existing index\n");
    fs::remove_file(index_path).unwrap();
}

#[test]
fn index_skips_leading_header_whitespace_and_supports_get() {
    let input = TemporaryFile::new(b">  seq description\nACGT\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let index = run_fasta_util(&["index", input.path()], b"");

    assert!(
        index.status.success(),
        "{}",
        String::from_utf8_lossy(&index.stderr)
    );
    assert_eq!(fs::read(&index_path).unwrap(), b"seq\t4\t19\t4\t5\n");
    let slice = run_fasta_util(
        &[
            "get",
            input.path(),
            "--fai-index",
            index_path.to_str().unwrap(),
            "2-3",
        ],
        b"",
    );
    assert!(
        slice.status.success(),
        "{}",
        String::from_utf8_lossy(&slice.stderr)
    );
    assert_eq!(slice.stdout, b">  seq description\nCG\n");
    fs::remove_file(index_path).unwrap();
}

#[test]
fn stats_reports_aggregate_lengths_n50_and_nucleotide_percentages() {
    let input = TemporaryFile::new(b">a\nACGTNN\n>b\nGGT\n");
    let output = run_fasta_util(&["stats", input.path()], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"sequences    2\ntotal_len    9\nmin_len      3\nmax_len      6\nmean_len     5\nN50          6\nGC           44.44%\nN            22.22%\ntype         DNA\n"
    );
}

#[test]
fn composition_reports_nucleotide_symbol_and_gc_percentages() {
    let input = TemporaryFile::new(b">a\nAaCcGgTtNnR-\n");
    let output = run_fasta_util(&["composition", input.path()], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"A\t16.67%\nC\t16.67%\nG\t16.67%\nT\t16.67%\nN\t16.67%\nR\t8.33%\n-\t8.33%\nGC\t33.33%\n"
    );
}

#[test]
fn composition_reports_protein_symbols_and_omits_gc() {
    let input = TemporaryFile::new(b">protein\nACDEX*\n");
    let output = run_fasta_util(
        &["composition", input.path(), "--sequence-type", "protein"],
        b"",
    );
    let stdout = String::from_utf8_lossy(&output.stdout);

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(stdout.starts_with("A\t16.67%\nC\t16.67%\nD\t16.67%\nE\t16.67%\n"));
    assert!(stdout.contains("Y\t0.00%\n"));
    assert!(stdout.contains("X\t16.67%\n"));
    assert!(stdout.contains("*\t16.67%\n"));
    assert!(!stdout.lines().any(|line| line.starts_with("GC\t")));
}

#[test]
fn composition_reads_rna_from_standard_input_and_rejects_invalid_symbols() {
    let output = run_fasta_util(&["composition", "-"], b">rna\r\nACGU\r\n");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"A\t25.00%\nC\t25.00%\nG\t25.00%\nU\t25.00%\nN\t0.00%\nGC\t50.00%\n"
    );

    let invalid = run_fasta_util(&["composition"], b">record\nAC?T\n");
    assert!(!invalid.status.success());
    assert!(
        String::from_utf8_lossy(&invalid.stderr).contains("line 2: invalid nucleotide symbol '?'")
    );
}

#[test]
fn stats_reads_fasta_from_standard_input_when_input_is_omitted_or_dash() {
    let contents = b">record\nACGT\n";
    for args in [&["stats"][..], &["stats", "-"][..]] {
        let output = run_fasta_util(args, contents);

        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("total_len    4"));
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn stats_each_reports_per_record_gc_and_n_percentages() {
    let input = TemporaryFile::new(b">a\nACGTNN\n>b\nGGT\n");
    let output = run_fasta_util(&["stats", input.path(), "--each"], b"");
    let stdout = String::from_utf8_lossy(&output.stdout);
    let rows = stdout
        .lines()
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .collect::<Vec<_>>();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(rows[0], ["id", "length", "gc", "n"]);
    assert_eq!(rows[1], ["a", "6", "33.33%", "33.33%"]);
    assert_eq!(rows[2], ["b", "3", "66.67%", "0.00%"]);
}

#[test]
fn stats_json_emits_a_summary_object_and_each_array() {
    let input = TemporaryFile::new(b">a\nACGTNN\n>b\nGGT\n");
    let summary = run_fasta_util(&["stats", input.path(), "--format", "json"], b"");
    let each = run_fasta_util(&["stats", input.path(), "--each", "--format", "json"], b"");

    assert!(summary.status.success());
    assert_eq!(
        summary.stdout,
        b"{\n  \"sequences\": 2,\n  \"total_len\": 9,\n  \"min_len\": 3,\n  \"max_len\": 6,\n  \"mean_len\": 4.5,\n  \"n50\": 6,\n  \"gc_percent\": 44.444444,\n  \"n_percent\": 22.222222,\n  \"type\": \"DNA\"\n}\n"
    );
    assert!(each.status.success());
    assert_eq!(
        each.stdout,
        b"[\n  {\"id\": \"a\", \"length\": 6, \"gc_percent\": 33.333333, \"n_percent\": 33.333333},\n  {\"id\": \"b\", \"length\": 3, \"gc_percent\": 66.666667, \"n_percent\": 0.000000}\n]\n"
    );
}

#[test]
fn stats_detects_protein_and_uses_null_gc_fields_in_json() {
    let input = TemporaryFile::new(b">protein\nACDE\n");
    let text = run_fasta_util(&["stats", input.path()], b"");
    let json = run_fasta_util(&["stats", input.path(), "--format", "json"], b"");

    assert!(text.status.success());
    assert!(
        String::from_utf8_lossy(&text.stdout)
            .contains("GC           n/a\nN            n/a\ntype         Protein\n")
    );
    assert!(json.status.success());
    assert!(String::from_utf8_lossy(&json.stdout).contains("\"gc_percent\": null"));
    assert!(String::from_utf8_lossy(&json.stdout).contains("\"n_percent\": null"));
}

#[test]
fn stats_rejects_whitespace_inside_sequence_lines() {
    let input = TemporaryFile::new(b">a\nAC GT\n");
    let output = run_fasta_util(&["stats", input.path()], b"");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("line 2: whitespace in sequence at column 3")
    );
}

#[test]
fn stats_rejects_sequence_data_before_first_record() {
    let input = TemporaryFile::new(b"ACGT\n>a\nACGT\n");
    let output = run_fasta_util(&["stats", input.path()], b"");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("line 1: sequence data appears before the first `>` record")
    );
}

#[test]
fn filter_combines_inclusive_length_bounds_and_preserves_record_order() {
    let input = TemporaryFile::new(
        b">short\nAC\n>first-match description\nACGT\n>long\nACGTACG\n>second-match\nACG\n",
    );
    let output = run_fasta_util(
        &["filter", input.path(), "--min-len", "3", "--max-len", "4"],
        b"",
    );

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b">first-match description\nACGT\n>second-match\nACG\n"
    );
}

#[test]
fn filter_reads_fasta_from_standard_input_when_input_is_omitted_or_dash() {
    let contents = b">short\nACG\n>keep\nACGT\n";
    for args in [
        &["filter", "--min-len", "4"][..],
        &["filter", "-", "--min-len", "4"][..],
    ] {
        let output = run_fasta_util(args, contents);

        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, b">keep\nACGT\n");
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn filter_combines_gc_and_n_fraction_bounds_case_insensitively() {
    let input =
        TemporaryFile::new(b">low-gc\nATAT\n>match\naCgTN\n>high-gc\nGGGCCN\n>high-n\nNNNN\n");
    let output = run_fasta_util(
        &[
            "filter",
            input.path(),
            "--min-gc",
            "0.4",
            "--max-gc",
            "0.4",
            "--max-n",
            "0.2",
        ],
        b"",
    );

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">match\naCgTN\n");
}

#[test]
fn filter_supports_protein_length_and_preserves_crlf_bytes() {
    let input = TemporaryFile::new(b">short protein\r\nACDE\r\n>long protein\r\nACDEFGHIK\r\n");
    let output = run_fasta_util(&["filter", input.path(), "--min-len", "5"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">long protein\r\nACDEFGHIK\r\n");
}

#[test]
fn filter_writes_selected_records_to_a_file() {
    let input = TemporaryFile::new(b">short\nAC\n>long\nACGT\n");
    let output_file = TemporaryFile::new(b"previous output");
    let output = run_fasta_util(
        &[
            "filter",
            input.path(),
            "--min-len",
            "3",
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(output_file.read(), b">long\nACGT\n");
}

#[test]
fn filter_and_grep_keep_selection_order_across_bitmap_word_boundaries() {
    let selected_indices = [0, 63, 64, 69];
    let mut contents = Vec::new();
    let mut expected = Vec::new();
    for index in 0..70 {
        let selected = selected_indices.contains(&index);
        let header = if selected {
            format!(">keep-{index}\n")
        } else {
            format!(">skip-{index}\n")
        };
        contents.extend_from_slice(header.as_bytes());
        if selected {
            contents.extend_from_slice(b"AC\n");
            expected.extend_from_slice(header.as_bytes());
            expected.extend_from_slice(b"AC\n");
        } else {
            contents.extend_from_slice(b"A\n");
        }
    }
    let input = TemporaryFile::new(&contents);

    let filtered = run_fasta_util(&["filter", input.path(), "--min-len", "2"], b"");
    let grepped = run_fasta_util(&["grep", input.path(), "keep-"], b"");

    assert!(
        filtered.status.success(),
        "{}",
        String::from_utf8_lossy(&filtered.stderr)
    );
    assert!(
        grepped.status.success(),
        "{}",
        String::from_utf8_lossy(&grepped.stderr)
    );
    assert_eq!(filtered.stdout, expected);
    let grep_expected = selected_indices
        .iter()
        .flat_map(|index| [format!(">keep-{index}\n").into_bytes(), b"AC\n".to_vec()])
        .flatten()
        .collect::<Vec<_>>();
    assert_eq!(grepped.stdout, grep_expected);
}

#[test]
fn filter_and_grep_write_nothing_when_no_records_are_selected() {
    let input = TemporaryFile::new(b">one\nAC\n>two\nGT\n");

    let filtered = run_fasta_util(&["filter", input.path(), "--min-len", "3"], b"");
    let grepped = run_fasta_util(&["grep", input.path(), "missing"], b"");

    assert!(
        filtered.status.success(),
        "{}",
        String::from_utf8_lossy(&filtered.stderr)
    );
    assert!(
        grepped.status.success(),
        "{}",
        String::from_utf8_lossy(&grepped.stderr)
    );
    assert!(filtered.stdout.is_empty());
    assert!(grepped.stdout.is_empty());
}

#[test]
fn filter_rejects_missing_or_contradictory_conditions() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let no_conditions = run_fasta_util(&["filter", input.path()], b"");
    let reversed_length = run_fasta_util(
        &["filter", input.path(), "--min-len", "5", "--max-len", "4"],
        b"",
    );
    let reversed_gc = run_fasta_util(
        &["filter", input.path(), "--min-gc", "0.8", "--max-gc", "0.2"],
        b"",
    );
    let protein_gc = run_fasta_util(
        &[
            "filter",
            input.path(),
            "--sequence-type",
            "protein",
            "--max-gc",
            "0.5",
        ],
        b"",
    );
    let auto_protein_gc_input = TemporaryFile::new(b">protein\nACDE\n");
    let auto_protein_gc = run_fasta_util(
        &["filter", auto_protein_gc_input.path(), "--max-gc", "0.5"],
        b"",
    );

    for (output, message) in [
        (no_conditions, "provide at least one filter option"),
        (
            reversed_length,
            "--min-len cannot be greater than --max-len",
        ),
        (reversed_gc, "--min-gc cannot be greater than --max-gc"),
        (
            protein_gc,
            "GC and N filters are only available for nucleotide sequences",
        ),
        (
            auto_protein_gc,
            "GC and N filters are only available for nucleotide sequences",
        ),
    ] {
        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "expected {message:?}, got {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn filter_rejects_out_of_range_fractions() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    for fraction in ["-0.1", "1.1", "NaN", "inf"] {
        let option = format!("--max-n={fraction}");
        let output = run_fasta_util(&["filter", input.path(), &option], b"");

        assert!(
            !output.status.success(),
            "fraction {fraction} unexpectedly passed"
        );
        assert!(String::from_utf8_lossy(&output.stderr).contains("fraction"));
    }
}

#[test]
fn filter_validation_failure_preserves_existing_output() {
    let input = TemporaryFile::new(b">valid\nACGT\n>invalid\nAC?Z\n");
    let output_file = TemporaryFile::new(b"keep this output");
    let output = run_fasta_util(
        &[
            "filter",
            input.path(),
            "--min-len",
            "2",
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid nucleotide symbol '?'"));
    assert_eq!(output_file.read(), b"keep this output");
}

#[test]
fn revcomp_preserves_headers_and_transforms_dna_sequences() {
    let input = TemporaryFile::new(b">seq description\nACGTTGCA\n");
    let output = run_fasta_util(&["revcomp", input.path()], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">seq description\nTGCAACGT\n");
}

#[test]
fn revcomp_reads_fasta_from_standard_input_when_input_is_omitted_or_dash() {
    let contents = b">record\nACGTTGCA\n";
    for args in [&["revcomp"][..], &["revcomp", "-"][..]] {
        let output = run_fasta_util(args, contents);

        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, b">record\nTGCAACGT\n");
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn filter_revcomp_stats_compose_through_standard_streams() {
    let input = TemporaryFile::new(b">short\nACG\n>keep\nACGTTGCA\n");
    let mut filter = Command::new(env!("CARGO_BIN_EXE_fasta-util"))
        .args(["filter", "--min-len", "5", input.path()])
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to start filter");

    let filter_stdout = filter.stdout.take().expect("filter stdout was not piped");
    let mut revcomp = Command::new(env!("CARGO_BIN_EXE_fasta-util"))
        .arg("revcomp")
        .stdin(Stdio::from(filter_stdout))
        .stdout(Stdio::piped())
        .spawn()
        .expect("failed to start revcomp");

    let revcomp_stdout = revcomp.stdout.take().expect("revcomp stdout was not piped");
    let stats = Command::new(env!("CARGO_BIN_EXE_fasta-util"))
        .arg("stats")
        .stdin(Stdio::from(revcomp_stdout))
        .output()
        .expect("failed to start stats");

    assert!(filter.wait().unwrap().success());
    assert!(revcomp.wait().unwrap().success());
    assert!(
        stats.status.success(),
        "{}",
        String::from_utf8_lossy(&stats.stderr)
    );
    assert!(String::from_utf8_lossy(&stats.stdout).contains("sequences    1"));
    assert!(String::from_utf8_lossy(&stats.stdout).contains("total_len    8"));
    assert!(stats.stderr.is_empty());
}

#[test]
fn remaining_fasta_commands_accept_standard_input_with_dash_or_omitted_path() {
    let fasta = b">record motif\nACGTTGCA\n";

    let formatted = run_fasta_util(&["format"], fasta);
    assert!(formatted.status.success());
    assert_eq!(formatted.stdout, fasta);

    let got = run_fasta_util(&["get", "-", "record"], fasta);
    assert!(got.status.success());
    assert_eq!(got.stdout, fasta);

    let grep = run_fasta_util(&["grep", "-", "motif"], fasta);
    assert!(grep.status.success());
    assert_eq!(grep.stdout, fasta);

    let located = run_fasta_util(&["locate", "-", "ACG"], fasta);
    assert!(located.status.success());
    assert_eq!(located.stdout, b"record\t1\t3\t+\nrecord\t2\t4\t-\n");

    let validated = run_fasta_util(&["validate", "-"], fasta);
    assert!(validated.status.success());
    assert!(String::from_utf8_lossy(&validated.stdout).contains("OK: 1 records"));

    let length = run_fasta_util(&["len", "-i", "-"], fasta);
    assert!(length.status.success());
    assert_eq!(length.stdout, b"8\n");
}

#[test]
fn validate_streams_multirecord_stdin_with_omitted_path_or_dash() {
    let fasta = b">first\r\nACGT\r\n>second\nAACC\n";
    for args in [&["validate"][..], &["validate", "-"][..]] {
        let output = run_fasta_util(args, fasta);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, b"OK: 2 records\ntype: DNA\n");
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn get_reads_record_regions_and_global_ranges_from_standard_input() {
    let fasta = b">first description\nACGT\n>second description\nTGCA\n";

    let record_region = run_fasta_util(&["get", "-", "first:2-4"], fasta);
    assert!(record_region.status.success());
    assert_eq!(record_region.stdout, b">first:2-4\nCGT\n");

    let global_range = run_fasta_util(&["get", "-", "3-6"], fasta);
    assert!(global_range.status.success());
    assert_eq!(
        global_range.stdout,
        b">first description\nGT\n>second description\nTG\n"
    );
}

#[test]
fn stdin_inputs_support_file_outputs_and_validation_diagnostics_name_stdin() {
    let output_file = TemporaryFile::new(b"old contents");
    let formatted = run_fasta_util(
        &["format", "--uppercase", "--output", output_file.path()],
        b">record\r\nacgt\r\n",
    );
    assert!(formatted.status.success());
    assert!(formatted.stdout.is_empty());
    assert_eq!(output_file.read(), b">record\nACGT\n");

    let valid = run_fasta_util(&["validate"], b">record\nACGT\n");
    assert!(valid.status.success());
    assert!(valid.stderr.is_empty());

    let invalid = run_fasta_util(&["validate", "-"], b">record\nACZ\n");
    assert!(!invalid.status.success());
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("stdin:2:3"));
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("invalid nucleotide 'Z'"));
}

#[test]
fn stdin_filter_failure_does_not_write_partial_fasta_to_stdout() {
    let output = run_fasta_util(
        &["filter", "--min-len", "1"],
        b">valid\nACGT\n>invalid\nAC?\n",
    );

    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid nucleotide symbol"));
}

#[test]
fn revcomp_complements_iupac_symbols_and_preserves_case() {
    let input = TemporaryFile::new(b">upper\nACGTRYKMBVDHSWN-\n>lower\nacgtrykmbvdhswn-\n");
    let output = run_fasta_util(&["revcomp", input.path(), "--chars-per-line", "8"], b"");

    assert!(output.status.success());
    assert_eq!(
        output.stdout,
        b">upper\n-NWSDHBV\nKMRYACGT\n>lower\n-nwsdhbv\nkmryacgt\n"
    );
}

#[test]
fn revcomp_preserves_rna_alphabet_and_normalizes_line_endings() {
    let input = TemporaryFile::new(b">rna\r\nAA\r\nCGU\r\n>dna\nAACGT\n");
    let output = run_fasta_util(&["revcomp", input.path()], b"");

    assert!(output.status.success());
    assert_eq!(output.stdout, b">rna\nACGUU\n>dna\nACGTT\n");
}

#[test]
fn revcomp_handles_sequences_larger_than_its_reverse_read_buffer() {
    let sequence = (0..65_539)
        .map(|index| b"ACGTN-"[index % 6])
        .collect::<Vec<_>>();
    let mut fasta = b">long\n".to_vec();
    fasta.extend_from_slice(&sequence);
    fasta.push(b'\n');
    let input = TemporaryFile::new(&fasta);
    let output = run_fasta_util(&["revcomp", input.path()], b"");

    assert!(output.status.success());
    let mut expected = b">long\n".to_vec();
    let transformed = sequence
        .iter()
        .rev()
        .map(|byte| match byte {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' | b'-' => *byte,
            _ => unreachable!(),
        })
        .collect::<Vec<_>>();
    for chunk in transformed.chunks(60) {
        expected.extend_from_slice(chunk);
        expected.push(b'\n');
    }
    assert_eq!(output.stdout, expected);
}

#[test]
fn revcomp_rejects_invalid_and_mixed_dna_rna_symbols_without_replacing_output() {
    for (sequence, message) in [
        (&b"ACGZ"[..], "invalid nucleotide 'Z'"),
        (&b"ACTU"[..], "contains both T and U"),
    ] {
        let input = TemporaryFile::new(&[b">record\n".as_slice(), sequence, b"\n"].concat());
        let output_file = TemporaryFile::new(b"previous output");
        let output = run_fasta_util(
            &["revcomp", input.path(), "--output", output_file.path()],
            b"",
        );

        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains(message));
        assert_eq!(output_file.read(), b"previous output");
    }
}

#[test]
fn revcomp_rejects_sequence_before_header_and_empty_header_ids() {
    for (contents, message) in [
        (
            &b"ACGT\n>record\nACGT\n"[..],
            "sequence data appears before",
        ),
        (&b">  \nACGT\n"[..], "record identifier is empty"),
    ] {
        let input = TemporaryFile::new(contents);
        let output = run_fasta_util(&["revcomp", input.path()], b"");

        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains(message));
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn revcomp_rejects_zero_output_line_width() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let output = run_fasta_util(&["revcomp", input.path(), "--chars-per-line", "0"], b"");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("greater than zero"));
}

#[test]
fn grep_selects_complete_records_by_header_text() {
    let input =
        TemporaryFile::new(b">BRCA1 description\r\nACGT\r\n>other\r\nBRCA\r\n>BRCA2\r\nTGCA\r\n");
    let output = run_fasta_util(&["grep", input.path(), "BRCA"], b"");

    assert!(output.status.success());
    assert_eq!(
        output.stdout,
        b">BRCA1 description\r\nACGT\r\n>BRCA2\r\nTGCA\r\n"
    );
}

#[test]
fn grep_supports_case_insensitive_and_inverted_matching() {
    let input = TemporaryFile::new(b">BRCA1\nACGT\n>other\nTTAA\n>brca2\nTGCA\n");
    let insensitive = run_fasta_util(&["grep", input.path(), "brca", "--ignore-case"], b"");
    let inverted = run_fasta_util(
        &[
            "grep",
            input.path(),
            "BRCA",
            "--ignore-case",
            "--invert-match",
        ],
        b"",
    );

    assert!(insensitive.status.success());
    assert_eq!(insensitive.stdout, b">BRCA1\nACGT\n>brca2\nTGCA\n");
    assert!(inverted.status.success());
    assert_eq!(inverted.stdout, b">other\nTTAA\n");
}

#[test]
fn grep_writes_matches_to_a_file_and_rejects_same_input_output() {
    let input = TemporaryFile::new(b">match\nACGT\n>skip\nTTAA\n");
    let output_file = TemporaryFile::new(b"old output");
    let output = run_fasta_util(
        &[
            "grep",
            input.path(),
            "match",
            "--output",
            output_file.path(),
        ],
        b"",
    );
    let same_file = run_fasta_util(
        &["grep", input.path(), "match", "--output", input.path()],
        b"",
    );

    assert!(output.status.success());
    assert_eq!(output_file.read(), b">match\nACGT\n");
    assert!(!same_file.status.success());
    assert_eq!(input.read(), b">match\nACGT\n>skip\nTTAA\n");
}

#[test]
fn grep_rejects_empty_pattern_and_data_before_first_record() {
    let input = TemporaryFile::new(b"ACGT\n>record\nACGT\n");
    let output = run_fasta_util(&["grep", input.path(), "record"], b"");
    let empty_pattern_input = TemporaryFile::new(b">record\nACGT\n");
    let empty_pattern = run_fasta_util(&["grep", empty_pattern_input.path(), ""], b"");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("sequence data appears before"));
    assert!(output.stdout.is_empty());
    assert!(!empty_pattern.status.success());
    assert!(String::from_utf8_lossy(&empty_pattern.stderr).contains("pattern cannot be empty"));
}

#[test]
fn format_wraps_sequences_and_normalizes_line_endings() {
    let input = TemporaryFile::new(b">first description\r\nACG\r\nTTA\n>second\nCCGGTT\n");
    let output = run_fasta_util(&["format", "--width", "4", input.path()], b"");

    assert!(output.status.success());
    assert_eq!(
        output.stdout,
        b">first description\nACGT\nTA\n>second\nCCGG\nTT\n"
    );
}

#[test]
fn format_width_zero_writes_one_sequence_line_per_record() {
    let input = TemporaryFile::new(b">first\nAC\nGT\n>second\nTT\nAA\n");
    let output = run_fasta_util(&["format", "--width", "0", input.path()], b"");

    assert!(output.status.success());
    assert_eq!(output.stdout, b">first\nACGT\n>second\nTTAA\n");
}

#[test]
fn format_changes_sequence_case_removes_gaps_and_trims_headers() {
    let input = TemporaryFile::new(b">  first description  \r\nac-g\r\nTt-a\n> second \nN--\n");
    let uppercase = run_fasta_util(
        &[
            "format",
            input.path(),
            "--width",
            "0",
            "--uppercase",
            "--remove-gaps",
            "--trim-header",
        ],
        b"",
    );
    let lowercase = run_fasta_util(&["format", input.path(), "--lowercase"], b"");

    assert!(uppercase.status.success());
    assert_eq!(
        uppercase.stdout,
        b">first description\nACGTTA\n>second\nN\n"
    );
    assert!(lowercase.status.success());
    assert_eq!(
        lowercase.stdout,
        b">  first description  \nac-gtt-a\n> second \nn--\n"
    );
}

#[test]
fn format_keeps_empty_records_and_emits_no_empty_sequence_lines() {
    let input = TemporaryFile::new(b">empty\n>gaps\n---\n>sequence\nA-C\n");
    let output = run_fasta_util(
        &["format", input.path(), "--remove-gaps", "--width", "0"],
        b"",
    );

    assert!(output.status.success());
    assert_eq!(output.stdout, b">empty\n>gaps\n>sequence\nAC\n");
}

#[test]
fn format_writes_to_a_file_and_rejects_using_input_as_output() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let output_file = TemporaryFile::new(b"previous output");
    let output = run_fasta_util(
        &[
            "format",
            input.path(),
            "--width",
            "2",
            "--output",
            output_file.path(),
        ],
        b"",
    );
    let same_file = run_fasta_util(&["format", input.path(), "--output", input.path()], b"");

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(output_file.read(), b">record\nAC\nGT\n");
    assert!(!same_file.status.success());
    assert_eq!(input.read(), b">record\nACGT\n");
}

#[test]
fn format_rejects_invalid_structure_without_replacing_output() {
    let input = TemporaryFile::new(b">valid\nACGT\n>invalid\nAC GT\n");
    let output_file = TemporaryFile::new(b"keep this output");
    let output = run_fasta_util(
        &["format", input.path(), "--output", output_file.path()],
        b"",
    );
    let before_header = TemporaryFile::new(b"ACGT\n>record\nACGT\n");
    let before_header_output = run_fasta_util(&["format", before_header.path()], b"");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("whitespace in sequence"));
    assert_eq!(output_file.read(), b"keep this output");
    assert!(!before_header_output.status.success());
    assert!(
        String::from_utf8_lossy(&before_header_output.stderr)
            .contains("sequence data appears before")
    );
    assert!(before_header_output.stdout.is_empty());
}

#[test]
fn format_rejects_conflicting_case_options_and_empty_ids() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let conflict = run_fasta_util(&["format", input.path(), "--uppercase", "--lowercase"], b"");
    let empty_id_input = TemporaryFile::new(b">  \nACGT\n");
    let empty_id = run_fasta_util(&["format", empty_id_input.path()], b"");

    assert!(!conflict.status.success());
    assert!(!empty_id.status.success());
    assert!(String::from_utf8_lossy(&empty_id.stderr).contains("record identifier is empty"));
}

#[test]
fn locate_reports_one_based_inclusive_coordinates_on_both_strands() {
    let input = TemporaryFile::new(b">chr1 description\nAAATA\nAAGG\n>chr3\nCTTTA\nTTG\n");
    let output = run_fasta_util(&["locate", input.path(), "AATAAA"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"chr1\t2\t7\t+\nchr3\t2\t7\t-\n");
}

#[test]
fn locate_streams_wrapped_stdin_and_reports_sequence_errors() {
    let output = run_fasta_util(&["locate", "-", "ACG"], b">first description\r\nA\r\nCGT\n");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"first\t1\t3\t+\nfirst\t2\t4\t-\n");

    let invalid = run_fasta_util(&["locate", "-", "ACG"], b">record\nAC?\n");
    assert!(!invalid.status.success());
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("line 2: invalid nucleotide '?'"));
}

#[test]
fn locate_supports_degenerate_iupac_motifs_across_wrapped_lines() {
    let input = TemporaryFile::new(b">seq description\natg\ngacTAA\n>other\nATGCCCTAA\n");
    let output = run_fasta_util(&["locate", input.path(), "ATGNNNTAA"], b"");

    assert!(output.status.success());
    assert_eq!(output.stdout, b"seq\t1\t9\t+\nother\t1\t9\t+\n");
}

#[test]
fn locate_finds_overlapping_hits_and_reports_palindromic_strands() {
    let input = TemporaryFile::new(b">repeat\nAAAA\n>palindrome\nAT\n");
    let overlapping = run_fasta_util(&["locate", input.path(), "AA"], b"");
    let palindrome = run_fasta_util(&["locate", input.path(), "AT"], b"");

    assert!(overlapping.status.success());
    assert_eq!(
        overlapping.stdout,
        b"repeat\t1\t2\t+\nrepeat\t2\t3\t+\nrepeat\t3\t4\t+\n"
    );
    assert!(palindrome.status.success());
    assert_eq!(
        palindrome.stdout,
        b"palindrome\t1\t2\t+\npalindrome\t1\t2\t-\n"
    );
}

#[test]
fn locate_reports_both_strands_for_palindromic_motif_with_one_mismatch() {
    let input = TemporaryFile::new(b">exact\nATAT\n>one-mismatch\nATCT\n");
    let output = run_fasta_util(
        &["locate", input.path(), "ATAT", "--max-mismatch", "1"],
        b"",
    );

    assert!(output.status.success());
    assert_eq!(
        output.stdout,
        b"exact\t1\t4\t+\nexact\t1\t4\t-\none-mismatch\t1\t4\t+\none-mismatch\t1\t4\t-\n"
    );
}

#[test]
fn locate_does_not_match_across_record_boundaries() {
    let input = TemporaryFile::new(b">first\nAAAA\n>second\nCCCC\n");
    let output = run_fasta_util(&["locate", input.path(), "AAAACCCC"], b"");

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
}

#[test]
fn locate_matches_long_motifs_without_crossing_record_boundaries() {
    let motif = "A".repeat(65);
    let short_sequence = "A".repeat(40);
    let input = TemporaryFile::new(
        format!(">first\n{short_sequence}\n>second\n{short_sequence}\n>exact\n{motif}\n")
            .as_bytes(),
    );
    let output = run_fasta_util(&["locate", input.path(), &motif], b"");

    assert!(output.status.success());
    assert_eq!(output.stdout, b"exact\t1\t65\t+\n");
}

#[test]
fn locate_uses_bounded_mismatch_search_for_long_motifs() {
    let motif = "A".repeat(65);
    let two_mismatches = format!("{}CC{}", "A".repeat(32), "A".repeat(31));
    let three_mismatches = format!("{}CCC{}", "A".repeat(32), "A".repeat(30));
    let input = TemporaryFile::new(
        format!(">exact\n{motif}\n>two\n{two_mismatches}\n>three\n{three_mismatches}\n").as_bytes(),
    );

    let two_allowed = run_fasta_util(
        &["locate", input.path(), &motif, "--max-mismatch", "2"],
        b"",
    );
    let three_allowed = run_fasta_util(
        &["locate", input.path(), &motif, "--max-mismatch", "3"],
        b"",
    );

    assert!(two_allowed.status.success());
    assert_eq!(two_allowed.stdout, b"exact\t1\t65\t+\ntwo\t1\t65\t+\n");
    assert!(three_allowed.status.success());
    assert_eq!(
        three_allowed.stdout,
        b"exact\t1\t65\t+\ntwo\t1\t65\t+\nthree\t1\t65\t+\n"
    );
}

#[test]
fn locate_matches_every_window_when_mismatch_limit_equals_motif_length() {
    let motif = "A".repeat(65);
    let input = TemporaryFile::new(
        format!(
            ">seq description\n{}\n{}\n>second\n{}\n",
            "C".repeat(32),
            "C".repeat(33),
            "G".repeat(65)
        )
        .as_bytes(),
    );
    let output = run_fasta_util(
        &["locate", input.path(), &motif, "--max-mismatch", "65"],
        b"",
    );

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"seq\t1\t65\t+\nseq\t1\t65\t-\nsecond\t1\t65\t+\nsecond\t1\t65\t-\n"
    );
}

#[test]
fn locate_allows_a_bounded_number_of_mismatches() {
    let input = TemporaryFile::new(b">seq\nAATCAA\n");
    let exact = run_fasta_util(&["locate", input.path(), "AATAAA"], b"");
    let one_mismatch = run_fasta_util(
        &["locate", input.path(), "AATAAA", "--max-mismatch", "1"],
        b"",
    );

    assert!(exact.status.success());
    assert!(exact.stdout.is_empty());
    assert!(one_mismatch.status.success());
    assert_eq!(one_mismatch.stdout, b"seq\t1\t6\t+\n");

    let excessive = run_fasta_util(
        &["locate", input.path(), "AATAAA", "--max-mismatch", "7"],
        b"",
    );
    assert!(!excessive.status.success());
    assert!(String::from_utf8_lossy(&excessive.stderr).contains("greater than the motif length"));
}

#[test]
fn locate_rejects_invalid_motifs_and_sequences_without_replacing_output() {
    let input = TemporaryFile::new(b">seq\nACGZ\n");
    let output_file = TemporaryFile::new(b"keep output");
    let invalid_sequence = run_fasta_util(
        &[
            "locate",
            input.path(),
            "ACG",
            "--output",
            output_file.path(),
        ],
        b"",
    );
    let invalid_motif = run_fasta_util(&["locate", input.path(), "ACZ"], b"");
    let empty_motif = run_fasta_util(&["locate", input.path(), ""], b"");

    assert!(!invalid_sequence.status.success());
    assert!(String::from_utf8_lossy(&invalid_sequence.stderr).contains("invalid nucleotide 'Z'"));
    assert_eq!(output_file.read(), b"keep output");
    assert!(!invalid_motif.status.success());
    assert!(String::from_utf8_lossy(&invalid_motif.stderr).contains("invalid IUPAC"));
    assert!(!empty_motif.status.success());
    assert!(String::from_utf8_lossy(&empty_motif.stderr).contains("motif cannot be empty"));
}

#[test]
fn get_extracts_requested_records_from_a_multi_fasta() {
    let input = TemporaryFile::new(b">chr1 description\nACGT\n>chr2\nTTAA\n>chr3 extra\nGGCC\n");
    let output = run_fasta_util(&["get", input.path(), "chr3", "chr1"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b">chr1 description\nACGT\n>chr3 extra\nGGCC\n"
    );
}

#[test]
fn get_reads_ids_from_a_file() {
    let input = TemporaryFile::new(b">chr1\nACGT\n>chr2\nTTAA\n");
    let ids = TemporaryFile::new(b"chr2\n\nchr1\n");
    let output = run_fasta_util(&["get", input.path(), "--ids", ids.path()], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">chr1\nACGT\n>chr2\nTTAA\n");
}

#[test]
fn get_extracts_one_based_inclusive_region_from_streaming_input() {
    let input = TemporaryFile::new(b">chr1 description\nAACCGGTA\n>chr2\nTTAA\n");
    let output = run_fasta_util(&["get", input.path(), "chr1:2-5"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">chr1:2-5\nACCG\n");
}

#[test]
fn get_uses_a_sidecar_fai_for_region_extraction() {
    let input = TemporaryFile::new(b">chr1 description\nAACCGG\nTTAA\n>chr2\nGGCC\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let index = run_fasta_util(&["index", input.path()], b"");
    assert!(index.status.success());

    let output = run_fasta_util(&["get", input.path(), "chr1:3-8"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">chr1:3-8\nCCGGTT\n");
    fs::remove_file(index_path).unwrap();
}

#[test]
fn get_can_bypass_an_invalid_adjacent_fai_for_record_and_global_ranges() {
    let input = TemporaryFile::new(b">chr1 description\nAACCGG\nTTAA\n>chr2\nGGCC\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let index = run_fasta_util(&["index", input.path()], b"");
    assert!(index.status.success());
    fs::write(&index_path, b"invalid index\n").unwrap();

    let automatic = run_fasta_util(&["get", input.path(), "chr1:3-8"], b"");
    assert!(!automatic.status.success());

    let record_region = run_fasta_util(&["get", input.path(), "chr1:3-8", "--no-fai-index"], b"");
    assert!(
        record_region.status.success(),
        "{}",
        String::from_utf8_lossy(&record_region.stderr)
    );
    assert_eq!(record_region.stdout, b">chr1:3-8\nCCGGTT\n");

    let global_range = run_fasta_util(&["get", input.path(), "8-13", "--no-fai-index"], b"");
    assert!(
        global_range.status.success(),
        "{}",
        String::from_utf8_lossy(&global_range.stderr)
    );
    assert_eq!(global_range.stdout, b">chr1 description\nTAA\n>chr2\nGGC\n");

    fs::remove_file(index_path).unwrap();
}

#[test]
fn get_rejects_conflicting_fai_options() {
    let input = TemporaryFile::new(b">chr1\nACGT\n");
    let output = run_fasta_util(
        &[
            "get",
            input.path(),
            "chr1",
            "--fai-index",
            "custom.fai",
            "--no-fai-index",
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot be used with"));
}

#[test]
fn get_uses_a_sidecar_fai_for_global_range_extraction() {
    let input = TemporaryFile::new(b">first description\nAACG\nTTGC\n>second\nCAAA\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let index = run_fasta_util(&["index", input.path()], b"");
    assert!(index.status.success());

    let output = run_fasta_util(&["get", input.path(), "3-10"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b">first description\nCGTTGC\n>second\nCA\n");
    fs::remove_file(index_path).unwrap();
}

#[test]
fn get_region_coordinates_include_both_ends_and_clip_at_sequence_length() {
    let input = TemporaryFile::new(b">chr1\nACGT\n");
    let first = run_fasta_util(&["get", input.path(), "chr1:1-1"], b"");
    let last = run_fasta_util(&["get", input.path(), "chr1:4-99"], b"");

    assert!(first.status.success());
    assert_eq!(first.stdout, b">chr1:1-1\nA\n");
    assert!(last.status.success());
    assert_eq!(last.stdout, b">chr1:4-99\nT\n");
}

#[test]
fn get_region_output_matches_between_streaming_and_explicit_indexed_paths() {
    let input = TemporaryFile::new(b">chr1 description\r\nAACG\r\nTTGC\r\nCA");
    let index_path = TemporaryFile::new(b"");
    let index = run_fasta_util(&["index", input.path()], b"");
    assert!(
        index.status.success(),
        "{}",
        String::from_utf8_lossy(&index.stderr)
    );
    let sidecar = PathBuf::from(format!("{}.fai", input.path()));
    fs::copy(&sidecar, index_path.path()).unwrap();
    fs::remove_file(sidecar).unwrap();

    let streamed = run_fasta_util(&["get", input.path(), "chr1:3-8"], b"");
    let indexed = run_fasta_util(
        &[
            "get",
            input.path(),
            "chr1:3-8",
            "--fai-index",
            index_path.path(),
        ],
        b"",
    );

    assert!(streamed.status.success());
    assert!(
        indexed.status.success(),
        "{}",
        String::from_utf8_lossy(&indexed.stderr)
    );
    assert_eq!(streamed.stdout, b">chr1:3-8\nCGTTGC\n");
    assert_eq!(indexed.stdout, streamed.stdout);
}

#[test]
fn get_wraps_protein_output_and_normalizes_crlf_in_both_paths() {
    let input = TemporaryFile::new(b">protein description\r\nACDE\r\nFG*\r\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let index = run_fasta_util(&["index", input.path()], b"");
    assert!(index.status.success());

    let streamed = run_fasta_util(
        &[
            "get",
            input.path(),
            "protein",
            "--sequence-type",
            "protein",
            "--chars-per-line",
            "3",
        ],
        b"",
    );
    let indexed = run_fasta_util(
        &[
            "get",
            input.path(),
            "protein",
            "--sequence-type",
            "protein",
            "--chars-per-line",
            "3",
            "--fai-index",
            index_path.to_str().unwrap(),
        ],
        b"",
    );

    assert!(streamed.status.success());
    assert!(indexed.status.success());
    assert_eq!(streamed.stdout, b">protein description\nACD\nEFG\n*\n");
    assert_eq!(indexed.stdout, streamed.stdout);
    fs::remove_file(index_path).unwrap();
}

#[test]
fn get_validates_only_the_selected_region_in_both_paths() {
    let input = TemporaryFile::new(b">chr1\nACGTZ\n");
    let index_path = PathBuf::from(format!("{}.fai", input.path()));
    let index = run_fasta_util(&["index", input.path()], b"");
    assert!(index.status.success());

    for extra_args in [
        Vec::<&str>::new(),
        vec!["--fai-index", index_path.to_str().unwrap()],
    ] {
        let valid = run_fasta_util(
            &[&["get", input.path(), "chr1:1-4"][..], &extra_args].concat(),
            b"",
        );
        let invalid = run_fasta_util(
            &[&["get", input.path(), "chr1:5-5"][..], &extra_args].concat(),
            b"",
        );
        assert!(valid.status.success());
        assert_eq!(valid.stdout, b">chr1:1-4\nACGT\n");
        assert!(!invalid.status.success());
        assert!(String::from_utf8_lossy(&invalid.stderr).contains("invalid nucleic acid: 'Z'"));
    }
    fs::remove_file(index_path).unwrap();
}

#[test]
fn get_rejects_bad_region_syntax_and_duplicate_ids_file_entries() {
    let input = TemporaryFile::new(b">chr1\nACGT\n");
    let bad_region = run_fasta_util(&["get", input.path(), "chr1:2-"], b"");
    let ids = TemporaryFile::new(b"chr1\nchr1:2-3\n");
    let duplicate_ids = run_fasta_util(&["get", input.path(), "--ids", ids.path()], b"");

    assert!(!bad_region.status.success());
    assert!(String::from_utf8_lossy(&bad_region.stderr).contains("invalid region `chr1:2-`"));
    assert!(!duplicate_ids.status.success());
    assert!(String::from_utf8_lossy(&duplicate_ids.stderr).contains("requested more than once"));
}

#[test]
fn get_rejects_missing_ids_and_preserves_existing_output() {
    let input = TemporaryFile::new(b">chr1\nACGT\n");
    let empty_ids = TemporaryFile::new(b"\n  \n");
    let output_file = TemporaryFile::new(b"existing output\n");
    let no_queries = run_fasta_util(&["get", input.path(), "--ids", empty_ids.path()], b"");
    let missing_region = run_fasta_util(
        &[
            "get",
            input.path(),
            "missing:1-2",
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(!no_queries.status.success());
    assert!(String::from_utf8_lossy(&no_queries.stderr).contains("provide one or more"));
    assert!(!missing_region.status.success());
    assert!(
        String::from_utf8_lossy(&missing_region.stderr).contains("record `missing` was not found")
    );
    assert_eq!(output_file.read(), b"existing output\n");
}

#[test]
fn get_missing_id_does_not_replace_an_output_file() {
    let input = TemporaryFile::new(b">chr1\nACGT\n");
    let output_file = TemporaryFile::new(b"existing output\n");
    let output = run_fasta_util(
        &[
            "get",
            input.path(),
            "missing",
            "--output",
            output_file.path(),
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("record `missing` was not found"));
    assert_eq!(output_file.read(), b"existing output\n");
}

#[test]
fn get_rejects_duplicate_record_ids() {
    let input = TemporaryFile::new(b">chr1\nACGT\n");
    let output = run_fasta_util(&["get", input.path(), "chr1", "chr1:2-3"], b"");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("record `chr1` is requested more than once")
    );
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
fn len_counts_protein_symbols_from_stdin() {
    let output = run_fasta_util(
        &["len", "--sequence-type", "protein"],
        b">protein\nACDEFGHIKL\nMNPQRSTVWY\nBZJXUO*-\n",
    );

    assert!(output.status.success());
    assert_eq!(output.stdout, b"28\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn get_handles_protein_symbols_and_preserves_case() {
    let input = TemporaryFile::new(b">protein\nacDEFGHIK*\n");
    let output = run_fasta_util(
        &[
            "get",
            input.path(),
            "--sequence-type",
            "protein",
            "3-9",
            "--chars-per-line",
            "4",
        ],
        b"",
    );

    assert!(output.status.success());
    assert_eq!(output.stdout, b">protein\nDEFG\nHIK\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn indexed_protein_get_matches_streaming_get() {
    let input = TemporaryFile::new(b">protein description\nACDE\nBZJU\nOX*-\n");
    let index = TemporaryFile::new(b"protein\t12\t21\t4\t5\n");
    let args = [
        "get",
        "--sequence-type",
        "protein",
        input.path(),
        "3-11",
        "--chars-per-line",
        "4",
    ];
    let normal = run_fasta_util(&args, b"");
    let indexed = run_fasta_util(
        &[
            "get",
            "--sequence-type",
            "protein",
            input.path(),
            "--fai-index",
            index.path(),
            "3-11",
            "--chars-per-line",
            "4",
        ],
        b"",
    );

    assert!(
        normal.status.success(),
        "{}",
        String::from_utf8_lossy(&normal.stderr)
    );
    assert!(
        indexed.status.success(),
        "{}",
        String::from_utf8_lossy(&indexed.stderr)
    );
    assert_eq!(indexed.stdout, normal.stdout);
    assert_eq!(normal.stdout, b">protein description\nDEBZ\nJUOX\n*\n");
}

#[test]
fn protein_mode_rejects_non_protein_symbols() {
    let output = run_fasta_util(&["len", "--sequence-type", "protein"], b">protein\nAC.D\n");

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid protein symbol"));
}

#[test]
fn get_writes_selected_global_range_to_stdout() {
    let input = TemporaryFile::new(b">record\nACGT\nNU-\n");
    let output = run_fasta_util(&["get", input.path(), "3-5", "--chars-per-line", "2"], b"");

    assert!(output.status.success());
    assert_eq!(output.stdout, b">record\nGT\nN\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn get_global_range_crosses_record_boundaries_and_preserves_each_header() {
    let input = TemporaryFile::new(b">first description\nACGT\n>second description\nTGCA\n");
    let output = run_fasta_util(&["get", input.path(), "3-6"], b"");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b">first description\nGT\n>second description\nTG\n"
    );
}

#[test]
fn get_rejects_invalid_or_combined_global_ranges() {
    let input = TemporaryFile::new(b">first\nACGT\n>second\nTGCA\n");

    for range in ["0-2", "4-3", "1-"] {
        let output = run_fasta_util(&["get", input.path(), range], b"");

        assert!(
            !output.status.success(),
            "range {range} unexpectedly passed"
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("invalid global range"),
            "range {range}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let combined = run_fasta_util(&["get", input.path(), "1-2", "first"], b"");
    assert!(!combined.status.success());
    assert!(
        String::from_utf8_lossy(&combined.stderr)
            .contains("a global range cannot be combined with record IDs")
    );
}

#[test]
fn get_reads_a_crlf_file_and_writes_normalized_fasta_to_a_file() {
    let input = TemporaryFile::new(b">first\r\nACGT\r\n>second\r\nNU-\r\n");
    let output_file = TemporaryFile::new(b"old contents");
    let output = run_fasta_util(
        &["get", input.path(), "1-7", "--output", output_file.path()],
        b"",
    );

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    assert_eq!(output_file.read(), b">first\nACGT\n>second\nNU-\n");
}

#[test]
fn indexed_get_matches_streaming_get_for_wrapped_multi_record_fasta() {
    let (input, index) = indexed_fasta_fixture();
    let cases = [
        ("1-14", "2"),
        ("4-10", "4"),
        ("9-9", "3"),
        ("11-14", "5"),
        ("15-15", "4"),
        ("100-100", "4"),
    ];

    for (range, chars_per_line) in cases {
        let normal = run_fasta_util(
            &[
                "get",
                input.path(),
                range,
                "--chars-per-line",
                chars_per_line,
            ],
            b"",
        );
        let indexed = run_fasta_util(
            &[
                "get",
                input.path(),
                "--fai-index",
                index.path(),
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
fn indexed_get_matches_streaming_across_record_layouts_and_line_endings() {
    let records = [
        ("empty", b"".as_slice(), 4),
        ("single", b"a".as_slice(), 4),
        ("short", b"ACGT".as_slice(), 2),
        ("odd", b"tgcaN".as_slice(), 3),
        ("long", b"ACGTagctN-".as_slice(), 4),
        ("last", b"u".as_slice(), 1),
    ];
    let ranges = ["1-21", "1-1", "2-5", "4-10", "9-18", "21-21"];
    let line_widths = ["1", "3", "8"];

    for line_ending in [b"\n".as_slice(), b"\r\n".as_slice()] {
        let (input, index) = indexed_fasta_with_layout(&records, line_ending);
        for range in ranges {
            for line_width in line_widths {
                let normal = run_fasta_util(
                    &["get", input.path(), range, "--chars-per-line", line_width],
                    b"",
                );
                let indexed = run_fasta_util(
                    &[
                        "get",
                        input.path(),
                        "--fai-index",
                        index.path(),
                        range,
                        "--chars-per-line",
                        line_width,
                    ],
                    b"",
                );

                assert!(normal.status.success(), "normal {range}: {normal:?}");
                assert!(indexed.status.success(), "indexed {range}: {indexed:?}");
                assert_eq!(indexed.stdout, normal.stdout, "{range}, width {line_width}");
            }
        }
    }
}

#[test]
fn indexed_get_reads_headers_longer_than_the_header_scan_buffer() {
    let first_header = b">previous\n";
    let second_name = "longname".repeat(800);
    let second_header = format!(">{second_name} long description\n");
    let first_offset = first_header.len();
    let second_offset = first_offset + 5 + second_header.len();
    let mut fasta = first_header.to_vec();
    fasta.extend_from_slice(b"ACGT\n");
    fasta.extend_from_slice(second_header.as_bytes());
    fasta.extend_from_slice(b"TT\n");
    let input = TemporaryFile::new(&fasta);
    let index = TemporaryFile::new(
        format!("previous\t4\t{first_offset}\t4\t5\n{second_name}\t2\t{second_offset}\t2\t3\n")
            .as_bytes(),
    );

    let normal = run_fasta_util(&["get", input.path(), "1-6"], b"");
    let indexed = run_fasta_util(
        &["get", input.path(), "--fai-index", index.path(), "1-6"],
        b"",
    );

    assert!(normal.status.success(), "{:?}", normal.stderr);
    assert!(indexed.status.success(), "{:?}", indexed.stderr);
    assert_eq!(indexed.stdout, normal.stdout);
}

#[test]
fn indexed_get_matches_streaming_get_across_read_buffer_chunks() {
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
    let cases = [("58-70031", "73"), ("65511-140008", "61")];

    for (range, chars_per_line) in cases {
        let normal = run_fasta_util(
            &[
                "get",
                input.path(),
                range,
                "--chars-per-line",
                chars_per_line,
            ],
            b"",
        );
        let indexed = run_fasta_util(
            &[
                "get",
                input.path(),
                "--fai-index",
                index.path(),
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
fn indexed_get_validates_only_selected_sequence_bases() {
    let input = TemporaryFile::new(b">record description\nACXT\n");
    let index = TemporaryFile::new(b"record\t4\t20\t4\t5\n");
    let output = run_fasta_util(
        &["get", input.path(), "--fai-index", index.path(), "1-2"],
        b"",
    );

    assert!(output.status.success(), "{:?}", output.stderr);
    assert_eq!(output.stdout, b">record description\nAC\n");
}

#[test]
fn indexed_get_rejects_mismatched_index_without_replacing_output() {
    let (input, _) = indexed_fasta_fixture();
    let index = TemporaryFile::new(b"wrong-name\t8\t35\t4\t6\n");
    let output_file = TemporaryFile::new(b"keep this output");
    let output = run_fasta_util(
        &[
            "get",
            input.path(),
            "--fai-index",
            index.path(),
            "--output",
            output_file.path(),
            "wrong-name:1-8",
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("FAI record name"));
    assert_eq!(output_file.read(), b"keep this output");
}

#[test]
fn indexed_get_rejects_index_as_output_without_changing_it() {
    let (input, index) = indexed_fasta_fixture();
    let original_index = index.read();
    let output = run_fasta_util(
        &[
            "get",
            input.path(),
            "--fai-index",
            index.path(),
            "--output",
            index.path(),
            "empty:1-1",
        ],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("same file"));
    assert_eq!(index.read(), original_index);
}

#[cfg(unix)]
#[test]
fn indexed_get_rejects_linked_index_outputs_without_changing_the_index() {
    let (input, index) = indexed_fasta_fixture();
    let original_index = index.read();

    for link_kind in ["hard link", "symbolic link"] {
        let output_alias = TemporaryFile::new(b"temporary alias path");
        fs::remove_file(&output_alias.0).unwrap();
        match link_kind {
            "hard link" => fs::hard_link(&index.0, &output_alias.0).unwrap(),
            "symbolic link" => std::os::unix::fs::symlink(&index.0, &output_alias.0).unwrap(),
            _ => unreachable!(),
        }

        let output = run_fasta_util(
            &[
                "get",
                input.path(),
                "--fai-index",
                index.path(),
                "--output",
                output_alias.path(),
                "empty:1-1",
            ],
            b"",
        );

        assert!(
            !output.status.success(),
            "accepted {link_kind} as the output path"
        );
        assert!(String::from_utf8_lossy(&output.stderr).contains("same file"));
        assert_eq!(index.read(), original_index, "via {link_kind}");
        assert_eq!(output_alias.read(), original_index, "via {link_kind}");
    }
}

#[test]
fn indexed_get_rejects_malformed_and_out_of_bounds_indexes() {
    let (input, _) = indexed_fasta_fixture();
    for contents in [
        &b"record\t4\t7\t4\n"[..],
        &b"\t4\t7\t4\t4\n"[..],
        &b"record\tlength\t0\t4\t4\n"[..],
        &b"record\t4\t999999\t4\t4\n"[..],
        &b"record\t4\t12\t0\t0\n"[..],
    ] {
        let index = TemporaryFile::new(contents);
        let output = run_fasta_util(
            &[
                "get",
                input.path(),
                "--fai-index",
                index.path(),
                "record:1-1",
            ],
            b"",
        );
        assert!(!output.status.success(), "index {contents:?} was accepted");
        assert!(!String::from_utf8_lossy(&output.stderr).contains("panicked"));
    }
}

#[test]
fn get_requires_a_fasta_input() {
    let index = TemporaryFile::new(b"record\t4\t7\t4\t4\n");
    let output = run_fasta_util(&["get", "--fai-index", index.path()], b"");

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("required arguments were not provided")
    );
}

#[cfg(unix)]
#[test]
fn get_accepts_a_non_utf8_output_path() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let output_file = temporary_non_utf8_file(b"old contents");
    let output = run_fasta_util_with(b"", |command| {
        command
            .arg("get")
            .arg(&input.0)
            .arg("record")
            .arg("--output")
            .arg(&output_file.0);
    });

    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    assert_eq!(output_file.read(), b">record\nACGT\n");
}

#[test]
fn get_keeps_existing_output_when_input_validation_fails_after_partial_output() {
    let input = TemporaryFile::new(b">record\nACGT\nACX\n");
    let output_file = TemporaryFile::new(b"previous output");
    let output = run_fasta_util(
        &["get", input.path(), "1-7", "--output", output_file.path()],
        b"",
    );

    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid nucleic acid"));
    assert_eq!(output_file.read(), b"previous output");
}

#[test]
fn bounded_global_get_stops_before_invalid_sequence_after_the_range() {
    let input = TemporaryFile::new(b">record\nACGT\nACX\n");
    let output = run_fasta_util(&["get", input.path(), "1-2"], b"");

    assert!(output.status.success(), "{output:?}");
    assert!(output.stderr.is_empty());
    assert_eq!(output.stdout, b">record\nAC\n");
}

#[test]
fn get_does_not_create_output_when_input_validation_fails() {
    let input = TemporaryFile::new(b">record\nACX\n");
    let output_file = TemporaryFile::new(b"remove me");
    fs::remove_file(&output_file.0).unwrap();

    let output = run_fasta_util(
        &["get", input.path(), "1-3", "--output", output_file.path()],
        b"",
    );

    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid nucleic acid"));
    assert!(!output_file.0.exists());
}

#[test]
fn get_rejects_same_input_and_output_file_without_changing_it() {
    let input = TemporaryFile::new(b">record\nACGT\n");
    let output = run_fasta_util(&["get", input.path(), "1-4", "--output", input.path()], b"");

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
