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
