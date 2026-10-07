use crate::{SequenceType, for_each_reader_line};
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

#[derive(Clone, Copy, PartialEq, Eq)]
enum LineEnding {
    Lf,
    CrLf,
}

struct FaiLine {
    line_number: usize,
    length: usize,
    source: Vec<u8>,
}

#[derive(Default)]
struct Record {
    header_line: usize,
    header: Vec<u8>,
    id: Vec<u8>,
    symbol_count: usize,
    first_sequence_line: Option<(usize, usize)>,
    pending_short_line: Option<FaiLine>,
    line_ending: Option<LineEnding>,
}

#[derive(Default)]
struct Reporter {
    errors: usize,
}

impl Reporter {
    fn at(
        &mut self,
        path: &Path,
        line_number: usize,
        column: usize,
        source: &[u8],
        message: &str,
        note: Option<&str>,
    ) {
        self.errors += 1;
        let source = String::from_utf8_lossy(source);
        let width = line_number.to_string().len();
        eprintln!("error: {message}");
        eprintln!(" --> {}:{line_number}:{column}", path.display());
        eprintln!("  |");
        eprintln!("{line_number:>width$} | {source}");
        eprintln!("  | {}^", " ".repeat(column.saturating_sub(1)));
        if let Some(note) = note {
            eprintln!("  = note: {note}");
        }
        eprintln!();
    }

    fn global(&mut self, message: &str) {
        self.errors += 1;
        eprintln!("error: {message}");
        eprintln!();
    }
}

pub(super) fn run(path: Option<&Path>, sequence_type: SequenceType) -> io::Result<bool> {
    let display_path = path
        .filter(|path| *path != Path::new("-"))
        .map(Path::to_path_buf)
        .unwrap_or_else(|| Path::new("stdin").to_path_buf());

    if let Some(path) = path.filter(|path| *path != Path::new("-")) {
        let file = File::open(path)?;
        run_reader(BufReader::new(file), &display_path, sequence_type)
    } else {
        let stdin = io::stdin();
        run_reader(stdin.lock(), &display_path, sequence_type)
    }
}

fn run_reader(
    reader: impl BufRead,
    display_path: &Path,
    sequence_type: SequenceType,
) -> io::Result<bool> {
    let mut reporter = Reporter::default();
    let mut records = 0usize;
    let mut current_record: Option<Record> = None;
    let mut ids = HashMap::<Vec<u8>, (usize, Vec<u8>)>::new();
    let mut first_t: Option<(usize, usize, Vec<u8>)> = None;
    let mut first_u: Option<(usize, usize, Vec<u8>)> = None;

    for_each_reader_line(reader, |line_number, buffer| {
        let (line, line_ending) = strip_line_ending(buffer);

        if line.first() == Some(&b'>') {
            if let Some(record) = current_record.take() {
                finish_record(display_path, record, &mut reporter);
            }
            records += 1;
            let header = line.to_vec();
            let id = line[1..]
                .trim_ascii_start()
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default()
                .to_vec();
            if id.is_empty() {
                reporter.at(
                    display_path,
                    line_number,
                    2,
                    line,
                    "record identifier is empty",
                    None,
                );
            } else if let Some((first_line, first_header)) = ids.get(&id) {
                let id_text = String::from_utf8_lossy(&id);
                let note = format!(
                    "identifier `{id_text}` first appeared on line {first_line}: {}",
                    String::from_utf8_lossy(first_header)
                );
                reporter.at(
                    display_path,
                    line_number,
                    2,
                    line,
                    &format!("duplicate record identifier `{id_text}`"),
                    Some(&note),
                );
            } else {
                ids.insert(id.clone(), (line_number, header.clone()));
            }
            current_record = Some(Record {
                header_line: line_number,
                header,
                id,
                ..Record::default()
            });
            return Ok(());
        }

        if current_record.is_none() {
            if !line.is_empty() {
                reporter.at(
                    display_path,
                    line_number,
                    1,
                    line,
                    "sequence data appears before the first `>` record",
                    None,
                );
            }
            return Ok(());
        }

        let record = current_record.as_mut().expect("record was checked above");
        if line.is_empty() {
            reporter.at(
                display_path,
                line_number,
                1,
                line,
                "blank sequence line is incompatible with `.fai` indexing",
                None,
            );
            return Ok(());
        }

        if let Some(ending) = line_ending {
            if record
                .line_ending
                .is_some_and(|previous| previous != ending)
            {
                reporter.at(
                    display_path,
                    line_number,
                    line.len().saturating_add(1),
                    line,
                    "mixed LF and CRLF line endings within a record",
                    None,
                );
            } else {
                record.line_ending = Some(ending);
            }
        }

        record.symbol_count = record.symbol_count.saturating_add(line.len());
        check_fai_line(display_path, record, line_number, line, &mut reporter);

        for (index, byte) in line.iter().copied().enumerate() {
            let column = index + 1;
            if byte.is_ascii_whitespace() {
                reporter.at(
                    display_path,
                    line_number,
                    column,
                    line,
                    "whitespace is not allowed in sequence lines",
                    None,
                );
                continue;
            }

            if sequence_type == SequenceType::Nucleotide {
                if byte == b'T' || byte == b't' {
                    first_t.get_or_insert((line_number, column, line.to_vec()));
                } else if byte == b'U' || byte == b'u' {
                    first_u.get_or_insert((line_number, column, line.to_vec()));
                }
            }

            let valid = match sequence_type {
                SequenceType::Nucleotide => crate::is_nucleic_acid(byte),
                SequenceType::Protein => crate::is_amino_acid(byte),
            };
            if !valid {
                let symbol = char::from(byte);
                let kind = match sequence_type {
                    SequenceType::Nucleotide => "nucleotide",
                    SequenceType::Protein => "protein",
                };
                reporter.at(
                    display_path,
                    line_number,
                    column,
                    line,
                    &format!("invalid {kind} '{symbol}'"),
                    None,
                );
            }
        }
        Ok(())
    })?;

    if let Some(record) = current_record.take() {
        finish_record(display_path, record, &mut reporter);
    }

    if records == 0 {
        reporter.global("no FASTA records found");
    }

    if sequence_type == SequenceType::Nucleotide
        && let (Some(t), Some(u)) = (&first_t, &first_u)
    {
        let (line, column, source) = if (t.0, t.1) <= (u.0, u.1) { u } else { t };
        reporter.at(
            display_path,
            *line,
            *column,
            source,
            "DNA and RNA symbols (`T` and `U`) are mixed",
            None,
        );
    }

    if reporter.errors == 0 {
        let kind = match sequence_type {
            SequenceType::Protein => "protein",
            SequenceType::Nucleotide => match (first_t.is_some(), first_u.is_some()) {
                (true, false) => "DNA",
                (false, true) => "RNA",
                (false, false) => "DNA/RNA ambiguous",
                (true, true) => unreachable!("mixed DNA/RNA symbols are reported above"),
            },
        };
        println!("OK: {records} records\ntype: {kind}");
    }

    Ok(reporter.errors == 0)
}

fn check_fai_line(
    path: &Path,
    record: &mut Record,
    line_number: usize,
    line: &[u8],
    reporter: &mut Reporter,
) {
    let Some((first_line, line_bases)) = record.first_sequence_line else {
        record.first_sequence_line = Some((line_number, line.len()));
        return;
    };

    if let Some(previous) = record.pending_short_line.take() {
        reporter.at(
            path,
            previous.line_number,
            previous.length.saturating_add(1),
            &previous.source,
            "non-final sequence line is shorter than the `.fai` line width",
            Some("all sequence lines except the final line of a record must have equal width"),
        );
    }

    if line.len() < line_bases {
        record.pending_short_line = Some(FaiLine {
            line_number,
            length: line.len(),
            source: line.to_vec(),
        });
    } else if line.len() > line_bases {
        reporter.at(
            path,
            line_number,
            line_bases.saturating_add(1),
            line,
            "sequence line is wider than the first line and cannot be represented by `.fai`",
            Some(&format!("first sequence line is line {first_line}")),
        );
    }
}

fn finish_record(path: &Path, record: Record, reporter: &mut Reporter) {
    if record.symbol_count == 0 {
        reporter.at(
            path,
            record.header_line,
            2,
            &record.header,
            "record has an empty sequence",
            Some(&format!(
                "record identifier: {}",
                String::from_utf8_lossy(&record.id)
            )),
        );
    }
}

fn strip_line_ending(line: &[u8]) -> (&[u8], Option<LineEnding>) {
    if let Some(without_lf) = line.strip_suffix(b"\n") {
        if let Some(without_cr) = without_lf.strip_suffix(b"\r") {
            (without_cr, Some(LineEnding::CrLf))
        } else {
            (without_lf, Some(LineEnding::Lf))
        }
    } else {
        (line, None)
    }
}
