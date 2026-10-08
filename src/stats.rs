use crate::{
    LinesInFile, StatsFormat, StatsSequenceType, for_each_reader_line, is_amino_acid,
    is_nucleic_acid, parallel, read_lines_from_file, strip_line_ending,
};
use std::{fs::File, io, num::NonZeroUsize, path::Path};

#[derive(Default)]
struct RecordStats {
    id: String,
    length: u64,
    gc: u64,
    n: u64,
}

#[derive(Default)]
struct Counts {
    length: u64,
    gc: u64,
    n: u64,
}

#[derive(Default)]
struct Summary {
    records: Vec<RecordStats>,
    record_count: u64,
    total_length: u64,
    total_gc: u64,
    total_n: u64,
    lengths: Vec<u64>,
    min_length: u64,
    max_length: u64,
    n50: u64,
    saw_t: bool,
    saw_u: bool,
    saw_protein_only: bool,
}

struct StatsCollector {
    summary: Summary,
    current: Option<(Option<String>, Counts)>,
    saw_record: bool,
    sequence_type: StatsSequenceType,
    each: bool,
}

pub(super) fn run(
    path: Option<&Path>,
    each: bool,
    format: StatsFormat,
    sequence_type: StatsSequenceType,
    threads: Option<NonZeroUsize>,
) -> io::Result<()> {
    let summary = match path {
        None => {
            let stdin = io::stdin();
            collect_stats_from_reader(stdin.lock(), sequence_type, each)?
        }
        Some(path) if path == Path::new("-") => {
            let stdin = io::stdin();
            collect_stats_from_reader(stdin.lock(), sequence_type, each)?
        }
        Some(path) => {
            let file = File::open(path)?;
            // SAFETY: This command only reads the input; the file must not be modified
            // while the memory map is alive.
            let lines = unsafe { read_lines_from_file(file)? };
            collect_stats(&lines, sequence_type, each, threads)?
        }
    };
    let kind = sequence_kind(&summary, sequence_type);

    match (format, each) {
        (StatsFormat::Text, false) => write_summary_text(&summary, kind),
        (StatsFormat::Text, true) => write_each_text(&summary.records, kind),
        (StatsFormat::Json, false) => write_summary_json(&summary, kind),
        (StatsFormat::Json, true) => write_each_json(&summary.records, kind),
    }
    Ok(())
}

fn collect_stats(
    lines: &LinesInFile,
    sequence_type: StatsSequenceType,
    each: bool,
    requested_threads: Option<NonZeroUsize>,
) -> io::Result<Summary> {
    let bytes = lines.as_bytes();
    let chunks = parallel::record_chunks(
        bytes,
        parallel::worker_count(bytes.len(), requested_threads),
    );
    let partials = parallel::process_chunks(&chunks, |chunk| {
        let mut collector = StatsCollector::new(sequence_type, each);
        for (line_index, line) in chunk
            .bytes
            .split_inclusive(|byte| *byte == b'\n')
            .enumerate()
        {
            collector.process_line(chunk.start_line + line_index, line)?;
        }
        collector.finish(false)
    });

    let mut summary = Summary::default();
    for partial in partials {
        merge_summary(&mut summary, partial?, each)?;
    }
    compute_n50(&mut summary);
    Ok(summary)
}

fn collect_stats_from_reader(
    reader: impl io::BufRead,
    sequence_type: StatsSequenceType,
    each: bool,
) -> io::Result<Summary> {
    let mut collector = StatsCollector::new(sequence_type, each);
    for_each_reader_line(reader, |line_number, line| {
        collector.process_line(line_number, line)
    })?;
    collector.finish(true)
}

impl StatsCollector {
    fn new(sequence_type: StatsSequenceType, each: bool) -> Self {
        Self {
            summary: Summary::default(),
            current: None,
            saw_record: false,
            sequence_type,
            each,
        }
    }

    fn process_line(&mut self, line_number: usize, raw_line: &[u8]) -> io::Result<()> {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            if let Some((record, counts)) = self.current.take() {
                add_record(&mut self.summary, record, counts, self.each)?;
            }
            let id = line[1..]
                .trim_ascii_start()
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            if id.is_empty() {
                return Err(stats_error(line_number, "record identifier is empty"));
            }
            let id = self.each.then(|| String::from_utf8_lossy(id).into_owned());
            self.current = Some((id, Counts::default()));
            self.saw_record = true;
            return Ok(());
        }

        let Some((_, counts)) = self.current.as_mut() else {
            if line.is_empty() {
                return Ok(());
            }
            return Err(stats_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        };

        let sequence = line;
        let mut line_gc = 0u64;
        let mut line_n = 0u64;
        for (index, byte) in sequence.iter().copied().enumerate() {
            if byte.is_ascii_whitespace() {
                return Err(stats_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", index + 1),
                ));
            }

            let (is_nucleic, valid) = match self.sequence_type {
                StatsSequenceType::Nucleotide => {
                    let is_nucleic = is_nucleic_acid(byte);
                    (is_nucleic, is_nucleic)
                }
                StatsSequenceType::Protein => (false, is_amino_acid(byte)),
                StatsSequenceType::Auto => {
                    let is_nucleic = is_nucleic_acid(byte);
                    (is_nucleic, is_nucleic || is_amino_acid(byte))
                }
            };
            if !valid {
                let kind = match self.sequence_type {
                    StatsSequenceType::Protein => "protein",
                    StatsSequenceType::Nucleotide | StatsSequenceType::Auto => "nucleotide",
                };
                return Err(stats_error(
                    line_number,
                    &format!("invalid {kind} symbol '{}'", char::from(byte)),
                ));
            }

            if matches!(byte, b'G' | b'g' | b'C' | b'c') {
                line_gc += 1;
            }
            if matches!(byte, b'N' | b'n') {
                line_n += 1;
            }
            if self.sequence_type != StatsSequenceType::Protein {
                if byte == b'T' || byte == b't' {
                    self.summary.saw_t = true;
                } else if byte == b'U' || byte == b'u' {
                    self.summary.saw_u = true;
                }
            }
            if self.sequence_type == StatsSequenceType::Auto && !is_nucleic {
                self.summary.saw_protein_only = true;
            }
        }
        let line_length = u64::try_from(sequence.len())
            .map_err(|_| stats_error(line_number, "sequence length overflow"))?;
        counts.length = checked_add(
            counts.length,
            line_length,
            line_number,
            "sequence length overflow",
        )?;
        counts.gc = checked_add(counts.gc, line_gc, line_number, "GC count overflow")?;
        counts.n = checked_add(counts.n, line_n, line_number, "N count overflow")?;
        Ok(())
    }

    fn finish(mut self, calculate_n50: bool) -> io::Result<Summary> {
        if let Some((record, counts)) = self.current {
            add_record(&mut self.summary, record, counts, self.each)?;
        }
        if !self.saw_record {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no FASTA records found",
            ));
        }
        if !self.each && calculate_n50 {
            compute_n50(&mut self.summary);
        }
        Ok(self.summary)
    }
}

fn add_record(
    summary: &mut Summary,
    id: Option<String>,
    counts: Counts,
    each: bool,
) -> io::Result<()> {
    summary.record_count = summary
        .record_count
        .checked_add(1)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "record count overflow"))?;
    summary.total_length = summary
        .total_length
        .checked_add(counts.length)
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "total sequence length overflow")
        })?;
    summary.total_gc = summary
        .total_gc
        .checked_add(counts.gc)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "GC count overflow"))?;
    summary.total_n = summary
        .total_n
        .checked_add(counts.n)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "N count overflow"))?;
    summary.min_length = if summary.record_count == 1 {
        counts.length
    } else {
        summary.min_length.min(counts.length)
    };
    summary.max_length = summary.max_length.max(counts.length);
    if each {
        summary.records.push(RecordStats {
            id: id.unwrap_or_default(),
            length: counts.length,
            gc: counts.gc,
            n: counts.n,
        });
    } else {
        summary.lengths.push(counts.length);
    }
    Ok(())
}

fn merge_summary(summary: &mut Summary, partial: Summary, each: bool) -> io::Result<()> {
    summary.record_count = summary
        .record_count
        .checked_add(partial.record_count)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "record count overflow"))?;
    summary.total_length = summary
        .total_length
        .checked_add(partial.total_length)
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "total sequence length overflow")
        })?;
    summary.total_gc = summary
        .total_gc
        .checked_add(partial.total_gc)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "GC count overflow"))?;
    summary.total_n = summary
        .total_n
        .checked_add(partial.total_n)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "N count overflow"))?;
    if summary.record_count == partial.record_count {
        summary.min_length = partial.min_length;
    } else {
        summary.min_length = summary.min_length.min(partial.min_length);
    }
    summary.max_length = summary.max_length.max(partial.max_length);
    summary.saw_t |= partial.saw_t;
    summary.saw_u |= partial.saw_u;
    summary.saw_protein_only |= partial.saw_protein_only;
    if each {
        summary.records.extend(partial.records);
    } else {
        summary.lengths.extend(partial.lengths);
    }
    Ok(())
}

fn compute_n50(summary: &mut Summary) {
    summary
        .lengths
        .sort_unstable_by(|left, right| right.cmp(left));
    let threshold = summary.total_length / 2 + summary.total_length % 2;
    let mut cumulative = 0u64;
    for &length in &summary.lengths {
        cumulative += length;
        if cumulative >= threshold {
            summary.n50 = length;
            break;
        }
    }
}

fn checked_add(value: u64, addition: u64, line_number: usize, message: &str) -> io::Result<u64> {
    value
        .checked_add(addition)
        .ok_or_else(|| stats_error(line_number, message))
}

fn stats_error(line_number: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("line {line_number}: {message}"),
    )
}

fn sequence_kind(summary: &Summary, sequence_type: StatsSequenceType) -> &'static str {
    match sequence_type {
        StatsSequenceType::Protein => "Protein",
        StatsSequenceType::Nucleotide => nucleotide_kind(summary),
        StatsSequenceType::Auto if summary.saw_protein_only => "Protein",
        StatsSequenceType::Auto => nucleotide_kind(summary),
    }
}

fn nucleotide_kind(summary: &Summary) -> &'static str {
    match (summary.saw_t, summary.saw_u) {
        (true, false) => "DNA",
        (false, true) => "RNA",
        (true, true) => "DNA/RNA mixed",
        (false, false) => "DNA/RNA ambiguous",
    }
}

fn percentage(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 * 100.0 / total as f64
    }
}

fn rounded_mean(summary: &Summary) -> u64 {
    let count = summary.record_count;
    let quotient = summary.total_length / count;
    let remainder = summary.total_length % count;
    if remainder >= count / 2 + count % 2 {
        quotient + 1
    } else {
        quotient
    }
}

fn is_protein(kind: &str) -> bool {
    kind == "Protein"
}

fn write_summary_text(summary: &Summary, kind: &str) {
    let gc = if is_protein(kind) {
        "n/a".to_owned()
    } else {
        format!("{:.2}%", percentage(summary.total_gc, summary.total_length))
    };
    let n = if is_protein(kind) {
        "n/a".to_owned()
    } else {
        format!("{:.2}%", percentage(summary.total_n, summary.total_length))
    };
    println!("sequences    {}", summary.record_count);
    println!("total_len    {}", with_grouping(summary.total_length));
    println!("min_len      {}", with_grouping(summary.min_length));
    println!("max_len      {}", with_grouping(summary.max_length));
    println!("mean_len     {}", with_grouping(rounded_mean(summary)));
    println!("N50          {}", with_grouping(summary.n50));
    println!("GC           {gc}");
    println!("N            {n}");
    println!("type         {kind}");
}

fn write_each_text(records: &[RecordStats], kind: &str) {
    let id_width = records
        .iter()
        .map(|record| record.id.len())
        .max()
        .unwrap_or(2)
        .max(2);
    println!(
        "{:<id_width$}  {:>10}  {:>8}  {:>8}",
        "id", "length", "gc", "n"
    );
    for record in records {
        let (gc, n) = if is_protein(kind) {
            ("n/a".to_owned(), "n/a".to_owned())
        } else {
            (
                format!("{:.2}%", percentage(record.gc, record.length)),
                format!("{:.2}%", percentage(record.n, record.length)),
            )
        };
        println!(
            "{:<id_width$}  {:>10}  {:>8}  {:>8}",
            record.id, record.length, gc, n
        );
    }
}

fn write_summary_json(summary: &Summary, kind: &str) {
    let gc = if is_protein(kind) {
        "null".to_owned()
    } else {
        format!("{:.6}", percentage(summary.total_gc, summary.total_length))
    };
    let n = if is_protein(kind) {
        "null".to_owned()
    } else {
        format!("{:.6}", percentage(summary.total_n, summary.total_length))
    };
    let mean = summary.total_length as f64 / summary.record_count as f64;
    println!(
        "{{\n  \"sequences\": {},\n  \"total_len\": {},\n  \"min_len\": {},\n  \"max_len\": {},\n  \"mean_len\": {mean},\n  \"n50\": {},\n  \"gc_percent\": {gc},\n  \"n_percent\": {n},\n  \"type\": {}\n}}",
        summary.record_count,
        summary.total_length,
        summary.min_length,
        summary.max_length,
        summary.n50,
        json_string(kind),
    );
}

fn write_each_json(records: &[RecordStats], kind: &str) {
    println!("[");
    for (index, record) in records.iter().enumerate() {
        let (gc, n) = if is_protein(kind) {
            ("null".to_owned(), "null".to_owned())
        } else {
            (
                format!("{:.6}", percentage(record.gc, record.length)),
                format!("{:.6}", percentage(record.n, record.length)),
            )
        };
        let comma = if index + 1 == records.len() { "" } else { "," };
        println!(
            "  {{\"id\": {}, \"length\": {}, \"gc_percent\": {gc}, \"n_percent\": {n}}}{comma}",
            json_string(&record.id),
            record.length,
        );
    }
    println!("]");
}

fn json_string(value: &str) -> String {
    let mut output = String::with_capacity(value.len() + 2);
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(output, "\\u{:04x}", character as u32);
            }
            character => output.push(character),
        }
    }
    output.push('"');
    output
}

fn with_grouping(value: u64) -> String {
    let digits = value.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, character) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(character);
    }
    grouped
}
