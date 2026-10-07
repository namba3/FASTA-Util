mod fasta_index;
mod filter;
mod format;
mod get;
mod grep;
mod locate;
mod output;
mod revcomp;
mod stats;
mod validate;

use clap::{Parser, Subcommand, ValueEnum};
#[cfg(test)]
use crossbeam::channel::Receiver;
use fasta_util::{
    LinesInFile, is_amino_acid, is_nucleic_acid, read_lines_from_file, read_lines_from_stdin,
};
use output::TemporaryOutput;
use std::{
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    #[command(subcommand)]
    sub: SubCommand,
}

#[derive(Subcommand)]
enum SubCommand {
    #[command(about = "Count the total length of the sequence")]
    Len(LenArgs),
    #[command(about = "Validate FASTA structure and sequence symbols")]
    Validate(ValidateArgs),
    #[command(about = "Create a FASTA .fai index")]
    Index(IndexArgs),
    #[command(about = "Summarize FASTA sequence statistics")]
    Stats(StatsArgs),
    #[command(about = "Get FASTA records, ID regions, or global ranges")]
    Get(GetArgs),
    #[command(about = "Filter FASTA records by sequence properties")]
    Filter(FilterArgs),
    #[command(about = "Reverse-complement nucleotide sequences")]
    Revcomp(RevcompArgs),
    #[command(about = "Search FASTA record headers")]
    Grep(GrepArgs),
    #[command(about = "Locate IUPAC nucleotide motifs in FASTA sequences")]
    Locate(LocateArgs),
    #[command(about = "Reformat FASTA records and sequence lines")]
    Format(FormatArgs),
}

#[derive(Parser)]
struct ValidateArgs {
    /// FASTA file to validate
    input: PathBuf,

    #[arg(
        long,
        value_enum,
        default_value_t = SequenceType::Nucleotide,
        help = "Sequence alphabet to validate (nucleotide or protein)"
    )]
    sequence_type: SequenceType,
}

#[derive(Parser)]
struct IndexArgs {
    /// FASTA file to index; writes <input>.fai
    input: PathBuf,
}

#[derive(Parser)]
struct StatsArgs {
    /// FASTA file to summarize
    input: PathBuf,

    /// Report statistics for each record instead of the whole file
    #[arg(long)]
    each: bool,

    /// Output format
    #[arg(long, value_enum, default_value_t = StatsFormat::Text)]
    format: StatsFormat,

    /// Sequence alphabet; auto selects protein if a protein-only symbol appears
    #[arg(long, value_enum, default_value_t = StatsSequenceType::Auto)]
    sequence_type: StatsSequenceType,
}

#[derive(Parser)]
struct GetArgs {
    /// FASTA file to read
    input: PathBuf,

    /// Record IDs, ID regions, or global ranges (coordinates are 1-based and inclusive)
    ids: Vec<String>,

    /// Read IDs or regions from a newline-delimited file
    #[arg(long = "ids", alias = "ids-file", conflicts_with = "ids")]
    ids_file: Option<PathBuf>,

    /// Use this FASTA .fai index; defaults to <input>.fai when present
    #[arg(long)]
    fai_index: Option<PathBuf>,

    /// Write output to a file instead of standard output
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Sequence alphabet to validate (nucleotide or protein)
    #[arg(long, value_enum, default_value_t = SequenceType::Nucleotide)]
    sequence_type: SequenceType,

    /// Number of sequence characters per output line
    #[arg(long, default_value_t = 60, value_parser = parse_positive_line_width)]
    chars_per_line: usize,
}

#[derive(Parser)]
struct FilterArgs {
    /// FASTA file to filter
    input: PathBuf,

    /// Keep records with at least this many sequence symbols
    #[arg(long)]
    min_len: Option<u64>,

    /// Keep records with at most this many sequence symbols
    #[arg(long)]
    max_len: Option<u64>,

    /// Keep records with at least this GC fraction (0.0 to 1.0)
    #[arg(long, value_parser = parse_fraction)]
    min_gc: Option<f64>,

    /// Keep records with at most this GC fraction (0.0 to 1.0)
    #[arg(long, value_parser = parse_fraction)]
    max_gc: Option<f64>,

    /// Keep records with at most this N fraction (0.0 to 1.0)
    #[arg(long, value_parser = parse_fraction)]
    max_n: Option<f64>,

    /// Sequence alphabet to validate (nucleotide or protein)
    #[arg(long, value_enum, default_value_t = StatsSequenceType::Auto)]
    sequence_type: StatsSequenceType,

    /// Write output to a file instead of standard output
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Parser)]
struct RevcompArgs {
    /// FASTA file containing nucleotide sequences
    input: PathBuf,

    /// Write output to a file instead of standard output
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Number of sequence characters per output line
    #[arg(long, default_value_t = 60, value_parser = parse_positive_line_width)]
    chars_per_line: usize,
}

#[derive(Parser)]
struct GrepArgs {
    /// FASTA file to search
    input: PathBuf,

    /// Literal text to search for in record headers
    pattern: String,

    /// Match ASCII letters without regard to case
    #[arg(short, long)]
    ignore_case: bool,

    /// Keep records whose headers do not match the pattern
    #[arg(short = 'v', long)]
    invert_match: bool,

    /// Write matching records to a file instead of standard output
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Parser)]
struct LocateArgs {
    /// FASTA file to search
    input: PathBuf,

    /// IUPAC nucleotide motif to locate
    pattern: String,

    /// Maximum number of mismatching positions (default: 0)
    #[arg(long, default_value_t = 0)]
    max_mismatch: usize,

    /// Write tab-separated matches to a file instead of standard output
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Parser)]
struct FormatArgs {
    /// FASTA file to format
    input: PathBuf,

    /// Sequence characters per line; 0 writes one sequence line per record
    #[arg(long, default_value_t = 60)]
    width: usize,

    /// Convert sequence symbols to uppercase
    #[arg(long, conflicts_with = "lowercase")]
    uppercase: bool,

    /// Convert sequence symbols to lowercase
    #[arg(long, conflicts_with = "uppercase")]
    lowercase: bool,

    /// Remove gap symbols (-) from sequences
    #[arg(long)]
    remove_gaps: bool,

    /// Trim surrounding ASCII whitespace from header text
    #[arg(long)]
    trim_header: bool,

    /// Write output to a file instead of standard output
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum StatsFormat {
    #[default]
    Text,
    Json,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum StatsSequenceType {
    #[default]
    Auto,
    Nucleotide,
    Protein,
}

#[derive(Parser)]
struct LenArgs {
    #[arg(
        short,
        long,
        help = "Specify input file\nIf omitted, read from standard input"
    )]
    input: Option<PathBuf>,

    #[arg(
        long,
        value_enum,
        default_value_t = SequenceType::Nucleotide,
        help = "Sequence alphabet to validate (nucleotide or protein)"
    )]
    sequence_type: SequenceType,
}

struct GlobalRangeArgs {
    input: PathBuf,
    sequence_type: SequenceType,
    output: Option<PathBuf>,
    fai_index: Option<PathBuf>,
    start: usize,
    end_exclusive: usize,
    chars_per_line: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum SequenceType {
    #[default]
    Nucleotide,
    Protein,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.sub {
        SubCommand::Len(args) => len(args)?,
        SubCommand::Validate(args) => {
            if !validate::run(&args.input, args.sequence_type)? {
                std::process::exit(1);
            }
        }
        SubCommand::Index(args) => {
            let index_path = fasta_index::index_path(&args.input);
            ensure_distinct_input_output(&args.input, &index_path)?;
            let records = fasta_index::create_index(&args.input, &index_path)?;
            println!("Indexed {records} records: {}", index_path.display());
        }
        SubCommand::Stats(args) => {
            stats::run(&args.input, args.each, args.format, args.sequence_type)?
        }
        SubCommand::Get(args) => get::run(args)?,
        SubCommand::Filter(args) => filter::run(args)?,
        SubCommand::Revcomp(args) => revcomp::run(args)?,
        SubCommand::Grep(args) => grep::run(args)?,
        SubCommand::Locate(args) => locate::run(args)?,
        SubCommand::Format(args) => format::run(args)?,
    }

    Ok(())
}

fn parse_positive_line_width(value: &str) -> Result<usize, String> {
    let width = value
        .parse::<usize>()
        .map_err(|error| format!("invalid characters-per-line value: {error}"))?;
    if width == 0 {
        return Err("characters per line must be greater than zero".to_owned());
    }
    Ok(width)
}

fn parse_fraction(value: &str) -> Result<f64, String> {
    let fraction = value
        .parse::<f64>()
        .map_err(|error| format!("invalid fraction: {error}"))?;
    if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
        return Err("fraction must be between 0.0 and 1.0".to_owned());
    }
    Ok(fraction)
}

fn len(args: LenArgs) -> Result<(), Box<dyn std::error::Error>> {
    let len = match args.input {
        Some(input) => {
            let input = std::fs::OpenOptions::new().read(true).open(input)?;
            // SAFETY: Input files must remain unchanged for the duration of this command;
            // this command only reads the file and never modifies it.
            let lines = unsafe { read_lines_from_file(input)? };
            count_sequence_bases_from_file(&lines, args.sequence_type)?
        }
        None => count_sequence_bases_for(read_lines_from_stdin(), args.sequence_type)?,
    };

    let stdout = io::stdout();
    write_length(&mut stdout.lock(), len)?;

    Ok(())
}

fn write_length(writer: &mut impl Write, length: u64) -> io::Result<()> {
    writeln!(writer, "{length}")
}

#[cfg(test)]
fn validated_sequence(line: &[u8]) -> io::Result<&[u8]> {
    validated_sequence_for(line, SequenceType::Nucleotide)
}

pub(crate) fn validated_sequence_for(
    line: &[u8],
    sequence_type: SequenceType,
) -> io::Result<&[u8]> {
    let sequence = line.trim_ascii_start().trim_ascii_end();
    if let Some(byte) = sequence
        .iter()
        .find(|byte| !is_sequence_symbol(**byte, sequence_type))
    {
        return Err(invalid_sequence_symbol(*byte, sequence_type));
    }
    Ok(sequence)
}

fn is_sequence_symbol(byte: u8, sequence_type: SequenceType) -> bool {
    match sequence_type {
        SequenceType::Nucleotide => is_nucleic_acid(byte),
        SequenceType::Protein => is_amino_acid(byte),
    }
}

fn with_line_context(line_number: usize, error: io::Error) -> io::Error {
    io::Error::new(error.kind(), format!("line {line_number}: {error}"))
}

fn sequence_length_overflow(line_number: usize) -> io::Error {
    with_line_context(
        line_number,
        io::Error::new(io::ErrorKind::InvalidData, "sequence length overflow"),
    )
}

#[cfg(test)]
fn count_sequence_bases<T, I>(iter: I) -> io::Result<u64>
where
    T: AsRef<[u8]>,
    I: IntoIterator<Item = Result<T, io::Error>>,
{
    count_sequence_bases_for(iter, SequenceType::Nucleotide)
}

fn count_sequence_bases_for<T, I>(iter: I, sequence_type: SequenceType) -> io::Result<u64>
where
    T: AsRef<[u8]>,
    I: IntoIterator<Item = Result<T, io::Error>>,
{
    let mut count = 0u64;
    for (line_index, line) in iter.into_iter().enumerate() {
        let line_number = line_index + 1;
        let line = line.map_err(|error| with_line_context(line_number, error))?;
        count_sequence_line(line_number, line.as_ref(), &mut count, sequence_type)?;
    }
    Ok(count)
}

fn count_sequence_bases_from_file(
    lines: &LinesInFile,
    sequence_type: SequenceType,
) -> io::Result<u64> {
    let mut count = 0u64;
    lines.try_for_each_line(|line_number, line| {
        count_sequence_line(line_number, line, &mut count, sequence_type)
    })?;
    Ok(count)
}

fn count_sequence_line(
    line_number: usize,
    line: &[u8],
    count: &mut u64,
    sequence_type: SequenceType,
) -> io::Result<()> {
    if line.first() == Some(&b'>') {
        return Ok(());
    }

    let sequence = validated_sequence_for(line, sequence_type)
        .map_err(|error| with_line_context(line_number, error))?;
    if !sequence.is_empty() {
        let length =
            u64::try_from(sequence.len()).map_err(|_| sequence_length_overflow(line_number))?;
        *count = (*count)
            .checked_add(length)
            .ok_or_else(|| sequence_length_overflow(line_number))?;
    }
    Ok(())
}

fn invalid_sequence_symbol(byte: u8, sequence_type: SequenceType) -> io::Error {
    let message = match sequence_type {
        SequenceType::Nucleotide => "invalid nucleic acid",
        SequenceType::Protein => "invalid protein symbol",
    };
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("{message}: {:?} (0x{byte:02x})", char::from(byte)),
    )
}

fn strip_line_ending(line: &[u8]) -> &[u8] {
    match line.strip_suffix(b"\n") {
        Some(line) => line.strip_suffix(b"\r").unwrap_or(line),
        None => line,
    }
}

#[derive(Debug, PartialEq, Eq)]
#[cfg(test)]
struct SequenceRange {
    start: usize,
    end_exclusive: Option<usize>,
}

fn ensure_distinct_input_output(input_path: &Path, output_path: &Path) -> io::Result<()> {
    let _output_metadata = match std::fs::metadata(output_path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error),
    };

    let same_file = std::fs::canonicalize(input_path)? == std::fs::canonicalize(output_path)?;
    #[cfg(unix)]
    let same_file = {
        use std::os::unix::fs::MetadataExt;
        let input_metadata = std::fs::metadata(input_path)?;
        same_file
            || (input_metadata.dev() == _output_metadata.dev()
                && input_metadata.ino() == _output_metadata.ino())
    };

    if same_file {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input and output refer to the same file",
        ));
    }

    Ok(())
}

fn write_global_range(args: GlobalRangeArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(output) = &args.output {
        ensure_distinct_input_output(&args.input, output)?;
    }
    if let (Some(index), Some(output)) = (&args.fai_index, &args.output) {
        ensure_distinct_input_output(index, output)?;
    }

    if let Some(index_path) = &args.fai_index {
        let mut temporary_output = args
            .output
            .as_deref()
            .map(TemporaryOutput::create)
            .transpose()?;
        let output: Box<dyn Write> = match temporary_output.as_mut() {
            Some(temporary_output) => Box::new(temporary_output.take_file()?),
            None => Box::new(std::io::stdout().lock()),
        };
        let mut writer = BufWriter::new(output);
        fasta_index::write_slice(
            &args.input,
            index_path,
            args.start,
            Some(args.end_exclusive),
            args.chars_per_line,
            args.sequence_type,
            &mut writer,
        )?;
        writer.flush()?;
        drop(writer);
        if let Some(temporary_output) = &mut temporary_output {
            temporary_output.commit()?;
        }
        return Ok(());
    }

    let writer_options = WriterOptions {
        chars_per_line: args.chars_per_line,
        start: args.start,
        end_exclusive: Some(args.end_exclusive),
        sequence_type: args.sequence_type,
    };

    let input = std::fs::OpenOptions::new().read(true).open(&args.input)?;
    // SAFETY: Input files must remain unchanged for the duration of this command;
    // this command only reads the file and never modifies it.
    let file_lines = unsafe { read_lines_from_file(input)? };

    let mut temporary_output = args
        .output
        .as_deref()
        .map(TemporaryOutput::create)
        .transpose()?;
    let output: Box<dyn Write> = match temporary_output.as_mut() {
        Some(temporary_output) => Box::new(temporary_output.take_file()?),
        None => Box::new(std::io::stdout().lock()),
    };

    let mut writer = Writer::new(output, writer_options);

    writer.run_file(&file_lines)?;
    drop(writer);
    if let Some(temporary_output) = &mut temporary_output {
        temporary_output.commit()?;
    }

    Ok(())
}

struct WriterOptions {
    chars_per_line: usize,
    start: usize,
    end_exclusive: Option<usize>,
    sequence_type: SequenceType,
}
struct Writer<T: std::io::Write> {
    inner: BufWriter<T>,
    options: WriterOptions,
    count: usize,
    written: usize,
}
impl<T: std::io::Write> Writer<T> {
    fn new(inner: T, options: WriterOptions) -> Self {
        Self {
            inner: BufWriter::new(inner),
            options,
            count: 0,
            written: 0,
        }
    }
    #[cfg(test)]
    fn run<Buf: AsRef<[u8]>>(
        &mut self,
        rx: Receiver<Result<Buf, io::Error>>,
    ) -> Result<(), std::io::Error> {
        let mut line_number = 0usize;

        while let Ok(line) = rx.recv() {
            line_number += 1;
            let line = line.map_err(|error| with_line_context(line_number, error))?;
            if !self.process_line(line_number, line.as_ref())? {
                break;
            }
        }

        self.inner.flush()
    }

    fn run_file(&mut self, lines: &LinesInFile) -> Result<(), io::Error> {
        lines.try_for_each_line_while(|line_number, line| self.process_line(line_number, line))?;
        self.inner.flush()
    }

    fn process_line(&mut self, line_number: usize, line: &[u8]) -> Result<bool, io::Error> {
        let buf = strip_line_ending(line);
        let writer = &mut self.inner;
        let chars_per_line = self.options.chars_per_line;

        if let Some(b'>') = buf.first() {
            if self.written > 0 && !self.written.is_multiple_of(chars_per_line) {
                writer.write_all(b"\n")?;
            }
            writer.write_all(buf)?;
            writer.write_all(b"\n")?;
            return Ok(true);
        }

        let buf = validated_sequence_for(buf, self.options.sequence_type)
            .map_err(|error| with_line_context(line_number, error))?;
        if buf.is_empty() {
            return Ok(true);
        }

        let line_end = self
            .count
            .checked_add(buf.len())
            .ok_or_else(|| sequence_length_overflow(line_number))?;
        let start_in_line = self.options.start.saturating_sub(self.count);
        let end_in_line = self
            .options
            .end_exclusive
            .map(|end| end.saturating_sub(self.count).min(buf.len()))
            .unwrap_or(buf.len());

        if start_in_line >= end_in_line {
            self.count = line_end;
            return Ok(!self
                .options
                .end_exclusive
                .is_some_and(|end| end <= self.count));
        }

        let mut bases = &buf[start_in_line..end_in_line];
        let line_written = if start_in_line == 0 {
            self.written % chars_per_line
        } else {
            0
        };
        let mut line_remain = chars_per_line - line_written;
        let written_end = self
            .written
            .checked_add(bases.len())
            .ok_or_else(|| sequence_length_overflow(line_number))?;

        while line_remain <= bases.len() {
            writer.write_all(&bases[..line_remain])?;
            writer.write_all(b"\n")?;
            bases = &bases[line_remain..];
            line_remain = chars_per_line;
        }

        writer.write_all(bases)?;
        self.count = line_end;
        self.written = written_end;
        if self
            .options
            .end_exclusive
            .is_some_and(|end| end <= self.count)
        {
            if self.written > 0 && !self.written.is_multiple_of(chars_per_line) {
                writer.write_all(b"\n")?;
            }
            return Ok(false);
        }

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Args, SequenceRange, SequenceType, Writer, WriterOptions, count_sequence_bases,
        count_sequence_line, strip_line_ending, validated_sequence, write_length,
    };
    use clap::Parser;
    use crossbeam::channel::unbounded;

    fn write_fasta(lines: &[&[u8]], options: WriterOptions) -> Vec<u8> {
        write_fasta_result(lines, options).unwrap()
    }

    fn write_fasta_result(
        lines: &[&[u8]],
        options: WriterOptions,
    ) -> Result<Vec<u8>, std::io::Error> {
        let (tx, rx) = unbounded();
        for line in lines {
            tx.send(Ok(*line)).unwrap();
        }
        drop(tx);

        let mut writer = Writer::new(Vec::new(), options);
        writer.run(rx)?;
        Ok(writer.inner.into_inner().unwrap())
    }

    fn write_fasta_for_range(
        lines: &[&[u8]],
        range: SequenceRange,
        chars_per_line: usize,
    ) -> Result<Vec<u8>, std::io::Error> {
        Ok(write_fasta(
            lines,
            WriterOptions {
                chars_per_line,
                start: range.start,
                end_exclusive: range.end_exclusive,
                sequence_type: SequenceType::Nucleotide,
            },
        ))
    }

    fn options(
        chars_per_line: usize,
        start: Option<usize>,
        end_exclusive: Option<usize>,
    ) -> WriterOptions {
        WriterOptions {
            chars_per_line,
            start: start.unwrap_or(0),
            end_exclusive,
            sequence_type: SequenceType::Nucleotide,
        }
    }

    #[test]
    fn get_requires_a_positive_line_width() {
        assert!(
            Args::try_parse_from([
                "fasta-util",
                "get",
                "input.fa",
                "1-2",
                "--chars-per-line",
                "0"
            ])
            .is_err()
        );
        assert!(
            Args::try_parse_from([
                "fasta-util",
                "get",
                "input.fa",
                "1-2",
                "--chars-per-line",
                "1"
            ])
            .is_ok()
        );
        assert!(
            Args::try_parse_from(["fasta-util", "get", "input.fa", "--range", "1..2"]).is_err()
        );
        assert!(Args::try_parse_from(["fasta-util", "slice"]).is_err());
    }

    #[test]
    fn len_counts_sequence_lines_and_ignores_headers_and_blank_lines() {
        let lines: [&[u8]; 5] = [
            b">record 1\n",
            b"ACGT\n",
            b" \t\n",
            b"NU-\r\n",
            b">record 2\n",
        ];

        assert_eq!(
            count_sequence_bases(lines.into_iter().map(Ok::<_, std::io::Error>)).unwrap(),
            7
        );
    }

    #[test]
    fn len_reports_sequence_count_overflow_with_line_context() {
        let mut count = u64::MAX;

        let error = count_sequence_line(8, b"A", &mut count, SequenceType::Nucleotide).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(error.to_string(), "line 8: sequence length overflow");
        assert_eq!(count, u64::MAX);
    }

    #[test]
    fn validated_sequence_trims_whitespace_and_accepts_empty_lines() {
        assert_eq!(validated_sequence(b" \tACGT-\r\n").unwrap(), b"ACGT-");
        assert!(validated_sequence(b" \t\r\n").unwrap().is_empty());
    }

    #[test]
    fn validated_sequence_accepts_and_preserves_lowercase_symbols() {
        assert_eq!(
            validated_sequence(b"acgtnuk-symwrbdhv").unwrap(),
            b"acgtnuk-symwrbdhv"
        );
    }

    #[test]
    fn len_counts_lowercase_soft_masked_symbols() {
        assert_eq!(
            count_sequence_bases([Ok::<_, std::io::Error>(b"aCgTn".as_slice())]).unwrap(),
            5
        );
    }

    #[test]
    fn len_writes_the_result_with_a_newline() {
        let mut output = Vec::new();

        write_length(&mut output, 42).unwrap();

        assert_eq!(output, b"42\n");
    }

    #[test]
    fn len_propagates_output_write_errors() {
        struct FailingWriter;

        impl std::io::Write for FailingWriter {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::new(
                    std::io::ErrorKind::BrokenPipe,
                    "output closed",
                ))
            }

            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }

        let error = write_length(&mut FailingWriter, 42).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }

    #[test]
    fn validated_sequence_reports_invalid_symbols_after_trimming() {
        let error = validated_sequence(b" \tACX\r\n").unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("invalid nucleic acid"));
    }

    #[test]
    fn len_rejects_invalid_sequence_symbols() {
        let error = count_sequence_bases([Ok::<_, std::io::Error>(b"ACX".as_slice())]).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("invalid nucleic acid"));
    }

    #[test]
    fn len_propagates_input_read_errors() {
        let error = count_sequence_bases([Err::<&[u8], _>(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "read failed",
        ))])
        .unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::PermissionDenied);
    }

    #[test]
    fn slice_wraps_sequence_at_requested_line_width() {
        let output = write_fasta(&[b">record\n", b"ACGT\n", b"NU\n"], options(3, None, None));

        assert_eq!(output, b">record\nACG\nTNU\n");
    }

    #[test]
    fn slice_preserves_lowercase_soft_masking() {
        let output = write_fasta(&[b">record\n", b"aCgTn\n"], options(10, None, None));

        assert_eq!(output, b">record\naCgTn");
    }

    #[test]
    fn slice_uses_inclusive_indices_across_sequence_lines() {
        let output = write_fasta(
            &[b">record\n", b"ACGT\n", b"NU-\n"],
            options(2, Some(2), Some(5)),
        );

        assert_eq!(output, b">record\nGT\nN\n");
    }

    #[test]
    fn slice_can_start_at_an_offset_without_an_end() {
        let output = write_fasta(
            &[b">record\n", b"ACGT\n", b"NU\n"],
            options(10, Some(3), None),
        );

        assert_eq!(output, b">record\nTNU");
    }

    #[test]
    fn slice_starting_at_a_line_boundary_uses_the_next_line() {
        let output = write_fasta(
            &[b">record\n", b"ACGT\n", b"NU\n"],
            options(10, Some(4), Some(6)),
        );

        assert_eq!(output, b">record\nNU\n");
    }

    #[test]
    fn empty_slice_range_writes_no_sequence_bases() {
        let output = write_fasta_for_range(
            &[b">record\n", b"ACGT\n"],
            SequenceRange {
                start: 0,
                end_exclusive: Some(0),
            },
            10,
        )
        .unwrap();

        assert_eq!(output, b">record\n");
    }

    #[test]
    fn exclusive_and_inclusive_range_text_selects_expected_bases() {
        let input = &[b">record\n".as_slice(), b"ACGTNU\n".as_slice()];
        let exclusive = write_fasta_for_range(
            input,
            SequenceRange {
                start: 2,
                end_exclusive: Some(4),
            },
            10,
        )
        .unwrap();
        let inclusive = write_fasta_for_range(
            input,
            SequenceRange {
                start: 2,
                end_exclusive: Some(5),
            },
            10,
        )
        .unwrap();

        assert_eq!(exclusive, b">record\nGT\n");
        assert_eq!(inclusive, b">record\nGTN\n");
    }

    #[test]
    fn slice_separates_headers_from_a_partial_sequence_line() {
        let output = write_fasta(
            &[b">first\n", b"AC\n", b">second\n", b"GT\n"],
            options(10, None, None),
        );

        assert_eq!(output, b">first\nAC\n>second\nGT");
    }

    #[test]
    fn slice_terminates_headers_without_a_line_ending() {
        let output = write_fasta(&[b">record", b"ACGT"], options(10, None, None));

        assert_eq!(output, b">record\nACGT");
    }

    #[test]
    fn slice_normalizes_lf_and_crlf_headers_to_the_same_output() {
        let file_style = write_fasta(&[b">record\r\n", b"ACGT\r\n"], options(10, None, None));
        let stdin_style = write_fasta(&[b">record", b"ACGT"], options(10, None, None));

        assert_eq!(file_style, b">record\nACGT");
        assert_eq!(file_style, stdin_style);
    }

    #[test]
    fn strip_line_ending_removes_lf_and_crlf_but_preserves_unterminated_cr() {
        assert_eq!(strip_line_ending(b"line\n"), b"line");
        assert_eq!(strip_line_ending(b"line\r\n"), b"line");
        assert_eq!(strip_line_ending(b"line\r"), b"line\r");
    }

    #[test]
    fn slice_rejects_invalid_sequence_symbols() {
        let error = write_fasta_result(&[b"ACX\n"], options(10, None, None)).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("invalid nucleic acid"));
    }

    #[test]
    fn slice_reports_sequence_position_overflow_without_panicking() {
        let mut writer = Writer::new(Vec::new(), options(10, None, None));
        writer.count = usize::MAX;

        let error = writer.process_line(4, b"A").unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(error.to_string(), "line 4: sequence length overflow");
    }

    #[test]
    fn slice_reports_written_base_overflow_without_panicking() {
        let mut writer = Writer::new(Vec::new(), options(10, None, None));
        writer.written = usize::MAX;

        let error = writer.process_line(5, b"A").unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(error.to_string(), "line 5: sequence length overflow");
    }

    #[test]
    fn slice_propagates_input_read_errors() {
        let (tx, rx) = unbounded();
        let input_error: Result<&[u8], std::io::Error> = Err(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "read failed",
        ));
        tx.send(input_error).unwrap();
        drop(tx);

        let mut writer = Writer::new(Vec::new(), options(10, None, None));
        let error = writer.run(rx).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::PermissionDenied);
    }
}
