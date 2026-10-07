mod fasta_index;
mod output;

use clap::{Parser, Subcommand, ValueEnum};
use crossbeam::channel::{Receiver, bounded};
use fasta_util::{
    LinesInFile, is_amino_acid, is_nucleic_acid, read_lines_from_file, read_lines_from_stdin,
};
use output::TemporaryOutput;
use std::{
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

const LINE_CHANNEL_CAPACITY: usize = 32;

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
    #[command(about = "Cut out a part of the sequence")]
    Slice(SliceArgs),
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

#[derive(Parser)]
struct SliceArgs {
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

    #[arg(
        short,
        long,
        help = "Specify output file\nIf omitted, write to standard output"
    )]
    output: Option<PathBuf>,

    #[arg(
        long,
        help = "Use a matching FASTA .fai index to read only the selected sequence region"
    )]
    fai_index: Option<PathBuf>,

    #[arg(
        long,
        default_value = "..",
        help = "Specify slice range\nexamples:\n\t2..10\tmeans [2,10)\n\t2..=10\tmeans [2,10]\n\t..10\tmeans [0,10)\n\t2..\tmeans [2,∞)\n\t..\tmeans [0,∞)\n"
    )]
    range: String,

    #[arg(
        long,
        default_value_t = 60,
        value_parser = parse_positive_line_width,
        help = "Specify the number of characters per line when exporting a sequence"
    )]
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
        SubCommand::Slice(args) => slice(args)?,
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
struct SequenceRange {
    start: usize,
    end_exclusive: Option<usize>,
}

fn invalid_range(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
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

fn parse_slice_range(range: &str) -> Result<SequenceRange, io::Error> {
    let parse_index = |value: &str| {
        value
            .parse::<usize>()
            .map_err(|error| invalid_range(format!("invalid range index: {error}")))
    };
    let (start, end) = range
        .split_once("..")
        .ok_or_else(|| invalid_range("range must contain `..`"))?;
    let start = if start.is_empty() {
        0
    } else {
        parse_index(start)?
    };
    let end_exclusive = if end.is_empty() {
        None
    } else if let Some(inclusive_end) = end.strip_prefix('=') {
        let inclusive_end = parse_index(inclusive_end)?;
        if inclusive_end < start {
            return Err(invalid_range("range end precedes range start"));
        }
        // `usize::MAX + 1` cannot be represented, and is equivalent to an open end
        // because sequence offsets cannot exceed the addressable slice length.
        inclusive_end.checked_add(1)
    } else {
        Some(parse_index(end)?)
    };

    if end_exclusive.is_some_and(|end| end < start) {
        return Err(invalid_range("range end precedes range start"));
    }

    Ok(SequenceRange {
        start,
        end_exclusive,
    })
}

fn slice(args: SliceArgs) -> Result<(), Box<dyn std::error::Error>> {
    if args.fai_index.is_some() && args.input.is_none() {
        return Err(invalid_range("--fai-index requires --input").into());
    }
    if let (Some(input), Some(output)) = (&args.input, &args.output) {
        ensure_distinct_input_output(input, output)?;
    }
    if let (Some(index), Some(output)) = (&args.fai_index, &args.output) {
        ensure_distinct_input_output(index, output)?;
    }

    let range = parse_slice_range(&args.range)?;
    if let (Some(index_path), Some(input_path)) = (&args.fai_index, &args.input) {
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
            input_path,
            index_path,
            range.start,
            range.end_exclusive,
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
        start: range.start,
        end_exclusive: range.end_exclusive,
        sequence_type: args.sequence_type,
    };

    let file_lines = match args.input {
        Some(input) => {
            let input = std::fs::OpenOptions::new().read(true).open(input)?;
            // SAFETY: Input files must remain unchanged for the duration of this command;
            // this command only reads the file and never modifies it.
            Some(unsafe { read_lines_from_file(input)? })
        }
        None => None,
    };

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

    let write_result = match file_lines {
        Some(lines) => writer.run_file(&lines),
        None => {
            let (tx, rx) = bounded(LINE_CHANNEL_CAPACITY);
            let hndl = std::thread::spawn(move || -> Result<(), std::io::Error> {
                let lines = read_lines_from_stdin();
                for line in lines {
                    if tx.send(line).is_err() {
                        break;
                    }
                }

                Ok(())
            });

            let write_result = writer.run(rx);
            let read_result = hndl
                .join()
                .unwrap_or_else(|_| Err(io::Error::other("input reader thread panicked")));
            write_result?;
            read_result?;
            Ok(())
        }
    };
    write_result?;
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
        count_sequence_line, parse_slice_range, strip_line_ending, validated_sequence,
        write_length,
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
        range: &str,
        chars_per_line: usize,
    ) -> Result<Vec<u8>, std::io::Error> {
        let range = parse_slice_range(range)?;
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
    fn slice_requires_a_positive_line_width() {
        assert!(Args::try_parse_from(["fasta-util", "slice", "--chars-per-line", "0"]).is_err());
        assert!(Args::try_parse_from(["fasta-util", "slice", "--chars-per-line", "1"]).is_ok());
    }

    #[test]
    fn parses_open_and_bounded_range_forms() {
        let cases = [
            (
                "..",
                SequenceRange {
                    start: 0,
                    end_exclusive: None,
                },
            ),
            (
                "2..",
                SequenceRange {
                    start: 2,
                    end_exclusive: None,
                },
            ),
            (
                "..10",
                SequenceRange {
                    start: 0,
                    end_exclusive: Some(10),
                },
            ),
            (
                "2..10",
                SequenceRange {
                    start: 2,
                    end_exclusive: Some(10),
                },
            ),
            (
                "2..=10",
                SequenceRange {
                    start: 2,
                    end_exclusive: Some(11),
                },
            ),
        ];

        for (input, expected) in cases {
            assert_eq!(parse_slice_range(input).unwrap(), expected, "range {input}");
        }
    }

    #[test]
    fn accepts_empty_ranges_without_underflow() {
        let cases = [
            (
                "..0",
                SequenceRange {
                    start: 0,
                    end_exclusive: Some(0),
                },
            ),
            (
                "2..2",
                SequenceRange {
                    start: 2,
                    end_exclusive: Some(2),
                },
            ),
        ];

        for (input, expected) in cases {
            assert_eq!(parse_slice_range(input).unwrap(), expected, "range {input}");
        }

        let max_inclusive = format!("..={}", usize::MAX);
        assert_eq!(
            parse_slice_range(&max_inclusive).unwrap(),
            SequenceRange {
                start: 0,
                end_exclusive: None
            }
        );
    }

    #[test]
    fn rejects_malformed_and_reversed_ranges_with_input_errors() {
        for input in [
            "10", "one..2", "2..=nope", "10..2", "10..=1", "2..=1", "..=",
        ] {
            let error = parse_slice_range(input).unwrap_err();
            assert_eq!(
                error.kind(),
                std::io::ErrorKind::InvalidInput,
                "range {input}"
            );
        }
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
        let output = write_fasta_for_range(&[b">record\n", b"ACGT\n"], "..0", 10).unwrap();

        assert_eq!(output, b">record\n");
    }

    #[test]
    fn exclusive_and_inclusive_range_text_selects_expected_bases() {
        let input = &[b">record\n".as_slice(), b"ACGTNU\n".as_slice()];
        let exclusive = write_fasta_for_range(input, "2..4", 10).unwrap();
        let inclusive = write_fasta_for_range(input, "2..=4", 10).unwrap();

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
