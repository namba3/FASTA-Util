use clap::{Parser, Subcommand};
use crossbeam::channel::{Receiver, bounded};
use fasta_util::{is_nucleic_acid, read_lines_from_file, read_lines_from_stdin};
use std::io::{self, BufWriter, Write};

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
    input: Option<String>,
}

#[derive(Parser)]
struct SliceArgs {
    #[arg(
        short,
        long,
        help = "Specify input file\nIf omitted, read from standard input"
    )]
    input: Option<String>,

    #[arg(
        short,
        long,
        help = "Specify output file\nIf omitted, write to standard output"
    )]
    output: Option<String>,

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
            let lines = read_lines_from_file(input)?;
            count_sequence_bases(lines.map(Ok::<_, io::Error>))?
        }
        None => count_sequence_bases(read_lines_from_stdin())?,
    };

    println!("{len}");

    Ok(())
}

fn count_sequence_bases<T, I>(iter: I) -> io::Result<u64>
where
    T: AsRef<[u8]>,
    I: IntoIterator<Item = Result<T, io::Error>>,
{
    let mut count = 0u64;
    for line in iter {
        let line = line?;
        let line = line.as_ref();

        if line.first() == Some(&b'>') {
            continue;
        }

        let sequence = line.trim_ascii_start().trim_ascii_end();
        if sequence.is_empty() {
            continue;
        }
        if let Some(x) = sequence.iter().find(|x| !is_nucleic_acid(**x)) {
            return Err(invalid_nucleic_acid(*x));
        }

        count += sequence.len() as u64;
    }
    Ok(count)
}

fn invalid_nucleic_acid(byte: u8) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "invalid nucleic acid: {:?} (0x{byte:02x})",
            char::from(byte)
        ),
    )
}

#[derive(Debug, PartialEq, Eq)]
struct SequenceRange {
    start: usize,
    end_exclusive: Option<usize>,
}

fn invalid_range(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
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
    let range = parse_slice_range(&args.range)?;
    let writer_options = WriterOptions {
        chars_per_line: args.chars_per_line,
        start: range.start,
        end_exclusive: range.end_exclusive,
    };

    let output: Box<dyn std::io::Write> = if let Some(path) = args.output {
        Box::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .open(path)?,
        )
    } else {
        Box::new(std::io::stdout().lock())
    };

    let mut writer = Writer::new(output, writer_options);

    let (hndl, write_result) = match args.input {
        Some(input) => {
            let input = std::fs::OpenOptions::new().read(true).open(input)?;
            let (tx, rx) = bounded(LINE_CHANNEL_CAPACITY);
            let hndl = std::thread::spawn(move || -> Result<(), std::io::Error> {
                let mut lines = read_lines_from_file(input)?;
                while let Some(line) = lines.next() {
                    if tx.send(Ok(line)).is_err() {
                        break;
                    }
                }
                Ok(())
            });

            let write_result = writer.run(rx);
            (hndl, write_result)
        }
        None => {
            let (tx, rx) = bounded(LINE_CHANNEL_CAPACITY);
            let hndl = std::thread::spawn(move || -> Result<(), std::io::Error> {
                let mut lines = read_lines_from_stdin();
                while let Some(line) = lines.next() {
                    if tx.send(line).is_err() {
                        break;
                    }
                }

                Ok(())
            });

            let write_result = writer.run(rx);
            (hndl, write_result)
        }
    };

    let read_result = hndl.join().unwrap_or_else(|_| {
        Err(io::Error::new(
            io::ErrorKind::Other,
            "input reader thread panicked",
        ))
    });
    write_result?;
    read_result?;

    Ok(())
}

struct WriterOptions {
    chars_per_line: usize,
    start: usize,
    end_exclusive: Option<usize>,
}
struct Writer<T: std::io::Write> {
    inner: BufWriter<T>,
    options: WriterOptions,
}
impl<T: std::io::Write> Writer<T> {
    fn new(inner: T, options: WriterOptions) -> Self {
        Self {
            inner: BufWriter::new(inner),
            options,
        }
    }
    fn run<Buf: AsRef<[u8]>>(
        &mut self,
        rx: Receiver<Result<Buf, io::Error>>,
    ) -> Result<(), std::io::Error> {
        let WriterOptions {
            chars_per_line,
            start,
            end_exclusive,
        } = &mut self.options;
        let writer = &mut self.inner;

        let mut cnt = 0usize;
        let mut written = 0usize;

        while let Ok(line) = rx.recv() {
            let line = line?;
            let buf = line.as_ref();

            if let Some(b'>') = buf.first() {
                if written > 0 && written % *chars_per_line != 0 {
                    writer.write_all(b"\n")?;
                }
                writer.write_all(&*buf)?;
                if !buf.ends_with(b"\n") {
                    writer.write_all(b"\n")?;
                }
                continue;
            }

            let buf = buf.trim_ascii_start().trim_ascii_end();
            if buf.len() == 0 {
                continue;
            }
            if let Some(x) = buf.iter().find(|x| !is_nucleic_acid(**x)) {
                return Err(invalid_nucleic_acid(*x));
            }

            let line_end = cnt + buf.len();
            let start_in_line = (*start).saturating_sub(cnt);
            let end_in_line = end_exclusive
                .map(|end| end.saturating_sub(cnt).min(buf.len()))
                .unwrap_or(buf.len());

            if start_in_line >= end_in_line {
                cnt = line_end;
                if end_exclusive.is_some_and(|end| end <= cnt) {
                    break;
                }
                continue;
            }

            let mut bases = &buf[start_in_line..end_in_line];
            let line_written = if start_in_line == 0 {
                written % *chars_per_line
            } else {
                0
            };
            let mut line_remain = *chars_per_line - line_written;
            written += bases.len();

            while line_remain <= bases.len() {
                writer.write_all(&bases[..line_remain as usize])?;
                writer.write_all(b"\n")?;
                bases = &bases[line_remain as usize..];
                line_remain = *chars_per_line;
            }

            writer.write_all(bases)?;

            cnt = line_end;
            if end_exclusive.is_some_and(|end| end <= cnt) {
                if written > 0 && written % *chars_per_line != 0 {
                    writer.write_all(b"\n")?;
                }
                break;
            }
        }

        writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Args, SequenceRange, Writer, WriterOptions, count_sequence_bases, parse_slice_range,
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
    fn slice_rejects_invalid_sequence_symbols() {
        let error = write_fasta_result(&[b"ACX\n"], options(10, None, None)).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("invalid nucleic acid"));
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
