use clap::{Parser, Subcommand};
use core::panic;
use crossbeam::channel::{unbounded, Receiver};
use fasta_util::{is_nucleic_acid, read_lines_from_file, read_lines_from_stdin};
use std::io::{BufWriter, Write};

#[derive(Parser)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    sub: SubCommand,
}

#[derive(Subcommand)]
enum SubCommand {
    #[clap(about = "Count the total length of the sequence")]
    Len(LenArgs),
    #[clap(about = "Cut out a part of the sequence")]
    Slice(SliceArgs),
}

#[derive(Parser)]
struct LenArgs {
    #[clap(
        short,
        long,
        help = "Specify input file\nIf omitted, read from standard input"
    )]
    input: Option<String>,
}

#[derive(Parser)]
struct SliceArgs {
    #[clap(
        short,
        long,
        help = "Specify input file\nIf omitted, read from standard input"
    )]
    input: Option<String>,

    #[clap(
        short,
        long,
        help = "Specify output file\nIf omitted, write to standard output"
    )]
    output: Option<String>,

    #[clap(
        long,
        default_value = "..",
        help = "Specify slice range\nexamples:\n\t2..10\tmeans [2,10)\n\t2..=10\tmeans [2,10]\n\t..10\tmeans [0,10)\n\t2..\tmeans [2,∞)\n\t..\tmeans [0,∞)\n"
    )]
    range: String,

    #[clap(
        long,
        default_value_t = 60,
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

fn len(args: LenArgs) -> Result<(), Box<dyn std::error::Error>> {
    let len = match args.input {
        Some(input) => {
            let input = std::fs::OpenOptions::new().read(true).open(input)?;
            let lines = read_lines_from_file(input)?;
            count_sequence_bases(lines)
        }
        None => {
            let lines = read_lines_from_stdin().filter_map(|x| x.ok());
            count_sequence_bases(lines)
        }
    };

    println!("{len}");

    Ok(())
}

fn count_sequence_bases<T: AsRef<[u8]>, I: IntoIterator<Item = T>>(iter: I) -> u64 {
    let mut count = 0u64;
    for line in iter {
        let line = line.as_ref();

        if line.first() == Some(&b'>') {
            continue;
        }

        let sequence = line.trim_ascii_start().trim_ascii_end();
        if sequence.is_empty() {
            continue;
        }
        if let Some(x) = sequence.iter().find(|x| !is_nucleic_acid(**x)) {
            panic!(
                "invalid nucleic acid: '{}' (0x{x:0x})",
                char::from_u32(*x as u32).unwrap()
            );
        }

        count += sequence.len() as u64;
    }
    count
}

fn slice(args: SliceArgs) -> Result<(), Box<dyn std::error::Error>> {
    let writer_options = {
        let (start, end_inclusive) = {
            let range = args.range;
            let (start, end) = range.split_once("..").expect("invalid range");
            let start = if start.len() == 0 {
                None
            } else {
                Some(start.parse::<usize>()?)
            };
            let end_inclusive = if end.len() == 0 {
                None
            } else {
                if end.starts_with("=") {
                    let (_, end) = end.split_once("=").unwrap();
                    Some(end.parse::<usize>()?)
                } else {
                    Some(end.parse::<usize>()? - 1)
                }
            };

            if let (Some(start), Some(end_inclusive)) = (start, end_inclusive) {
                if end_inclusive < start {
                    panic!("invalid range");
                }
            }

            (start, end_inclusive)
        };
        let chars_per_line = args.chars_per_line.max(1);

        WriterOptions {
            chars_per_line,
            start,
            end_inclusive,
        }
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

    let hndl = match args.input {
        Some(input) => {
            let input = std::fs::OpenOptions::new().read(true).open(input)?;
            let (tx, rx) = unbounded();
            let hndl = std::thread::spawn(move || -> Result<(), std::io::Error> {
                let mut lines = read_lines_from_file(input)?;
                while let Some(line) = lines.next() {
                    if let Err(_) = tx.send(line) {
                        break;
                    };
                }
                Ok(())
            });

            writer.run(rx)?;
            hndl
        }
        None => {
            let (tx, rx) = unbounded();
            let hndl = std::thread::spawn(move || -> Result<(), std::io::Error> {
                let mut lines = read_lines_from_stdin();
                while let Some(Ok(line)) = lines.next() {
                    if let Err(_) = tx.send(line) {
                        break;
                    };
                }

                Ok(())
            });

            writer.run(rx)?;
            hndl
        }
    };

    match hndl.join() {
        Ok(x) => x?,
        Err(_why) => {
            panic!("error: read thread panicked.");
        }
    }

    Ok(())
}

struct WriterOptions {
    chars_per_line: usize,
    start: Option<usize>,
    end_inclusive: Option<usize>,
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
    fn run<Buf: AsRef<[u8]>>(&mut self, rx: Receiver<Buf>) -> Result<(), std::io::Error> {
        let WriterOptions {
            chars_per_line,
            start,
            end_inclusive,
        } = &mut self.options;
        let writer = &mut self.inner;

        let mut cnt = 0usize;
        let mut written = 0usize;

        while let Ok(buf) = rx.recv() {
            let buf = buf.as_ref();

            if let Some(b'>') = buf.first() {
                if written > 0 && written % *chars_per_line != 0 {
                    writer.write_all(b"\n")?;
                }
                writer.write_all(&*buf)?;
                continue;
            }

            let buf = buf.trim_ascii_start().trim_ascii_end();
            if buf.len() == 0 {
                continue;
            }
            if let Some(x) = buf.iter().find(|x| !is_nucleic_acid(**x)) {
                panic!(
                    "invalid nucleic acid: '{}' (0x{x:0x})",
                    char::from_u32(*x as u32).unwrap()
                );
            }

            let s = match start.take() {
                Some(n) if cnt + buf.len() <= n => {
                    cnt += buf.len();
                    *start = n.into();
                    continue;
                }
                Some(n) => {
                    *start = None;
                    (n - cnt) as usize
                }
                None => 0,
            };
            let e = match end_inclusive {
                Some(n) if *n < cnt + buf.len() => *n - cnt,
                _ => buf.len() - 1,
            };

            let mut bases = &buf[s..=e];
            let line_written = if s == 0 { written % *chars_per_line } else { 0 };
            let mut line_remain = *chars_per_line - line_written;
            written += bases.len();

            while line_remain <= bases.len() {
                writer.write_all(&bases[..line_remain as usize])?;
                writer.write_all(b"\n")?;
                bases = &bases[line_remain as usize..];
                line_remain = *chars_per_line;
            }

            writer.write_all(bases)?;

            cnt += buf.len();
            if let Some(n) = end_inclusive {
                if *n < cnt {
                    writer.write_all(b"\n")?;
                    break;
                }
            }
        }

        writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::{count_sequence_bases, Writer, WriterOptions};
    use crossbeam::channel::unbounded;

    fn write_fasta(lines: &[&[u8]], options: WriterOptions) -> Vec<u8> {
        let (tx, rx) = unbounded();
        for line in lines {
            tx.send(*line).unwrap();
        }
        drop(tx);

        let mut writer = Writer::new(Vec::new(), options);
        writer.run(rx).unwrap();
        writer.inner.into_inner().unwrap()
    }

    fn options(
        chars_per_line: usize,
        start: Option<usize>,
        end_inclusive: Option<usize>,
    ) -> WriterOptions {
        WriterOptions {
            chars_per_line,
            start,
            end_inclusive,
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

        assert_eq!(count_sequence_bases(lines), 7);
    }

    #[test]
    #[should_panic(expected = "invalid nucleic acid")]
    fn len_rejects_invalid_sequence_symbols() {
        count_sequence_bases([b"ACX".as_slice()]);
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
            options(2, Some(2), Some(4)),
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
            options(10, Some(4), Some(5)),
        );

        assert_eq!(output, b">record\nNU\n");
    }

    #[test]
    fn slice_stops_at_the_inclusive_end() {
        let output = write_fasta(&[b">record\n", b"ACGT\n"], options(10, None, Some(2)));

        assert_eq!(output, b">record\nACG\n");
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
    #[should_panic(expected = "invalid nucleic acid")]
    fn slice_rejects_invalid_sequence_symbols() {
        let _ = write_fasta(&[b"ACX\n"], options(10, None, None));
    }
}
