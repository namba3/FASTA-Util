#![feature(byte_slice_trim_ascii)]

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
    chars_per_line: u64,
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
    fn count<T: AsRef<[u8]>, I: Iterator<Item = T>>(mut iter: I) -> u64 {
        let mut cnt = 0u64;
        while let Some(buf) = iter.next() {
            let buf = buf.as_ref();

            if buf[0] == b'>' {
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

            cnt += buf.len() as u64;
        }
        cnt
    }

    let len = match args.input {
        Some(input) => {
            let input = std::fs::OpenOptions::new().read(true).open(input)?;
            let lines = read_lines_from_file(input)?;
            count(lines)
        }
        None => {
            let lines = read_lines_from_stdin().filter_map(|x| x.ok());
            count(lines)
        }
    };

    println!("{len}");

    Ok(())
}

fn slice(args: SliceArgs) -> Result<(), Box<dyn std::error::Error>> {
    let writer_options = {
        let (start, end_inclusive) = {
            let range = args.range;
            let (start, end) = range.split_once("..").expect("invalid range");
            let start = if start.len() == 0 {
                None
            } else {
                Some(start.parse::<u64>()?)
            };
            let end_inclusive = if end.len() == 0 {
                None
            } else {
                if end.starts_with("=") {
                    let (_, end) = end.split_once("=").unwrap();
                    Some(end.parse::<u64>()?)
                } else {
                    Some(end.parse::<u64>()? - 1)
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
    chars_per_line: u64,
    start: Option<u64>,
    end_inclusive: Option<u64>,
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

        let mut cnt = 0u64;
        let mut written = 0u64;

        while let Ok(buf) = rx.recv() {
            let buf = buf.as_ref();

            if buf[0] == b'>' {
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

            let s = if let Some(n) = start.take() {
                if (cnt + buf.len() as u64) < n {
                    cnt += buf.len() as u64;
                    *start = n.into();
                    continue;
                }
                *start = None;
                (n - cnt) as usize
            } else {
                0
            };
            let e = if let Some(n) = end_inclusive {
                if *n <= (cnt + buf.len() as u64) {
                    (*n - cnt) as usize
                } else {
                    buf.len() - 1
                }
            } else {
                buf.len() - 1
            };

            let mut bases = &buf[s..=e];
            let line_written = if s == 0 { written % *chars_per_line } else { 0 };
            let mut line_remain = *chars_per_line - line_written;
            written += bases.len() as u64;

            while line_remain <= (bases.len() as u64) {
                writer.write_all(&bases[..line_remain as usize])?;
                writer.write_all(b"\n")?;
                bases = &bases[line_remain as usize..];
                line_remain = *chars_per_line;
            }

            writer.write_all(bases)?;

            cnt += buf.len() as u64;
            if let Some(n) = end_inclusive {
                if *n <= cnt {
                    writer.write_all(b"\n")?;
                    break;
                }
            }
        }

        writer.flush()
    }
}
