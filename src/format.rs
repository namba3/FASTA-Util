use crate::{
    FormatArgs, ensure_distinct_input_output, for_each_reader_line, line_error,
    output::with_output, read_lines_from_file, strip_line_ending,
};
use fasta_util::LinesInFile;
use std::{
    fs::File,
    io::{self, BufRead, Write},
    path::Path,
};

pub(super) fn run(args: FormatArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let (Some(input), Some(output)) = (&args.input, &args.output)
        && input != std::path::Path::new("-")
    {
        ensure_distinct_input_output(input, output)?;
    }

    match args.input.as_deref() {
        None => {
            let stdin = io::stdin();
            with_output(args.output.as_deref(), |writer| {
                format_reader(stdin.lock(), &args, writer)
            })?;
        }
        Some(path) if path == Path::new("-") => {
            let stdin = io::stdin();
            with_output(args.output.as_deref(), |writer| {
                format_reader(stdin.lock(), &args, writer)
            })?;
        }
        Some(path) => {
            let file = File::open(path)?;
            // SAFETY: The input file must not change while its memory map is alive.
            let lines = unsafe { read_lines_from_file(file)? };
            with_output(args.output.as_deref(), |writer| {
                format_records(&lines, &args, writer)
            })?;
        }
    }
    Ok(())
}

fn format_records(
    lines: &LinesInFile,
    args: &FormatArgs,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut formatter = Formatter::new(args, writer);
    lines.try_for_each_line(|line_number, line| formatter.process_line(line_number, line))?;
    formatter.finish()
}

fn format_reader(
    reader: impl BufRead,
    args: &FormatArgs,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut formatter = Formatter::new(args, writer);
    for_each_reader_line(reader, |line_number, line| {
        formatter.process_line(line_number, line)
    })?;
    formatter.finish()
}

struct Formatter<'a, W: Write> {
    args: &'a FormatArgs,
    writer: &'a mut W,
    saw_record: bool,
    record_has_sequence: bool,
    sequence_line: Vec<u8>,
}

impl<'a, W: Write> Formatter<'a, W> {
    fn new(args: &'a FormatArgs, writer: &'a mut W) -> Self {
        Self {
            args,
            writer,
            saw_record: false,
            record_has_sequence: false,
            sequence_line: Vec::new(),
        }
    }

    fn process_line(&mut self, line_number: usize, raw_line: &[u8]) -> io::Result<()> {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            if self.saw_record {
                finish_record(
                    self.args.width,
                    self.record_has_sequence,
                    &mut self.sequence_line,
                    self.writer,
                )?;
            }
            let header = &line[1..];
            let normalized_header = if self.args.trim_header {
                header.trim_ascii()
            } else {
                header
            };
            if normalized_header.trim_ascii().is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            self.writer.write_all(b">")?;
            self.writer.write_all(normalized_header)?;
            self.writer.write_all(b"\n")?;
            self.saw_record = true;
            self.record_has_sequence = false;
            return Ok(());
        }

        if !self.saw_record {
            if line.is_empty() {
                return Ok(());
            }
            return Err(line_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        }

        if !self.args.uppercase && !self.args.lowercase && !self.args.remove_gaps {
            if let Some(column) = line.iter().position(u8::is_ascii_whitespace) {
                return Err(line_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", column + 1),
                ));
            }
            if !line.is_empty() {
                self.record_has_sequence = true;
                if self.args.width == 0 {
                    self.writer.write_all(line)?;
                } else {
                    self.append_sequence(line)?;
                }
            }
            return Ok(());
        }

        for (column, mut byte) in line.iter().copied().enumerate() {
            if byte.is_ascii_whitespace() {
                return Err(line_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", column + 1),
                ));
            }
            if self.args.remove_gaps && byte == b'-' {
                continue;
            }
            if self.args.uppercase {
                byte.make_ascii_uppercase();
            } else if self.args.lowercase {
                byte.make_ascii_lowercase();
            }
            self.record_has_sequence = true;
            if self.args.width == 0 {
                self.writer.write_all(&[byte])?;
            } else {
                self.sequence_line.push(byte);
                if self.sequence_line.len() == self.args.width {
                    self.writer.write_all(&self.sequence_line)?;
                    self.writer.write_all(b"\n")?;
                    self.sequence_line.clear();
                }
            }
        }
        Ok(())
    }

    fn append_sequence(&mut self, mut sequence: &[u8]) -> io::Result<()> {
        while !sequence.is_empty() {
            let available = self.args.width - self.sequence_line.len();
            let count = available.min(sequence.len());
            self.sequence_line.extend_from_slice(&sequence[..count]);
            sequence = &sequence[count..];
            if self.sequence_line.len() == self.args.width {
                self.writer.write_all(&self.sequence_line)?;
                self.writer.write_all(b"\n")?;
                self.sequence_line.clear();
            }
        }
        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        if !self.saw_record {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no FASTA records found",
            ));
        }
        finish_record(
            self.args.width,
            self.record_has_sequence,
            &mut self.sequence_line,
            self.writer,
        )
    }
}

fn finish_record(
    width: usize,
    has_sequence: bool,
    sequence_line: &mut Vec<u8>,
    writer: &mut impl Write,
) -> io::Result<()> {
    if width == 0 {
        if has_sequence {
            writer.write_all(b"\n")?;
        }
    } else if !sequence_line.is_empty() {
        writer.write_all(sequence_line)?;
        writer.write_all(b"\n")?;
        sequence_line.clear();
    }
    Ok(())
}
