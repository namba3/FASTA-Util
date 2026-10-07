use crate::{
    FormatArgs, ensure_distinct_input_output, output::TemporaryOutput, read_lines_from_file,
};
use fasta_util::LinesInFile;
use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

pub(super) fn run(args: FormatArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(output) = &args.output {
        ensure_distinct_input_output(&args.input, output)?;
    }

    let file = File::open(&args.input)?;
    // SAFETY: The input file must not change while its memory map is alive.
    let lines = unsafe { read_lines_from_file(file)? };
    let mut temporary_output = args
        .output
        .as_deref()
        .map(TemporaryOutput::create)
        .transpose()?;
    let output: Box<dyn Write> = match temporary_output.as_mut() {
        Some(temporary_output) => Box::new(temporary_output.take_file()?),
        None => Box::new(io::stdout().lock()),
    };
    let mut writer = BufWriter::new(output);
    format_records(&lines, &args, &mut writer)?;
    writer.flush()?;
    drop(writer);
    if let Some(temporary_output) = &mut temporary_output {
        temporary_output.commit()?;
    }
    Ok(())
}

fn format_records(
    lines: &LinesInFile,
    args: &FormatArgs,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut saw_record = false;
    let mut record_has_sequence = false;
    let mut sequence_line = Vec::new();

    lines.try_for_each_line(|line_number, raw_line| {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            if saw_record {
                finish_record(args.width, record_has_sequence, &mut sequence_line, writer)?;
            }
            let header = &line[1..];
            let normalized_header = if args.trim_header {
                header.trim_ascii()
            } else {
                header
            };
            if normalized_header.trim_ascii().is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            writer.write_all(b">")?;
            writer.write_all(normalized_header)?;
            writer.write_all(b"\n")?;
            saw_record = true;
            record_has_sequence = false;
            return Ok(());
        }

        if !saw_record {
            if line.is_empty() {
                return Ok(());
            }
            return Err(line_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        }

        for (column, mut byte) in line.iter().copied().enumerate() {
            if byte.is_ascii_whitespace() {
                return Err(line_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", column + 1),
                ));
            }
            if args.remove_gaps && byte == b'-' {
                continue;
            }
            if args.uppercase {
                byte.make_ascii_uppercase();
            } else if args.lowercase {
                byte.make_ascii_lowercase();
            }
            record_has_sequence = true;
            if args.width == 0 {
                writer.write_all(&[byte])?;
            } else {
                sequence_line.push(byte);
                if sequence_line.len() == args.width {
                    writer.write_all(&sequence_line)?;
                    writer.write_all(b"\n")?;
                    sequence_line.clear();
                }
            }
        }
        Ok(())
    })?;

    if !saw_record {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no FASTA records found",
        ));
    }
    finish_record(args.width, record_has_sequence, &mut sequence_line, writer)
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

fn strip_line_ending(line: &[u8]) -> &[u8] {
    match line.strip_suffix(b"\n") {
        Some(line) => line.strip_suffix(b"\r").unwrap_or(line),
        None => line,
    }
}

fn line_error(line_number: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("line {line_number}: {message}"),
    )
}
