use crate::{
    GrepArgs, ensure_distinct_input_output, output::TemporaryOutput, read_lines_from_file,
};
use fasta_util::LinesInFile;
use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

pub(super) fn run(args: GrepArgs) -> Result<(), Box<dyn std::error::Error>> {
    if args.pattern.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "pattern cannot be empty").into());
    }
    if let Some(output) = &args.output {
        ensure_distinct_input_output(&args.input, output)?;
    }

    let file = File::open(&args.input)?;
    // SAFETY: The input file must not change while its memory map is alive.
    let lines = unsafe { read_lines_from_file(file)? };
    let selected = select_records(&lines, &args)?;

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
    write_selected_records(&lines, &selected, &mut writer)?;
    writer.flush()?;
    drop(writer);
    if let Some(temporary_output) = &mut temporary_output {
        temporary_output.commit()?;
    }
    Ok(())
}

fn select_records(lines: &LinesInFile, args: &GrepArgs) -> io::Result<Vec<bool>> {
    let pattern = args.pattern.as_bytes();
    let mut selected = Vec::new();
    let mut saw_record = false;

    lines.try_for_each_line(|line_number, raw_line| {
        if raw_line.first() == Some(&b'>') {
            if raw_line[1..].trim_ascii().is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            let matched = contains(strip_line_ending(raw_line), pattern, args.ignore_case);
            selected.push(matched != args.invert_match);
            saw_record = true;
        } else if !saw_record && !strip_line_ending(raw_line).is_empty() {
            return Err(line_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        }
        Ok(())
    })?;

    if !saw_record {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no FASTA records found",
        ));
    }
    Ok(selected)
}

fn contains(haystack: &[u8], needle: &[u8], ignore_case: bool) -> bool {
    if needle.len() > haystack.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|window| {
        window
            .iter()
            .copied()
            .zip(needle.iter().copied())
            .all(|(left, right)| {
                if ignore_case {
                    left.eq_ignore_ascii_case(&right)
                } else {
                    left == right
                }
            })
    })
}

fn write_selected_records(
    lines: &LinesInFile,
    selected: &[bool],
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut record_index = None;
    lines.try_for_each_line(|_, raw_line| {
        if raw_line.first() == Some(&b'>') {
            record_index = Some(record_index.map_or(0, |index: usize| index + 1));
        }
        if record_index.is_some_and(|index| selected[index]) {
            writer.write_all(raw_line)?;
        }
        Ok(())
    })
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
