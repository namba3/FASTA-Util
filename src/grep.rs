use crate::{
    GrepArgs, ensure_distinct_input_output, line_error,
    output::{InputSource, with_output},
    read_lines_from_file,
    selection_bitmap::{SelectionBitmap, write_selected_records},
    strip_line_ending,
};
use fasta_util::LinesInFile;
use std::{fs::File, io};

pub(super) fn run(args: GrepArgs) -> Result<(), Box<dyn std::error::Error>> {
    if args.pattern.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "pattern cannot be empty").into());
    }
    if let Some(output) = &args.output
        && args.input != std::path::Path::new("-")
    {
        ensure_distinct_input_output(&args.input, output)?;
    }

    let input = InputSource::from_optional_path(Some(&args.input))?;
    let file = File::open(input.path())?;
    // SAFETY: The input file must not change while its memory map is alive.
    let lines = unsafe { read_lines_from_file(file)? };
    let selected = select_records(&lines, &args)?;

    with_output(args.output.as_deref(), |writer| {
        write_selected_records(&lines, &selected, writer)
    })?;
    Ok(())
}

fn select_records(lines: &LinesInFile, args: &GrepArgs) -> io::Result<SelectionBitmap> {
    let pattern = args.pattern.as_bytes();
    let mut selected = SelectionBitmap::default();
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
