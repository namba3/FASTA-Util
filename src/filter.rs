use crate::{
    FilterArgs, StatsSequenceType, ensure_distinct_input_output, is_amino_acid, is_nucleic_acid,
    line_error,
    output::{InputSource, with_output},
    read_lines_from_file,
    selection_bitmap::{SelectionBitmap, write_selected_records},
    strip_line_ending,
};
use fasta_util::LinesInFile;
use std::{fs::File, io};

#[derive(Default)]
struct Counts {
    length: u64,
    gc: u64,
    n: u64,
}

pub(super) fn run(args: FilterArgs) -> Result<(), Box<dyn std::error::Error>> {
    if args.min_len.is_none()
        && args.max_len.is_none()
        && args.min_gc.is_none()
        && args.max_gc.is_none()
        && args.max_n.is_none()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "provide at least one filter option",
        )
        .into());
    }
    if args.sequence_type == StatsSequenceType::Protein
        && (args.min_gc.is_some() || args.max_gc.is_some() || args.max_n.is_some())
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "GC and N filters are only available for nucleotide sequences",
        )
        .into());
    }
    if let (Some(minimum), Some(maximum)) = (args.min_len, args.max_len)
        && minimum > maximum
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-len cannot be greater than --max-len",
        )
        .into());
    }
    if let (Some(minimum), Some(maximum)) = (args.min_gc, args.max_gc)
        && minimum > maximum
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-gc cannot be greater than --max-gc",
        )
        .into());
    }
    if let (Some(input), Some(output)) = (&args.input, &args.output)
        && input != std::path::Path::new("-")
    {
        ensure_distinct_input_output(input, output)?;
    }

    let input = InputSource::from_optional_path(args.input.as_deref())?;
    let file = File::open(input.path())?;
    // SAFETY: The input file must not change while its memory map is alive.
    let lines = unsafe { read_lines_from_file(file)? };
    let selected = select_records(&lines, &args)?;

    with_output(args.output.as_deref(), |writer| {
        write_selected_records(&lines, &selected, writer)
    })?;
    Ok(())
}

fn select_records(lines: &LinesInFile, args: &FilterArgs) -> io::Result<SelectionBitmap> {
    let mut selected = SelectionBitmap::default();
    let mut current: Option<Counts> = None;
    let mut saw_record = false;
    let mut saw_protein_only = false;

    lines.try_for_each_line(|line_number, raw_line| {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            if let Some(counts) = current.take() {
                selected.push(matches_filter(&counts, args));
            }
            if line[1..].trim_ascii().is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            current = Some(Counts::default());
            saw_record = true;
            return Ok(());
        }

        let Some(counts) = current.as_mut() else {
            if line.is_empty() {
                return Ok(());
            }
            return Err(line_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        };

        for (index, byte) in line.iter().copied().enumerate() {
            if byte.is_ascii_whitespace() {
                return Err(line_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", index + 1),
                ));
            }
            let valid = match args.sequence_type {
                StatsSequenceType::Nucleotide => is_nucleic_acid(byte),
                StatsSequenceType::Protein | StatsSequenceType::Auto => is_amino_acid(byte),
            };
            if args.sequence_type == StatsSequenceType::Auto && !is_nucleic_acid(byte) {
                saw_protein_only = true;
            }
            if !valid {
                let kind = match args.sequence_type {
                    StatsSequenceType::Nucleotide | StatsSequenceType::Auto => "nucleotide",
                    StatsSequenceType::Protein => "protein",
                };
                return Err(line_error(
                    line_number,
                    &format!("invalid {kind} symbol '{}'", char::from(byte)),
                ));
            }
            counts.length = counts
                .length
                .checked_add(1)
                .ok_or_else(|| line_error(line_number, "sequence length overflow"))?;
            if matches!(byte, b'G' | b'g' | b'C' | b'c') {
                counts.gc = counts
                    .gc
                    .checked_add(1)
                    .ok_or_else(|| line_error(line_number, "GC count overflow"))?;
            }
            if matches!(byte, b'N' | b'n') {
                counts.n = counts
                    .n
                    .checked_add(1)
                    .ok_or_else(|| line_error(line_number, "N count overflow"))?;
            }
        }
        Ok(())
    })?;

    if let Some(counts) = current {
        selected.push(matches_filter(&counts, args));
    }
    if !saw_record {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no FASTA records found",
        ));
    }
    let is_protein = args.sequence_type == StatsSequenceType::Protein
        || (args.sequence_type == StatsSequenceType::Auto && saw_protein_only);
    if is_protein && (args.min_gc.is_some() || args.max_gc.is_some() || args.max_n.is_some()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "GC and N filters are only available for nucleotide sequences",
        ));
    }
    Ok(selected)
}

fn matches_filter(counts: &Counts, args: &FilterArgs) -> bool {
    if args.min_len.is_some_and(|minimum| counts.length < minimum)
        || args.max_len.is_some_and(|maximum| counts.length > maximum)
    {
        return false;
    }
    let gc_fraction = fraction(counts.gc, counts.length);
    let n_fraction = fraction(counts.n, counts.length);
    !(args.min_gc.is_some_and(|minimum| gc_fraction < minimum)
        || args.max_gc.is_some_and(|maximum| gc_fraction > maximum)
        || args.max_n.is_some_and(|maximum| n_fraction > maximum))
}

fn fraction(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64
    }
}
