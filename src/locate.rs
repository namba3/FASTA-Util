use crate::{
    LocateArgs, ensure_distinct_input_output, is_nucleic_acid, output::TemporaryOutput,
    read_lines_from_file,
};
use fasta_util::LinesInFile;
use std::{
    collections::VecDeque,
    fs::File,
    io::{self, BufWriter, Write},
};

pub(super) fn run(args: LocateArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(output) = &args.output {
        ensure_distinct_input_output(&args.input, output)?;
    }
    let (pattern, reverse_pattern) = parse_pattern(args.pattern.as_bytes())?;
    if args.max_mismatch > pattern.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--max-mismatch cannot be greater than the motif length",
        )
        .into());
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
    locate_matches(
        &lines,
        &pattern,
        &reverse_pattern,
        args.max_mismatch,
        &mut writer,
    )?;
    writer.flush()?;
    drop(writer);
    if let Some(temporary_output) = &mut temporary_output {
        temporary_output.commit()?;
    }
    Ok(())
}

fn parse_pattern(pattern: &[u8]) -> io::Result<(Vec<u8>, Vec<u8>)> {
    if pattern.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "motif cannot be empty",
        ));
    }
    let masks = pattern
        .iter()
        .copied()
        .map(|byte| {
            iupac_mask(byte).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "invalid IUPAC nucleotide symbol '{}' in motif",
                        char::from(byte)
                    ),
                )
            })
        })
        .collect::<io::Result<Vec<_>>>()?;
    let reverse = masks.iter().rev().copied().map(complement_mask).collect();
    Ok((masks, reverse))
}

fn locate_matches(
    lines: &LinesInFile,
    pattern: &[u8],
    reverse_pattern: &[u8],
    max_mismatch: usize,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut record_id = None::<String>;
    let mut saw_record = false;
    let mut position = 0u64;
    let mut window = VecDeque::with_capacity(pattern.len());

    lines.try_for_each_line(|line_number, raw_line| {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            let id = line[1..]
                .trim_ascii_start()
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            if id.is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            record_id = Some(String::from_utf8_lossy(id).into_owned());
            saw_record = true;
            position = 0;
            window.clear();
            return Ok(());
        }
        if record_id.is_none() {
            if line.is_empty() {
                return Ok(());
            }
            return Err(line_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        }

        for (column, byte) in line.iter().copied().enumerate() {
            if byte.is_ascii_whitespace() {
                return Err(line_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", column + 1),
                ));
            }
            if !is_nucleic_acid(byte) {
                return Err(line_error(
                    line_number,
                    &format!("invalid nucleotide '{}'", char::from(byte)),
                ));
            }
            let mask = iupac_mask(byte).expect("validated nucleotide must have an IUPAC mask");
            position = position
                .checked_add(1)
                .ok_or_else(|| line_error(line_number, "sequence position overflow"))?;
            window.push_back(mask);
            if window.len() > pattern.len() {
                window.pop_front();
            }
            if window.len() != pattern.len() {
                continue;
            }

            let start = position - pattern.len() as u64 + 1;
            let id = record_id.as_deref().expect("record ID was validated");
            if window_mismatches(&window, pattern) <= max_mismatch {
                writeln!(writer, "{id}\t{start}\t{position}\t+")?;
            }
            if window_mismatches(&window, reverse_pattern) <= max_mismatch {
                writeln!(writer, "{id}\t{start}\t{position}\t-")?;
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
    Ok(())
}

fn window_mismatches(window: &VecDeque<u8>, pattern: &[u8]) -> usize {
    window
        .iter()
        .copied()
        .zip(pattern.iter().copied())
        .filter(|(sequence_mask, motif_mask)| sequence_mask & motif_mask == 0)
        .count()
}

fn iupac_mask(byte: u8) -> Option<u8> {
    Some(match byte.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        b'-' => 0b1_0000,
        _ => return None,
    })
}

fn complement_mask(mask: u8) -> u8 {
    ((mask & 0b0001) << 3)
        | ((mask & 0b0010) << 1)
        | ((mask & 0b0100) >> 1)
        | ((mask & 0b1000) >> 3)
        | (mask & 0b1_0000)
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

#[cfg(test)]
mod tests {
    use super::{complement_mask, iupac_mask};

    #[test]
    fn maps_standard_iupac_symbols_to_base_sets() {
        let symbols = b"ACGTURYSWKMBDHVN-";
        let masks = symbols
            .iter()
            .map(|symbol| iupac_mask(*symbol).unwrap())
            .collect::<Vec<_>>();

        assert_eq!(
            masks,
            [
                0b0001, 0b0010, 0b0100, 0b1000, 0b1000, 0b0101, 0b1010, 0b0110, 0b1001, 0b1100,
                0b0011, 0b1110, 0b1101, 0b1011, 0b0111, 0b1111, 0b1_0000,
            ]
        );
    }

    #[test]
    fn complementing_an_iupac_mask_twice_returns_the_original_set() {
        for symbol in b"ACGTURYSWKMBDHVN-" {
            let mask = iupac_mask(*symbol).unwrap();
            assert_eq!(complement_mask(complement_mask(mask)), mask);
        }
    }
}
