use crate::{
    LocateArgs, ensure_distinct_input_output, is_nucleic_acid, output::with_output,
    read_lines_from_file,
};
use fasta_util::LinesInFile;
use std::{
    collections::VecDeque,
    fs::File,
    io::{self, BufRead, Write},
    path::Path,
};

pub(super) fn run(args: LocateArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(output) = &args.output
        && args.input != std::path::Path::new("-")
    {
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

    if args.input == Path::new("-") {
        let stdin = io::stdin();
        with_output(args.output.as_deref(), |writer| {
            locate_matches_reader(
                stdin.lock(),
                &pattern,
                &reverse_pattern,
                args.max_mismatch,
                writer,
            )
        })?;
    } else {
        let file = File::open(&args.input)?;
        // SAFETY: The input file must not be modified while the memory map is alive.
        let lines = unsafe { read_lines_from_file(file)? };
        with_output(args.output.as_deref(), |writer| {
            locate_matches(
                &lines,
                &pattern,
                &reverse_pattern,
                args.max_mismatch,
                writer,
            )
        })?;
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
    let mut scanner = LocateScanner::new(pattern, reverse_pattern, max_mismatch, writer);
    lines.try_for_each_line(|line_number, raw_line| scanner.process_line(line_number, raw_line))?;
    scanner.finish()
}

fn locate_matches_reader(
    mut reader: impl BufRead,
    pattern: &[u8],
    reverse_pattern: &[u8],
    max_mismatch: usize,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut scanner = LocateScanner::new(pattern, reverse_pattern, max_mismatch, writer);
    let mut line = Vec::new();
    let mut line_number = 0usize;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_number = line_number
            .checked_add(1)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "line number overflow"))?;
        scanner.process_line(line_number, &line)?;
    }
    scanner.finish()
}

struct LocateScanner<'a, W: Write> {
    pattern: &'a [u8],
    reverse_pattern: &'a [u8],
    max_mismatch: usize,
    writer: &'a mut W,
    record_id: Option<String>,
    saw_record: bool,
    position: u64,
    window: VecDeque<u8>,
}

impl<'a, W: Write> LocateScanner<'a, W> {
    fn new(
        pattern: &'a [u8],
        reverse_pattern: &'a [u8],
        max_mismatch: usize,
        writer: &'a mut W,
    ) -> Self {
        Self {
            pattern,
            reverse_pattern,
            max_mismatch,
            writer,
            record_id: None,
            saw_record: false,
            position: 0,
            window: VecDeque::with_capacity(pattern.len()),
        }
    }

    fn process_line(&mut self, line_number: usize, raw_line: &[u8]) -> io::Result<()> {
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
            self.record_id = Some(String::from_utf8_lossy(id).into_owned());
            self.saw_record = true;
            self.position = 0;
            self.window.clear();
            return Ok(());
        }
        if self.record_id.is_none() {
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
            self.position = self
                .position
                .checked_add(1)
                .ok_or_else(|| line_error(line_number, "sequence position overflow"))?;
            self.window.push_back(mask);
            if self.window.len() > self.pattern.len() {
                self.window.pop_front();
            }
            if self.window.len() != self.pattern.len() {
                continue;
            }

            let start = self.position - self.pattern.len() as u64 + 1;
            let id = self.record_id.as_deref().expect("record ID was validated");
            let (matches_forward, matches_reverse) = window_matches(
                &self.window,
                self.pattern,
                self.reverse_pattern,
                self.max_mismatch,
            );
            if matches_forward {
                writeln!(self.writer, "{id}\t{start}\t{}\t+", self.position)?;
            }
            if matches_reverse {
                writeln!(self.writer, "{id}\t{start}\t{}\t-", self.position)?;
            }
        }
        Ok(())
    }

    fn finish(self) -> io::Result<()> {
        if !self.saw_record {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no FASTA records found",
            ));
        }
        Ok(())
    }
}

fn window_matches(
    window: &VecDeque<u8>,
    forward_pattern: &[u8],
    reverse_pattern: &[u8],
    max_mismatch: usize,
) -> (bool, bool) {
    let mut forward_mismatches = 0;
    let mut reverse_mismatches = 0;
    for ((sequence_mask, forward_mask), reverse_mask) in window
        .iter()
        .copied()
        .zip(forward_pattern.iter().copied())
        .zip(reverse_pattern.iter().copied())
    {
        forward_mismatches += usize::from(sequence_mask & forward_mask == 0);
        reverse_mismatches += usize::from(sequence_mask & reverse_mask == 0);
        if forward_mismatches > max_mismatch && reverse_mismatches > max_mismatch {
            break;
        }
    }
    (
        forward_mismatches <= max_mismatch,
        reverse_mismatches <= max_mismatch,
    )
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
