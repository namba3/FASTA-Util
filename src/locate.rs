use crate::{
    LocateArgs, ensure_distinct_input_output, for_each_reader_line, line_error,
    output::with_output, read_lines_from_file, strip_line_ending,
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
    reader: impl BufRead,
    pattern: &[u8],
    reverse_pattern: &[u8],
    max_mismatch: usize,
    writer: &mut impl Write,
) -> io::Result<()> {
    let mut scanner = LocateScanner::new(pattern, reverse_pattern, max_mismatch, writer);
    for_each_reader_line(reader, |line_number, line| {
        scanner.process_line(line_number, line)
    })?;
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
    bit_parallel: Option<BitParallelMatcher>,
}

impl<'a, W: Write> LocateScanner<'a, W> {
    fn new(
        pattern: &'a [u8],
        reverse_pattern: &'a [u8],
        max_mismatch: usize,
        writer: &'a mut W,
    ) -> Self {
        let bit_parallel = (pattern.len() <= u64::BITS as usize && max_mismatch <= 1)
            .then(|| BitParallelMatcher::new(pattern, reverse_pattern));
        Self {
            pattern,
            reverse_pattern,
            max_mismatch,
            writer,
            record_id: None,
            saw_record: false,
            position: 0,
            window: VecDeque::with_capacity(if bit_parallel.is_some() {
                0
            } else {
                pattern.len()
            }),
            bit_parallel,
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
            if let Some(matcher) = &mut self.bit_parallel {
                matcher.reset();
            }
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
            let Some(mask) = iupac_mask(byte) else {
                return Err(line_error(
                    line_number,
                    &format!("invalid nucleotide '{}'", char::from(byte)),
                ));
            };
            self.position = self
                .position
                .checked_add(1)
                .ok_or_else(|| line_error(line_number, "sequence position overflow"))?;
            let (matches_forward, matches_reverse) = if let Some(matcher) = &mut self.bit_parallel {
                let matches = matcher.advance(mask, self.max_mismatch);
                if self.position < self.pattern.len() as u64 {
                    continue;
                }
                matches
            } else {
                self.window.push_back(mask);
                if self.window.len() > self.pattern.len() {
                    self.window.pop_front();
                }
                if self.window.len() != self.pattern.len() {
                    continue;
                }
                window_matches(
                    &self.window,
                    self.pattern,
                    self.reverse_pattern,
                    self.max_mismatch,
                )
            };
            let start = self.position - self.pattern.len() as u64 + 1;
            let id = self.record_id.as_deref().expect("record ID was validated");
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

struct BitParallelMatcher {
    forward: ShiftAndState,
    reverse: Option<ShiftAndState>,
}

impl BitParallelMatcher {
    fn new(forward_pattern: &[u8], reverse_pattern: &[u8]) -> Self {
        Self {
            forward: ShiftAndState::new(forward_pattern),
            reverse: (forward_pattern != reverse_pattern)
                .then(|| ShiftAndState::new(reverse_pattern)),
        }
    }

    fn reset(&mut self) {
        self.forward.reset();
        if let Some(reverse) = &mut self.reverse {
            reverse.reset();
        }
    }

    fn advance(&mut self, sequence_mask: u8, max_mismatch: usize) -> (bool, bool) {
        let forward = self.forward.advance(sequence_mask, max_mismatch);
        let reverse = self.reverse.as_mut().map_or(forward, |matcher| {
            matcher.advance(sequence_mask, max_mismatch)
        });
        (forward, reverse)
    }
}

struct ShiftAndState {
    matching_positions: [u64; 17],
    exact_state: u64,
    one_mismatch_state: u64,
    final_position: u64,
}

impl ShiftAndState {
    fn new(pattern: &[u8]) -> Self {
        let mut matching_positions = [0u64; 17];
        for (index, pattern_mask) in pattern.iter().copied().enumerate() {
            let position = 1u64 << index;
            for (sequence_mask, matches) in matching_positions.iter_mut().enumerate().skip(1) {
                if pattern_mask & sequence_mask as u8 != 0 {
                    *matches |= position;
                }
            }
        }
        Self {
            matching_positions,
            exact_state: 0,
            one_mismatch_state: 0,
            final_position: 1u64 << (pattern.len() - 1),
        }
    }

    fn reset(&mut self) {
        self.exact_state = 0;
        self.one_mismatch_state = 0;
    }

    fn advance(&mut self, sequence_mask: u8, max_mismatch: usize) -> bool {
        let matching_positions = self.matching_positions[sequence_mask as usize];
        let previous_exact = self.exact_state;
        self.exact_state = ((previous_exact << 1) | 1) & matching_positions;
        if max_mismatch == 0 {
            return self.exact_state & self.final_position != 0;
        }

        let previous_one_mismatch = self.one_mismatch_state;
        self.one_mismatch_state =
            (((previous_one_mismatch << 1) | 1) & matching_positions) | ((previous_exact << 1) | 1);
        self.one_mismatch_state & self.final_position != 0
    }
}

fn window_matches(
    window: &VecDeque<u8>,
    forward_pattern: &[u8],
    reverse_pattern: &[u8],
    max_mismatch: usize,
) -> (bool, bool) {
    if forward_pattern == reverse_pattern {
        let mut mismatches = 0;
        for (sequence_mask, pattern_mask) in window.iter().copied().zip(forward_pattern.iter()) {
            mismatches += usize::from(sequence_mask & pattern_mask == 0);
            if mismatches > max_mismatch {
                return (false, false);
            }
        }
        let matched = mismatches <= max_mismatch;
        return (matched, matched);
    }

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

#[cfg(test)]
mod tests {
    use super::{BitParallelMatcher, complement_mask, iupac_mask, window_matches};
    use crate::is_nucleic_acid;
    use std::collections::VecDeque;

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
    fn iupac_mask_coverage_matches_the_nucleotide_validator() {
        for byte in 0..=u8::MAX {
            assert_eq!(
                iupac_mask(byte).is_some(),
                is_nucleic_acid(byte),
                "byte {byte}"
            );
        }
    }

    #[test]
    fn complementing_an_iupac_mask_twice_returns_the_original_set() {
        for symbol in b"ACGTURYSWKMBDHVN-" {
            let mask = iupac_mask(*symbol).unwrap();
            assert_eq!(complement_mask(complement_mask(mask)), mask);
        }
    }

    #[test]
    fn palindromic_patterns_match_both_strands_with_one_mismatch_count() {
        let pattern = [0b0001, 0b0010, 0b0010, 0b0001];
        let exact_window = VecDeque::from(pattern);
        let one_mismatch_window = VecDeque::from([0b0001, 0b0010, 0b0100, 0b0001]);
        let two_mismatch_window = VecDeque::from([0b0100, 0b0010, 0b0100, 0b0001]);

        assert_eq!(
            window_matches(&exact_window, &pattern, &pattern, 0),
            (true, true)
        );
        assert_eq!(
            window_matches(&one_mismatch_window, &pattern, &pattern, 1),
            (true, true)
        );
        assert_eq!(
            window_matches(&two_mismatch_window, &pattern, &pattern, 1),
            (false, false)
        );
    }

    #[test]
    fn bit_parallel_matching_agrees_with_window_matching_for_iupac_masks() {
        let patterns = [
            vec![0b0001],
            vec![0b0001, 0b1000],
            vec![0b0101, 0b1010],
            vec![0b0001, 0b0010, 0b0100],
            vec![0b0101, 0b1111, 0b1010, 0b1_0000],
            vec![0b0001; 64],
        ];
        let sequence = (0..257)
            .map(|index| [1, 2, 4, 8, 5, 10, 15, 16, 3, 12][index % 10])
            .collect::<Vec<_>>();

        for pattern in patterns {
            let reverse_pattern = pattern
                .iter()
                .rev()
                .copied()
                .map(complement_mask)
                .collect::<Vec<_>>();
            for max_mismatch in 0..=1 {
                let mut matcher = BitParallelMatcher::new(&pattern, &reverse_pattern);
                let mut window = VecDeque::with_capacity(pattern.len());

                for (position, mask) in sequence.iter().copied().enumerate() {
                    let actual = matcher.advance(mask, max_mismatch);
                    window.push_back(mask);
                    if window.len() > pattern.len() {
                        window.pop_front();
                    }
                    let expected = if window.len() == pattern.len() {
                        window_matches(&window, &pattern, &reverse_pattern, max_mismatch)
                    } else {
                        (false, false)
                    };
                    assert_eq!(
                        actual,
                        expected,
                        "pattern length {}, position {position}, mismatch {max_mismatch}",
                        pattern.len()
                    );
                }
            }
        }
    }
}
