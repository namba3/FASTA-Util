mod matcher;

use crate::{
    LocateArgs, ensure_distinct_input_output, for_each_reader_line, line_error,
    output::with_output, read_lines_from_file, strip_line_ending,
};
use fasta_util::LinesInFile;
use matcher::{
    BitParallelMatcher, MAX_BIT_PARALLEL_MISMATCHES, iupac_mask, parse_pattern, window_matches,
};
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
        let bit_parallel = (max_mismatch <= MAX_BIT_PARALLEL_MISMATCHES)
            .then(|| BitParallelMatcher::new(pattern, reverse_pattern, max_mismatch));
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
                let matches = matcher.advance(mask);
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
