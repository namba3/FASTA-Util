use crate::parallel::{first_header_offset, line_chunks, process_chunks, worker_count};
use crate::{
    StatsSequenceType, for_each_reader_line, is_amino_acid, is_nucleic_acid, line_error,
    read_lines_from_file, strip_line_ending,
};
use fasta_util::LinesInFile;
use std::{
    fs::File,
    io::{self, Write},
    num::NonZeroUsize,
    path::Path,
};

const PROTEIN_SYMBOLS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const PROTEIN_EXTENDED_SYMBOLS: &[u8] = b"BJOUXZ*-";
const NUCLEOTIDE_AMBIGUOUS_SYMBOLS: &[u8] = b"KSYMWRBDHV-";

struct Composition {
    counts: [u64; 256],
    total: u64,
    saw_record: bool,
    saw_t: bool,
    saw_u: bool,
    saw_protein_only: bool,
}

impl Default for Composition {
    fn default() -> Self {
        Self {
            counts: [0; 256],
            total: 0,
            saw_record: false,
            saw_t: false,
            saw_u: false,
            saw_protein_only: false,
        }
    }
}

pub(super) fn run(
    path: Option<&Path>,
    sequence_type: StatsSequenceType,
    threads: Option<NonZeroUsize>,
) -> io::Result<()> {
    let composition = match path {
        None => {
            let stdin = io::stdin();
            collect_from_reader(stdin.lock(), sequence_type)?
        }
        Some(path) if path == Path::new("-") => {
            let stdin = io::stdin();
            collect_from_reader(stdin.lock(), sequence_type)?
        }
        Some(path) => {
            let file = File::open(path)?;
            // SAFETY: This command only reads the input; the file must not be modified
            // while the memory map is alive.
            let lines = unsafe { read_lines_from_file(file)? };
            collect_from_lines(&lines, sequence_type, threads)?
        }
    };
    let stdout = io::stdout();
    write_composition(&composition, sequence_type, &mut stdout.lock())
}

fn collect_from_lines(
    lines: &LinesInFile,
    sequence_type: StatsSequenceType,
    requested_threads: Option<NonZeroUsize>,
) -> io::Result<Composition> {
    let bytes = lines.as_bytes();
    let chunks = line_chunks(bytes, worker_count(bytes.len(), requested_threads));
    let first_header = first_header_offset(bytes);
    let partials = process_chunks(&chunks, |chunk| {
        let mut composition = Composition {
            saw_record: first_header.is_some_and(|offset| chunk.range.start > offset),
            ..Composition::default()
        };
        for (line_index, raw_line) in chunk
            .bytes
            .split_inclusive(|byte| *byte == b'\n')
            .enumerate()
        {
            composition.process_line(
                chunk.start_line + line_index,
                strip_line_ending(raw_line),
                sequence_type,
            )?;
        }
        Ok(composition)
    });

    let mut composition = Composition::default();
    for partial in partials {
        let partial = partial?;
        composition.total = composition
            .total
            .checked_add(partial.total)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "sequence length overflow")
            })?;
        for (total, count) in composition.counts.iter_mut().zip(partial.counts) {
            *total = total.checked_add(count).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "symbol count overflow")
            })?;
        }
        composition.saw_record |= partial.saw_record;
        composition.saw_t |= partial.saw_t;
        composition.saw_u |= partial.saw_u;
        composition.saw_protein_only |= partial.saw_protein_only;
    }
    composition.finish()
}

fn collect_from_reader(
    reader: impl io::BufRead,
    sequence_type: StatsSequenceType,
) -> io::Result<Composition> {
    let mut composition = Composition::default();
    for_each_reader_line(reader, |line_number, line| {
        composition.process_line(line_number, strip_line_ending(line), sequence_type)
    })?;
    composition.finish()
}

impl Composition {
    fn process_line(
        &mut self,
        line_number: usize,
        line: &[u8],
        sequence_type: StatsSequenceType,
    ) -> io::Result<()> {
        if line.first() == Some(&b'>') {
            if line[1..].trim_ascii().is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            self.saw_record = true;
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

        for (index, byte) in line.iter().copied().enumerate() {
            if byte.is_ascii_whitespace() {
                return Err(line_error(
                    line_number,
                    &format!("whitespace in sequence at column {}", index + 1),
                ));
            }
            let valid = match sequence_type {
                StatsSequenceType::Nucleotide => is_nucleic_acid(byte),
                StatsSequenceType::Protein | StatsSequenceType::Auto => is_amino_acid(byte),
            };
            if !valid {
                let kind = if sequence_type == StatsSequenceType::Protein {
                    "protein"
                } else {
                    "nucleotide"
                };
                return Err(line_error(
                    line_number,
                    &format!("invalid {kind} symbol '{}'", char::from(byte)),
                ));
            }

            let symbol = byte.to_ascii_uppercase();
            self.counts[symbol as usize] = self.counts[symbol as usize]
                .checked_add(1)
                .ok_or_else(|| line_error(line_number, "symbol count overflow"))?;
            self.total = self
                .total
                .checked_add(1)
                .ok_or_else(|| line_error(line_number, "sequence length overflow"))?;
            self.saw_t |= symbol == b'T';
            self.saw_u |= symbol == b'U';
            self.saw_protein_only |=
                sequence_type == StatsSequenceType::Auto && !is_nucleic_acid(byte);
        }
        Ok(())
    }

    fn finish(self) -> io::Result<Self> {
        if !self.saw_record {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no FASTA records found",
            ));
        }
        Ok(self)
    }
}

fn write_composition(
    composition: &Composition,
    sequence_type: StatsSequenceType,
    writer: &mut impl Write,
) -> io::Result<()> {
    let protein = match sequence_type {
        StatsSequenceType::Protein => true,
        StatsSequenceType::Nucleotide => false,
        StatsSequenceType::Auto => composition.saw_protein_only,
    };

    if protein {
        for symbol in PROTEIN_SYMBOLS {
            write_symbol(composition, *symbol, writer)?;
        }
        for symbol in PROTEIN_EXTENDED_SYMBOLS {
            if composition.counts[*symbol as usize] > 0 {
                write_symbol(composition, *symbol, writer)?;
            }
        }
        return Ok(());
    }

    for symbol in b"ACG" {
        write_symbol(composition, *symbol, writer)?;
    }
    if composition.saw_t || !composition.saw_u {
        write_symbol(composition, b'T', writer)?;
    }
    if composition.saw_u {
        write_symbol(composition, b'U', writer)?;
    }
    write_symbol(composition, b'N', writer)?;
    for symbol in NUCLEOTIDE_AMBIGUOUS_SYMBOLS {
        if composition.counts[*symbol as usize] > 0 {
            write_symbol(composition, *symbol, writer)?;
        }
    }
    let gc = composition.counts[b'G' as usize] + composition.counts[b'C' as usize];
    writeln!(writer, "GC\t{:.2}%", percentage(gc, composition.total))
}

fn write_symbol(composition: &Composition, symbol: u8, writer: &mut impl Write) -> io::Result<()> {
    writeln!(
        writer,
        "{}\t{:.2}%",
        char::from(symbol),
        percentage(composition.counts[symbol as usize], composition.total)
    )
}

fn percentage(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 * 100.0 / total as f64
    }
}
