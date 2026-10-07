use crate::{
    RevcompArgs, ensure_distinct_input_output, is_nucleic_acid, output::TemporaryOutput,
    read_lines_from_file,
};
use fasta_util::LinesInFile;
use std::{
    fs::File,
    io::{self, BufWriter, Read, Seek, SeekFrom, Write},
};

const REVERSE_READ_SIZE: usize = 64 * 1024;

struct Record {
    header: Vec<u8>,
    sequence_start: u64,
    sequence_end: u64,
    saw_t: bool,
    saw_u: bool,
}

pub(super) fn run(args: RevcompArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(output) = &args.output {
        ensure_distinct_input_output(&args.input, output)?;
    }

    let file = File::open(&args.input)?;
    // SAFETY: The input file must not change while its memory map is alive.
    let lines = unsafe { read_lines_from_file(file)? };
    let records = scan_records(&lines)?;

    let mut input = File::open(&args.input)?;
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
    let mut block = vec![0; REVERSE_READ_SIZE];
    let mut sequence_line = Vec::new();

    for record in &records {
        reverse_complement_record(
            &mut input,
            record,
            args.chars_per_line,
            &mut block,
            &mut sequence_line,
            &mut writer,
        )?;
    }
    writer.flush()?;
    drop(writer);
    if let Some(temporary_output) = &mut temporary_output {
        temporary_output.commit()?;
    }
    Ok(())
}

fn scan_records(lines: &LinesInFile) -> io::Result<Vec<Record>> {
    let mut records = Vec::new();
    let mut current: Option<Record> = None;
    let mut offset = 0u64;

    lines.try_for_each_line(|line_number, raw_line| {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            if let Some(mut record) = current.take() {
                record.sequence_end = offset;
                records.push(record);
            }
            if line[1..].trim_ascii().is_empty() {
                return Err(line_error(line_number, "record identifier is empty"));
            }
            let sequence_start = offset
                .checked_add(raw_line.len() as u64)
                .ok_or_else(|| line_error(line_number, "file offset overflow"))?;
            current = Some(Record {
                header: line.to_vec(),
                sequence_start,
                sequence_end: 0,
                saw_t: false,
                saw_u: false,
            });
        } else if let Some(record) = current.as_mut() {
            for (index, byte) in line.iter().copied().enumerate() {
                if byte.is_ascii_whitespace() {
                    return Err(line_error(
                        line_number,
                        &format!("whitespace in sequence at column {}", index + 1),
                    ));
                }
                if !is_nucleic_acid(byte) {
                    return Err(line_error(
                        line_number,
                        &format!("invalid nucleotide '{}'", char::from(byte)),
                    ));
                }
                record.saw_t |= matches!(byte, b'T' | b't');
                record.saw_u |= matches!(byte, b'U' | b'u');
                if record.saw_t && record.saw_u {
                    return Err(line_error(
                        line_number,
                        "sequence contains both T and U; use a consistent DNA or RNA alphabet",
                    ));
                }
            }
        } else if !line.is_empty() {
            return Err(line_error(
                line_number,
                "sequence data appears before the first `>` record",
            ));
        }

        offset = offset
            .checked_add(raw_line.len() as u64)
            .ok_or_else(|| line_error(line_number, "file offset overflow"))?;
        Ok(())
    })?;

    if let Some(mut record) = current {
        record.sequence_end = offset;
        records.push(record);
    }
    if records.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no FASTA records found",
        ));
    }
    Ok(records)
}

fn reverse_complement_record(
    input: &mut File,
    record: &Record,
    chars_per_line: usize,
    block: &mut [u8],
    sequence_line: &mut Vec<u8>,
    writer: &mut impl Write,
) -> io::Result<()> {
    writer.write_all(&record.header)?;
    writer.write_all(b"\n")?;

    let is_rna = record.saw_u && !record.saw_t;
    let mut remaining = record.sequence_end.saturating_sub(record.sequence_start);
    while remaining > 0 {
        let block_length = remaining.min(block.len() as u64) as usize;
        let block_start = record.sequence_start + remaining - block_length as u64;
        input.seek(SeekFrom::Start(block_start))?;
        input.read_exact(&mut block[..block_length])?;
        for byte in block[..block_length].iter().rev().copied() {
            if matches!(byte, b'\n' | b'\r') {
                continue;
            }
            sequence_line.push(complement(byte, is_rna));
            if sequence_line.len() == chars_per_line {
                writer.write_all(sequence_line)?;
                writer.write_all(b"\n")?;
                sequence_line.clear();
            }
        }
        remaining -= block_length as u64;
    }
    if !sequence_line.is_empty() {
        writer.write_all(sequence_line)?;
        writer.write_all(b"\n")?;
        sequence_line.clear();
    }
    Ok(())
}

fn complement(byte: u8, is_rna: bool) -> u8 {
    match byte {
        b'A' => {
            if is_rna {
                b'U'
            } else {
                b'T'
            }
        }
        b'a' => {
            if is_rna {
                b'u'
            } else {
                b't'
            }
        }
        b'T' | b'U' => b'A',
        b't' | b'u' => b'a',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'R' => b'Y',
        b'r' => b'y',
        b'Y' => b'R',
        b'y' => b'r',
        b'K' => b'M',
        b'k' => b'm',
        b'M' => b'K',
        b'm' => b'k',
        b'B' => b'V',
        b'b' => b'v',
        b'V' => b'B',
        b'v' => b'b',
        b'D' => b'H',
        b'd' => b'h',
        b'H' => b'D',
        b'h' => b'd',
        b'S' | b'W' | b'N' | b'-' => byte,
        b's' | b'w' | b'n' => byte,
        _ => unreachable!("sequence symbols were validated before transformation"),
    }
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
