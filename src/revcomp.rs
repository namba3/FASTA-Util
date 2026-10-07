use crate::{
    RevcompArgs, ensure_distinct_input_output, is_nucleic_acid, line_error,
    output::{InputSource, with_output},
    read_lines_from_file, strip_line_ending,
};
use fasta_util::LinesInFile;
use std::{
    fs::File,
    io::{self, Read, Seek, SeekFrom, Write},
};

const REVERSE_READ_SIZE: usize = 64 * 1024;
const DNA_COMPLEMENT: [u8; 256] = complement_table(false);
const RNA_COMPLEMENT: [u8; 256] = complement_table(true);

const fn complement_table(is_rna: bool) -> [u8; 256] {
    let mut table = [0; 256];
    let mut byte = 0;
    while byte < table.len() {
        table[byte] = byte as u8;
        byte += 1;
    }

    table[b'A' as usize] = if is_rna { b'U' } else { b'T' };
    table[b'a' as usize] = if is_rna { b'u' } else { b't' };
    table[b'T' as usize] = b'A';
    table[b'U' as usize] = b'A';
    table[b't' as usize] = b'a';
    table[b'u' as usize] = b'a';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'g';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'c';
    table[b'R' as usize] = b'Y';
    table[b'r' as usize] = b'y';
    table[b'Y' as usize] = b'R';
    table[b'y' as usize] = b'r';
    table[b'K' as usize] = b'M';
    table[b'k' as usize] = b'm';
    table[b'M' as usize] = b'K';
    table[b'm' as usize] = b'k';
    table[b'B' as usize] = b'V';
    table[b'b' as usize] = b'v';
    table[b'V' as usize] = b'B';
    table[b'v' as usize] = b'b';
    table[b'D' as usize] = b'H';
    table[b'd' as usize] = b'h';
    table[b'H' as usize] = b'D';
    table[b'h' as usize] = b'd';
    table
}

struct Record {
    header: Vec<u8>,
    sequence_start: u64,
    sequence_end: u64,
    saw_t: bool,
    saw_u: bool,
}

pub(super) fn run(args: RevcompArgs) -> Result<(), Box<dyn std::error::Error>> {
    if let (Some(input), Some(output)) = (&args.input, &args.output)
        && input != std::path::Path::new("-")
    {
        ensure_distinct_input_output(input, output)?;
    }

    let input_source = InputSource::from_optional_path(args.input.as_deref())?;
    let file = File::open(input_source.path())?;
    // SAFETY: The input file must not change while its memory map is alive.
    let lines = unsafe { read_lines_from_file(file)? };
    let records = scan_records(&lines)?;

    let mut input = File::open(input_source.path())?;
    with_output(args.output.as_deref(), |writer| {
        let mut block = vec![0; REVERSE_READ_SIZE];
        let mut sequence_line = Vec::new();
        for record in &records {
            reverse_complement_record(
                &mut input,
                record,
                args.chars_per_line,
                &mut block,
                &mut sequence_line,
                writer,
            )?;
        }
        Ok(())
    })?;
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
    if is_rna {
        RNA_COMPLEMENT[byte as usize]
    } else {
        DNA_COMPLEMENT[byte as usize]
    }
}
