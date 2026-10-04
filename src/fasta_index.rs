use crate::{strip_line_ending, validated_sequence};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write},
    path::Path,
};

const COPY_BUFFER_SIZE: usize = 64 * 1024;
const HEADER_SCAN_BUFFER_SIZE: u64 = 4 * 1024;

struct FaiRecord {
    name: Vec<u8>,
    length: usize,
    sequence_offset: u64,
    line_bases: usize,
    line_width: usize,
}

pub(super) fn write_slice<W: Write>(
    input_path: &Path,
    index_path: &Path,
    start: usize,
    end_exclusive: Option<usize>,
    chars_per_line: usize,
    writer: &mut W,
) -> io::Result<()> {
    let records = read_index(index_path)?;
    let mut input = File::open(input_path)?;
    let file_length = input.metadata()?.len();
    let total_length = validate_index_bounds(&records, file_length)?;
    let range_end = end_exclusive.unwrap_or(total_length).min(total_length);
    let mut record_start = 0usize;
    let mut written = 0usize;

    for record in &records {
        let record_end = record_start + record.length;
        let header = read_header(&mut input, record.sequence_offset, &record.name)?;
        if written > 0 && written % chars_per_line != 0 {
            writer.write_all(b"\n")?;
        }
        writer.write_all(strip_line_ending(&header))?;
        writer.write_all(b"\n")?;

        let selected_start = start.max(record_start).min(record_end);
        let selected_end = range_end.min(record_end);
        if selected_start < selected_end {
            write_sequence_range(
                &mut input,
                record,
                selected_start - record_start,
                selected_end - record_start,
                chars_per_line,
                &mut written,
                writer,
            )?;
        }

        record_start = record_end;
        if record.length > 0 && end_exclusive.is_some_and(|end| end <= record_end) {
            break;
        }
    }

    if end_exclusive.is_some_and(|end| end <= total_length)
        && written > 0
        && written % chars_per_line != 0
    {
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn read_index(path: &Path) -> io::Result<Vec<FaiRecord>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut line = Vec::new();
    let mut line_number = 0usize;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_number += 1;
        let line = strip_line_ending(&line);
        if line.is_empty() {
            return Err(invalid_index_line(line_number, "empty index row"));
        }

        let fields = line.split(|byte| *byte == b'\t').collect::<Vec<_>>();
        if fields.len() < 5 || fields[0].is_empty() {
            return Err(invalid_index_line(
                line_number,
                "expected at least five tab-separated fields",
            ));
        }

        let length = parse_index_number::<usize>(fields[1], line_number, "sequence length")?;
        let sequence_offset = parse_index_number::<u64>(fields[2], line_number, "sequence offset")?;
        let line_bases = parse_index_number::<usize>(fields[3], line_number, "line bases")?;
        let line_width = parse_index_number::<usize>(fields[4], line_number, "line width")?;
        if length > 0 {
            let line_ending_width = line_width.saturating_sub(line_bases);
            let has_invalid_line_ending =
                line_ending_width > 2 || (length > line_bases && line_ending_width == 0);
            if line_bases == 0 || line_width < line_bases || has_invalid_line_ending {
                return Err(invalid_index_line(
                    line_number,
                    "invalid line bases or line width for a FASTA record",
                ));
            }
        }

        records.push(FaiRecord {
            name: fields[0].to_vec(),
            length,
            sequence_offset,
            line_bases,
            line_width,
        });
    }

    Ok(records)
}

fn parse_index_number<T>(field: &[u8], line_number: usize, description: &str) -> io::Result<T>
where
    T: std::str::FromStr,
{
    let value = std::str::from_utf8(field)
        .map_err(|_| invalid_index_line(line_number, &format!("invalid {description}")))?;
    value
        .parse()
        .map_err(|_| invalid_index_line(line_number, &format!("invalid {description}")))
}

fn invalid_index_line(line_number: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("FAI line {line_number}: {message}"),
    )
}

fn validate_index_bounds(records: &[FaiRecord], file_length: u64) -> io::Result<usize> {
    records.iter().try_fold(0usize, |total, record| {
        let record_end = total.checked_add(record.length).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FAI sequence lengths overflow")
        })?;
        if record.length > 0 {
            let last_base_index = record.length - 1;
            let line_width = u64::try_from(record.line_width).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "FAI line width overflow")
            })?;
            let line_offset = u64::try_from(last_base_index / record.line_bases)
                .ok()
                .and_then(|line| line.checked_mul(line_width))
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
                })?;
            let column_offset =
                u64::try_from(last_base_index % record.line_bases).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
                })?;
            let last_base = record
                .sequence_offset
                .checked_add(line_offset)
                .and_then(|offset| offset.checked_add(column_offset))
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
                })?;
            if last_base >= file_length {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FAI sequence extends beyond the input file",
                ));
            }
        } else if record.sequence_offset > file_length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FAI sequence offset is beyond the input file",
            ));
        }
        Ok(record_end)
    })
}

fn read_header(file: &mut File, sequence_offset: u64, expected_name: &[u8]) -> io::Result<Vec<u8>> {
    let header_end = sequence_offset.checked_sub(1).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "FAI sequence offset has no header",
        )
    })?;
    file.seek(SeekFrom::Start(header_end))?;
    let mut terminator = [0];
    file.read_exact(&mut terminator)?;
    if terminator[0] != b'\n' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FAI sequence offset does not follow a FASTA header",
        ));
    }

    let mut search_end = header_end;
    let mut header_start = 0u64;
    let mut block = [0; HEADER_SCAN_BUFFER_SIZE as usize];
    while search_end > 0 {
        let block_start = search_end.saturating_sub(HEADER_SCAN_BUFFER_SIZE);
        let block_length = usize::try_from(search_end - block_start)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "FASTA header is too long"))?;
        file.seek(SeekFrom::Start(block_start))?;
        let chunk = &mut block[..block_length];
        file.read_exact(chunk)?;
        if let Some(newline) = chunk.iter().rposition(|byte| *byte == b'\n') {
            header_start = block_start + newline as u64 + 1;
            break;
        }
        search_end = block_start;
    }

    let header_length = usize::try_from(sequence_offset - header_start)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "FASTA header is too long"))?;
    let mut header = vec![0; header_length];
    file.seek(SeekFrom::Start(header_start))?;
    file.read_exact(&mut header)?;
    let header_body = header
        .strip_prefix(b">")
        .and_then(|header| header.strip_suffix(b"\n"))
        .map(|header| header.strip_suffix(b"\r").unwrap_or(header))
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FAI record header is malformed")
        })?;
    let actual_name = header_body
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    if actual_name != expected_name {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FAI record name does not match the input FASTA header",
        ));
    }
    Ok(header)
}

fn write_sequence_range<W: Write>(
    file: &mut File,
    record: &FaiRecord,
    start: usize,
    end: usize,
    chars_per_line: usize,
    written: &mut usize,
    writer: &mut W,
) -> io::Result<()> {
    let mut position = start;
    let mut raw = Vec::with_capacity(COPY_BUFFER_SIZE);
    while position < end {
        let column = position % record.line_bases;
        let lines_per_chunk = (COPY_BUFFER_SIZE / record.line_width).max(1);
        let bases_per_chunk = lines_per_chunk
            .checked_mul(record.line_bases)
            .unwrap_or(usize::MAX)
            .saturating_sub(column)
            .min(COPY_BUFFER_SIZE);
        let count = (end - position).min(bases_per_chunk);
        let line_index = position / record.line_bases;
        let ending_line_count = column.checked_add(count - 1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
        })? / record.line_bases;
        let skipped_line_bytes = ending_line_count
            .checked_mul(record.line_width - record.line_bases)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
            })?;
        let raw_count = count.checked_add(skipped_line_bytes).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
        })?;
        let raw_count = u64::try_from(raw_count).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
        })?;
        let line_width = u64::try_from(record.line_width)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "FAI line width overflow"))?;
        let column_offset = u64::try_from(column).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
        })?;
        let line_offset = u64::try_from(line_index)
            .ok()
            .and_then(|line| line.checked_mul(line_width))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
            })?;
        let byte_offset = record
            .sequence_offset
            .checked_add(line_offset)
            .and_then(|offset| offset.checked_add(column_offset))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "FAI sequence offset overflow")
            })?;
        let raw_count = usize::try_from(raw_count).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "FAI sequence range is too large",
            )
        })?;
        raw.resize(raw_count, 0);
        file.seek(SeekFrom::Start(byte_offset))?;
        file.read_exact(&mut raw)?;
        // Compact wrapped sequence lines in place and reuse this buffer next iteration.
        let mut raw_position = 0usize;
        let mut compacted_position = 0usize;
        let mut remaining = count;
        let mut line_column = column;
        while remaining > 0 {
            let line_bases = (record.line_bases - line_column).min(remaining);
            raw.copy_within(raw_position..raw_position + line_bases, compacted_position);
            raw_position += line_bases;
            compacted_position += line_bases;
            remaining -= line_bases;
            if remaining > 0 {
                raw_position += record.line_width - record.line_bases;
                line_column = 0;
            }
        }
        raw.truncate(count);
        let bases = validated_sequence(&raw)?;
        if bases.len() != count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FAI sequence layout does not match the input FASTA",
            ));
        }
        write_wrapped_bases(writer, bases, chars_per_line, written)?;
        position += count;
    }
    Ok(())
}

fn write_wrapped_bases<W: Write>(
    writer: &mut W,
    mut bases: &[u8],
    chars_per_line: usize,
    written: &mut usize,
) -> io::Result<()> {
    while !bases.is_empty() {
        let line_remaining = chars_per_line - (*written % chars_per_line);
        let count = line_remaining.min(bases.len());
        writer.write_all(&bases[..count])?;
        *written += count;
        bases = &bases[count..];
        if *written % chars_per_line == 0 {
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}
