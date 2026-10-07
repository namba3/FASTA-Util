use crate::{SequenceType, output::with_output, strip_line_ending, validated_sequence_for};
use std::{
    collections::HashSet,
    ffi::OsString,
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

const COPY_BUFFER_SIZE: usize = 64 * 1024;
const HEADER_SCAN_BUFFER_SIZE: u64 = 4 * 1024;

#[derive(Clone, Copy, PartialEq, Eq)]
enum FaiLineEnding {
    Lf,
    CrLf,
}

struct IndexRecord {
    name: Vec<u8>,
    sequence_offset: u64,
    length: usize,
    line_bases: Option<usize>,
    line_width: Option<usize>,
    line_ending: Option<FaiLineEnding>,
    pending_short_line: Option<usize>,
}

pub(super) fn index_path(input_path: &Path) -> PathBuf {
    let mut path = OsString::from(input_path.as_os_str());
    path.push(".fai");
    PathBuf::from(path)
}

pub(super) fn create_index(input_path: &Path, index_path: &Path) -> io::Result<usize> {
    let input = File::open(input_path)?;
    with_output(Some(index_path), |writer| write_index(input, writer))
}

fn write_index<W: Write>(input: File, writer: &mut W) -> io::Result<usize> {
    let mut reader = BufReader::new(input);
    let mut line = Vec::new();
    let mut byte_offset = 0u64;
    let mut line_number = 0usize;
    let mut records = 0usize;
    let mut current_record: Option<IndexRecord> = None;
    let mut names = HashSet::<Vec<u8>>::new();

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_number += 1;
        let raw_line_length = u64::try_from(line.len()).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, "FASTA byte offset overflow")
        })?;
        let next_byte_offset = byte_offset.checked_add(raw_line_length).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FASTA byte offset overflow")
        })?;
        let (content, ending) = strip_index_line_ending(&line);

        if content.first() == Some(&b'>') {
            if let Some(record) = current_record.take() {
                write_index_record(record, writer)?;
            }

            let name = content[1..]
                .trim_ascii_start()
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            if name.is_empty() {
                return Err(index_error(line_number, "record identifier is empty"));
            }
            if !names.insert(name.to_vec()) {
                return Err(index_error(
                    line_number,
                    &format!(
                        "duplicate record identifier `{}`",
                        String::from_utf8_lossy(name)
                    ),
                ));
            }
            current_record = Some(IndexRecord {
                name: name.to_vec(),
                sequence_offset: next_byte_offset,
                length: 0,
                line_bases: None,
                line_width: None,
                line_ending: None,
                pending_short_line: None,
            });
            records += 1;
            byte_offset = next_byte_offset;
            continue;
        }

        let Some(record) = current_record.as_mut() else {
            if !content.is_empty() {
                return Err(index_error(
                    line_number,
                    "sequence data appears before the first `>` record",
                ));
            }
            return Err(index_error(
                line_number,
                "blank line appears outside a FASTA record",
            ));
        };

        if content.is_empty() {
            return Err(index_error(
                line_number,
                "blank sequence line cannot be represented by `.fai`",
            ));
        }
        if let Some(column) = content.iter().position(u8::is_ascii_whitespace) {
            return Err(index_error(
                line_number,
                &format!("whitespace in sequence line at column {}", column + 1),
            ));
        }

        if let Some(previous_short_line) = record.pending_short_line.take() {
            return Err(index_error(
                previous_short_line,
                "non-final sequence line is shorter than the `.fai` line width",
            ));
        }

        match record.line_ending {
            Some(previous) if ending.is_some_and(|ending| ending != previous) => {
                return Err(index_error(
                    line_number,
                    "line endings must be consistent within each FASTA record",
                ));
            }
            None => record.line_ending = ending,
            _ => {}
        }

        if let Some(line_bases) = record.line_bases {
            if content.len() < line_bases {
                record.pending_short_line = Some(line_number);
            } else if content.len() > line_bases {
                return Err(index_error(
                    line_number,
                    "sequence line is wider than the first line",
                ));
            }
        } else {
            record.line_bases = Some(content.len());
            record.line_width = Some(line.len());
        }
        record.length = record
            .length
            .checked_add(content.len())
            .ok_or_else(|| index_error(line_number, "sequence length overflow"))?;
        byte_offset = next_byte_offset;
    }

    if let Some(record) = current_record {
        write_index_record(record, writer)?;
    }
    if records == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no FASTA records found",
        ));
    }
    Ok(records)
}

fn write_index_record<W: Write>(record: IndexRecord, writer: &mut W) -> io::Result<()> {
    let line_bases = record.line_bases.unwrap_or(0);
    let line_width = record.line_width.unwrap_or(0);
    writer.write_all(&record.name)?;
    writeln!(
        writer,
        "\t{}\t{}\t{}\t{}",
        record.length, record.sequence_offset, line_bases, line_width
    )
}

fn index_error(line_number: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("line {line_number}: {message}"),
    )
}

fn strip_index_line_ending(line: &[u8]) -> (&[u8], Option<FaiLineEnding>) {
    if let Some(without_lf) = line.strip_suffix(b"\n") {
        if let Some(without_cr) = without_lf.strip_suffix(b"\r") {
            (without_cr, Some(FaiLineEnding::CrLf))
        } else {
            (without_lf, Some(FaiLineEnding::Lf))
        }
    } else {
        (line, None)
    }
}

struct FaiRecord {
    name: Vec<u8>,
    length: usize,
    sequence_offset: u64,
    line_bases: usize,
    line_width: usize,
}

#[derive(Clone, Copy)]
struct SequenceOutputOptions {
    chars_per_line: usize,
    sequence_type: SequenceType,
}

pub(super) struct NamedRange<'a> {
    pub(super) name: &'a [u8],
    pub(super) start: usize,
    pub(super) end_exclusive: Option<usize>,
    pub(super) region_header: Option<&'a [u8]>,
}

pub(super) fn write_named_ranges<W: Write>(
    input_path: &Path,
    index_path: &Path,
    queries: &[NamedRange<'_>],
    chars_per_line: usize,
    sequence_type: SequenceType,
    writer: &mut W,
) -> io::Result<Vec<bool>> {
    let records = read_index(index_path)?;
    let mut input = File::open(input_path)?;
    let file_length = input.metadata()?.len();
    validate_index_bounds(&records, file_length)?;
    let options = SequenceOutputOptions {
        chars_per_line,
        sequence_type,
    };
    let mut matched = vec![false; queries.len()];

    for record in &records {
        let matching_queries = queries
            .iter()
            .enumerate()
            .filter(|(_, query)| query.name == record.name);
        for (query_index, query) in matching_queries {
            matched[query_index] = true;
            let source_header = read_header(&mut input, record.sequence_offset, &record.name)?;
            let header = match query.region_header {
                Some(region_header) => {
                    let mut output_header = Vec::with_capacity(region_header.len() + 2);
                    output_header.push(b'>');
                    output_header.extend_from_slice(region_header);
                    output_header.push(b'\n');
                    output_header
                }
                None => source_header,
            };
            if header.first() != Some(&b'>') {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "FAI record header is malformed",
                ));
            }
            writer.write_all(strip_line_ending(&header))?;
            writer.write_all(b"\n")?;

            let start = query.start.min(record.length);
            let end = query
                .end_exclusive
                .unwrap_or(record.length)
                .min(record.length);
            let mut written = 0usize;
            if start < end {
                write_sequence_range(
                    &mut input,
                    record,
                    start,
                    end,
                    options,
                    &mut written,
                    writer,
                )?;
            }
            if written > 0 && !written.is_multiple_of(chars_per_line) {
                writer.write_all(b"\n")?;
            }
        }
    }
    Ok(matched)
}

pub(super) fn write_slice<W: Write>(
    input_path: &Path,
    index_path: &Path,
    start: usize,
    end_exclusive: Option<usize>,
    chars_per_line: usize,
    sequence_type: SequenceType,
    writer: &mut W,
) -> io::Result<()> {
    let records = read_index(index_path)?;
    let mut input = File::open(input_path)?;
    let file_length = input.metadata()?.len();
    let total_length = validate_index_bounds(&records, file_length)?;
    let range_end = end_exclusive.unwrap_or(total_length).min(total_length);
    let sequence_options = SequenceOutputOptions {
        chars_per_line,
        sequence_type,
    };
    let mut record_start = 0usize;
    let mut written = 0usize;

    for record in &records {
        let record_end = record_start + record.length;
        let header = read_header(&mut input, record.sequence_offset, &record.name)?;
        if written > 0 && !written.is_multiple_of(chars_per_line) {
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
                sequence_options,
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
        && !written.is_multiple_of(chars_per_line)
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

        let mut fields = line.split(|byte| *byte == b'\t');
        let (Some(name), Some(length), Some(sequence_offset), Some(line_bases), Some(line_width)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            return Err(invalid_index_line(
                line_number,
                "expected at least five tab-separated fields",
            ));
        };
        if name.is_empty() {
            return Err(invalid_index_line(
                line_number,
                "expected at least five tab-separated fields",
            ));
        }

        let length = parse_index_number::<usize>(length, line_number, "sequence length")?;
        let sequence_offset =
            parse_index_number::<u64>(sequence_offset, line_number, "sequence offset")?;
        let line_bases = parse_index_number::<usize>(line_bases, line_number, "line bases")?;
        let line_width = parse_index_number::<usize>(line_width, line_number, "line width")?;
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
            name: name.to_vec(),
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
        .trim_ascii_start()
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
    options: SequenceOutputOptions,
    written: &mut usize,
    writer: &mut W,
) -> io::Result<()> {
    let mut position = start;
    let mut raw = Vec::with_capacity(COPY_BUFFER_SIZE);
    while position < end {
        let column = position % record.line_bases;
        let lines_per_chunk = (COPY_BUFFER_SIZE / record.line_width).max(1);
        let bases_per_chunk = lines_per_chunk
            .saturating_mul(record.line_bases)
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
        let bases = validated_sequence_for(&raw, options.sequence_type)?;
        if bases.len() != count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FAI sequence layout does not match the input FASTA",
            ));
        }
        write_wrapped_bases(writer, bases, options.chars_per_line, written)?;
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
        if (*written).is_multiple_of(chars_per_line) {
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}
