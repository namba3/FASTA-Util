use crate::{
    GetArgs, SequenceType, fasta_index,
    output::{InputSource, TemporaryOutput},
    read_lines_from_file, validated_sequence_for,
};
use std::{
    collections::HashSet,
    fs::{self, File},
    io::{self, BufWriter, Write},
    path::Path,
};

struct Query {
    id: Vec<u8>,
    start: usize,
    end_exclusive: Option<usize>,
    region_header: Option<Vec<u8>>,
}

enum Request {
    Record(Query),
    GlobalRange { start: usize, end: usize },
}

pub(super) fn run(args: GetArgs) -> Result<(), Box<dyn std::error::Error>> {
    let input_source = InputSource::from_optional_path(Some(&args.input))?;
    let input_path = input_source.path();
    let is_stdin = args.input == Path::new("-");
    let requests = load_requests(args.ids_file.as_deref(), &args.ids)?;
    if requests.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "provide one or more record IDs or use --ids",
        )
        .into());
    }

    let global_ranges = requests
        .iter()
        .filter(|request| matches!(request, Request::GlobalRange { .. }))
        .count();
    if global_ranges > 0 {
        if requests.len() != 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "a global range cannot be combined with record IDs or another range",
            )
            .into());
        }
        let Request::GlobalRange { start, end } = &requests[0] else {
            unreachable!();
        };
        let index_path = args.fai_index.or_else(|| {
            if is_stdin {
                return None;
            }
            let path = fasta_index::index_path(input_path);
            path.exists().then_some(path)
        });
        return crate::write_global_range(crate::GlobalRangeArgs {
            input: input_path.to_path_buf(),
            sequence_type: args.sequence_type,
            output: args.output,
            fai_index: index_path,
            start: start - 1,
            end_exclusive: *end,
            chars_per_line: args.chars_per_line,
        });
    }
    let queries = requests
        .into_iter()
        .map(|request| match request {
            Request::Record(query) => query,
            Request::GlobalRange { .. } => unreachable!(),
        })
        .collect::<Vec<_>>();

    if let Some(output) = &args.output
        && !is_stdin
    {
        crate::ensure_distinct_input_output(input_path, output)?;
    }

    let index_path = match args.fai_index {
        Some(index_path) => Some(index_path),
        None => {
            if is_stdin {
                None
            } else {
                let index_path = fasta_index::index_path(input_path);
                index_path.exists().then_some(index_path)
            }
        }
    };
    if let (Some(index), Some(output)) = (&index_path, &args.output) {
        crate::ensure_distinct_input_output(index, output)?;
    }

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

    let matched = if let Some(index_path) = &index_path {
        let named_ranges = queries
            .iter()
            .map(|query| fasta_index::NamedRange {
                name: &query.id,
                start: query.start,
                end_exclusive: query.end_exclusive,
                region_header: query.region_header.as_deref(),
            })
            .collect::<Vec<_>>();
        fasta_index::write_named_ranges(
            input_path,
            index_path,
            &named_ranges,
            args.chars_per_line,
            args.sequence_type,
            &mut writer,
        )?
    } else {
        write_from_stream(
            input_path,
            &queries,
            args.chars_per_line,
            args.sequence_type,
            &mut writer,
        )?
    };
    if let Some(index) = matched.iter().position(|matched| !matched) {
        let id = String::from_utf8_lossy(&queries[index].id);
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("record `{id}` was not found"),
        )
        .into());
    }

    writer.flush()?;
    drop(writer);
    if let Some(temporary_output) = &mut temporary_output {
        temporary_output.commit()?;
    }
    Ok(())
}

fn load_requests(ids_file: Option<&Path>, ids: &[String]) -> io::Result<Vec<Request>> {
    let values = if let Some(path) = ids_file {
        let contents = fs::read_to_string(path)?;
        contents
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_owned)
            .collect::<Vec<_>>()
    } else {
        ids.to_vec()
    };
    let requests = values
        .iter()
        .map(|value| parse_request(value))
        .collect::<io::Result<Vec<_>>>()?;
    let mut seen = HashSet::with_capacity(requests.len());
    for query in requests.iter().filter_map(|request| match request {
        Request::Record(query) => Some(query),
        Request::GlobalRange { .. } => None,
    }) {
        if !seen.insert(query.id.as_slice()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "record `{}` is requested more than once",
                    String::from_utf8_lossy(&query.id)
                ),
            ));
        }
    }
    Ok(requests)
}

fn parse_request(value: &str) -> io::Result<Request> {
    if value.contains('-')
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || byte == b'-')
    {
        return parse_global_range(value);
    }
    parse_query(value).map(Request::Record)
}

fn parse_global_range(value: &str) -> io::Result<Request> {
    let Some((start, end)) = value.split_once('-') else {
        return Err(invalid_global_range(value));
    };
    if !start.bytes().all(|byte| byte.is_ascii_digit())
        || !end.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Err(invalid_global_range(value));
    }
    let start = start
        .parse::<usize>()
        .map_err(|_| invalid_global_range(value))?;
    let end = end
        .parse::<usize>()
        .map_err(|_| invalid_global_range(value))?;
    if start == 0 || end < start {
        return Err(invalid_global_range(value));
    }
    Ok(Request::GlobalRange { start, end })
}

fn invalid_global_range(value: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "invalid global range `{value}`; expected START-END with 1-based inclusive coordinates"
        ),
    )
}

fn parse_query(value: &str) -> io::Result<Query> {
    if let Some((id, coordinates)) = value.rsplit_once(':')
        && coordinates
            .bytes()
            .next()
            .is_some_and(|byte| byte.is_ascii_digit() || byte == b'-')
    {
        let Some((start, end)) = coordinates.split_once('-') else {
            return Err(invalid_region(value));
        };
        if !start.bytes().all(|byte| byte.is_ascii_digit())
            || !end.bytes().all(|byte| byte.is_ascii_digit())
        {
            return Err(invalid_region(value));
        }
        let start = start.parse::<usize>().map_err(|_| invalid_region(value))?;
        let end = end.parse::<usize>().map_err(|_| invalid_region(value))?;
        if id.is_empty() || start == 0 || end < start {
            return Err(invalid_region(value));
        }
        let start = start - 1;
        let end_exclusive = end;
        return Ok(Query {
            id: id.as_bytes().to_vec(),
            start,
            end_exclusive: Some(end_exclusive),
            region_header: Some(value.as_bytes().to_vec()),
        });
    }

    Ok(Query {
        id: value.as_bytes().to_vec(),
        start: 0,
        end_exclusive: None,
        region_header: None,
    })
}

fn invalid_region(value: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "invalid region `{value}`; expected ID:START-END with 1-based inclusive coordinates"
        ),
    )
}

fn write_from_stream<W: Write>(
    input_path: &Path,
    queries: &[Query],
    chars_per_line: usize,
    sequence_type: SequenceType,
    writer: &mut W,
) -> io::Result<Vec<bool>> {
    let input = File::open(input_path)?;
    // SAFETY: This command only reads the input; it must not be modified while mapped.
    let lines = unsafe { read_lines_from_file(input)? };
    let mut matched = vec![false; queries.len()];
    let mut active = Vec::new();
    let mut position = 0usize;
    let mut written = vec![0usize; queries.len()];
    let mut saw_header = false;

    lines.try_for_each_line(|_, raw_line| {
        let line = strip_line_ending(raw_line);
        if line.first() == Some(&b'>') {
            saw_header = true;
            finish_active(&active, &written, chars_per_line, writer)?;
            active.clear();
            position = 0;
            let id = line[1..]
                .trim_ascii_start()
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            for (index, query) in queries.iter().enumerate() {
                if query.id != id {
                    continue;
                }
                matched[index] = true;
                active.push(index);
                if let Some(region_header) = query.region_header.as_deref() {
                    writer.write_all(b">")?;
                    writer.write_all(region_header)?;
                } else {
                    writer.write_all(line)?;
                }
                writer.write_all(b"\n")?;
            }
            return Ok(());
        }

        if !saw_header && !line.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sequence data appears before the first `>` record",
            ));
        }
        if active.is_empty() || line.is_empty() {
            return Ok(());
        }
        let sequence = line.trim_ascii_start().trim_ascii_end();
        let line_end = position.checked_add(sequence.len()).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "sequence length overflow")
        })?;
        for index in active.iter().copied() {
            let query = &queries[index];
            let start = position.max(query.start);
            let end = line_end.min(query.end_exclusive.unwrap_or(usize::MAX));
            if start >= end {
                continue;
            }
            let bases =
                validated_sequence_for(&sequence[start - position..end - position], sequence_type)?;
            write_wrapped(writer, bases, chars_per_line, &mut written[index])?;
        }
        position = line_end;
        Ok(())
    })?;
    finish_active(&active, &written, chars_per_line, writer)?;
    Ok(matched)
}

fn finish_active<W: Write>(
    active: &[usize],
    written: &[usize],
    chars_per_line: usize,
    writer: &mut W,
) -> io::Result<()> {
    for index in active {
        if written[*index] > 0 && !written[*index].is_multiple_of(chars_per_line) {
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn write_wrapped<W: Write>(
    writer: &mut W,
    mut bases: &[u8],
    chars_per_line: usize,
    written: &mut usize,
) -> io::Result<()> {
    while !bases.is_empty() {
        let count = (chars_per_line - (*written % chars_per_line)).min(bases.len());
        writer.write_all(&bases[..count])?;
        *written += count;
        bases = &bases[count..];
        if (*written).is_multiple_of(chars_per_line) {
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn strip_line_ending(line: &[u8]) -> &[u8] {
    match line.strip_suffix(b"\n") {
        Some(line) => line.strip_suffix(b"\r").unwrap_or(line),
        None => line,
    }
}

#[cfg(test)]
mod tests {
    use super::{Request, parse_query, parse_request};

    #[test]
    fn parses_ids_and_one_based_inclusive_regions() {
        let Request::Record(id) = parse_request("chr1").unwrap() else {
            panic!("expected record query");
        };
        assert_eq!(id.id, b"chr1");
        assert_eq!(id.start, 0);
        assert_eq!(id.end_exclusive, None);

        let Request::Record(region) = parse_request("chr1:1000-2000").unwrap() else {
            panic!("expected record region");
        };
        assert_eq!(region.id, b"chr1");
        assert_eq!(region.start, 999);
        assert_eq!(region.end_exclusive, Some(2000));

        let Request::GlobalRange { start, end } = parse_request("1000-2000").unwrap() else {
            panic!("expected global range");
        };
        assert_eq!((start, end), (1000, 2000));
    }

    #[test]
    fn rejects_invalid_region_coordinates() {
        assert!(parse_query("chr1:0-3").is_err());
        assert!(parse_query("chr1:4-3").is_err());
        assert!(parse_query("chr1:999999999999999999999-3").is_err());
        assert!(parse_request("0-3").is_err());
        assert!(parse_request("4-3").is_err());
        assert!(parse_request("1-").is_err());
    }
}
