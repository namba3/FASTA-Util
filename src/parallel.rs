use std::{io, num::NonZeroUsize, ops::Range, thread};

const AUTO_PARALLEL_THRESHOLD: usize = 8 * 1024 * 1024;
const MAX_AUTO_THREADS: usize = 4;
const MAX_REQUESTED_THREADS: usize = 64;

#[derive(Clone)]
pub(crate) struct LineChunk<'a> {
    pub(crate) range: Range<usize>,
    pub(crate) start_line: usize,
    pub(crate) bytes: &'a [u8],
}

pub(crate) fn worker_count(file_size: usize, requested: Option<NonZeroUsize>) -> usize {
    if let Some(requested) = requested {
        return requested.get().min(MAX_REQUESTED_THREADS);
    }
    if file_size < AUTO_PARALLEL_THRESHOLD {
        return 1;
    }
    thread::available_parallelism()
        .map(NonZeroUsize::get)
        .unwrap_or(1)
        .min(MAX_AUTO_THREADS)
}

pub(crate) fn line_chunks(bytes: &[u8], workers: usize) -> Vec<LineChunk<'_>> {
    if bytes.is_empty() {
        return vec![LineChunk {
            range: 0..0,
            start_line: 1,
            bytes,
        }];
    }

    let workers = workers.max(1).min(bytes.len());
    let mut chunks = Vec::with_capacity(workers);
    let mut start = 0;
    let mut start_line = 1;
    let quotient = bytes.len() / workers;
    let remainder = bytes.len() % workers;

    for index in 1..workers {
        let target = quotient * index + remainder * index / workers;
        if target <= start || target >= bytes.len() {
            continue;
        }
        let boundary = bytes[target..]
            .iter()
            .position(|byte| *byte == b'\n')
            .map_or(bytes.len(), |offset| target + offset + 1);
        if boundary <= start || boundary >= bytes.len() {
            continue;
        }

        let range = start..boundary;
        let chunk = &bytes[range.clone()];
        chunks.push(LineChunk {
            range,
            start_line,
            bytes: chunk,
        });
        start_line += chunk.iter().filter(|byte| **byte == b'\n').count();
        start = boundary;
    }

    let range = start..bytes.len();
    chunks.push(LineChunk {
        bytes: &bytes[range.clone()],
        range,
        start_line,
    });
    chunks
}

/// Splits input at FASTA record headers so a record is parsed by one worker.
pub(crate) fn record_chunks(bytes: &[u8], workers: usize) -> Vec<LineChunk<'_>> {
    if bytes.is_empty() {
        return line_chunks(bytes, 1);
    }

    let mut headers = Vec::new();
    let mut offset = 0;
    let mut line_number = 1;
    for line in bytes.split_inclusive(|byte| *byte == b'\n') {
        if line.first() == Some(&b'>') {
            headers.push((offset, line_number));
        }
        offset += line.len();
        line_number += line.iter().filter(|byte| **byte == b'\n').count();
    }
    if headers.len() < 2 || workers <= 1 {
        return line_chunks(bytes, 1);
    }

    let workers = workers.min(headers.len());
    let mut chunks = Vec::with_capacity(workers);
    let mut start = 0;
    let mut start_line = 1;
    for index in 1..workers {
        let target = bytes.len() / workers * index;
        if let Some(&(boundary, boundary_line)) = headers
            .iter()
            .find(|(header_offset, _)| *header_offset > start && *header_offset >= target)
        {
            chunks.push(LineChunk {
                range: start..boundary,
                start_line,
                bytes: &bytes[start..boundary],
            });
            start = boundary;
            start_line = boundary_line;
        }
    }
    chunks.push(LineChunk {
        range: start..bytes.len(),
        start_line,
        bytes: &bytes[start..],
    });
    chunks
}

pub(crate) fn first_header_offset(bytes: &[u8]) -> Option<usize> {
    let mut offset = 0;
    for line in bytes.split_inclusive(|byte| *byte == b'\n') {
        if line.first() == Some(&b'>') {
            return Some(offset);
        }
        offset += line.len();
    }
    None
}

pub(crate) fn process_chunks<T: Send>(
    chunks: &[LineChunk<'_>],
    process: impl Fn(LineChunk<'_>) -> io::Result<T> + Sync + Send,
) -> Vec<io::Result<T>> {
    if chunks.len() == 1 {
        return vec![process(chunks[0].clone())];
    }

    thread::scope(|scope| {
        let process = &process;
        let handles = chunks
            .iter()
            .cloned()
            .map(|chunk| scope.spawn(move || process(chunk)))
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .map(|handle| {
                handle
                    .join()
                    .unwrap_or_else(|_| Err(io::Error::other("FASTA worker thread panicked")))
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::{first_header_offset, line_chunks, record_chunks, worker_count};
    use std::num::NonZeroUsize;

    #[test]
    fn line_chunks_cover_input_once_and_keep_line_boundaries() {
        let input = b">one\r\nACGT\n>two\nTGCA\n>three\nNN\n";
        let chunks = line_chunks(input, 3);

        assert_eq!(chunks.first().unwrap().range.start, 0);
        assert_eq!(chunks.last().unwrap().range.end, input.len());
        assert!(chunks.windows(2).all(|pair| {
            pair[0].range.end == pair[1].range.start
                && (pair[0].range.end == input.len() || input[pair[0].range.end - 1] == b'\n')
        }));
        assert_eq!(
            chunks
                .iter()
                .flat_map(|chunk| chunk.bytes.split_inclusive(|byte| *byte == b'\n'))
                .collect::<Vec<_>>(),
            input
                .split_inclusive(|byte| *byte == b'\n')
                .collect::<Vec<_>>()
        );
        assert_eq!(
            chunks
                .iter()
                .map(|chunk| chunk.start_line)
                .collect::<Vec<_>>(),
            [1, 3, 5]
        );
    }

    #[test]
    fn empty_and_unterminated_inputs_have_valid_chunks() {
        let empty = line_chunks(b"", 4);
        let unterminated = line_chunks(b">record\nACGT", 4);

        assert_eq!(empty.len(), 1);
        assert_eq!(empty[0].start_line, 1);
        assert_eq!(unterminated.last().unwrap().bytes, b"ACGT");
    }

    #[test]
    fn finds_first_header_only_at_a_line_start() {
        assert_eq!(first_header_offset(b"\n>record\nACGT\n"), Some(1));
        assert_eq!(first_header_offset(b"AC>GT\n"), None);
    }

    #[test]
    fn record_chunks_keep_complete_records_and_global_line_numbers() {
        let input = b"\n>one\nAC\n>two\r\nGT\r\n>three\nNN";
        let chunks = record_chunks(input, 3);

        assert_eq!(chunks.first().unwrap().start_line, 1);
        assert_eq!(
            chunks
                .iter()
                .map(|chunk| chunk.start_line)
                .collect::<Vec<_>>(),
            [1, 4, 6]
        );
        assert!(chunks.iter().all(|chunk| {
            chunk
                .bytes
                .split_inclusive(|byte| *byte == b'\n')
                .filter(|line| line.first() == Some(&b'>'))
                .count()
                >= 1
        }));
        assert_eq!(
            chunks.iter().map(|chunk| chunk.bytes.len()).sum::<usize>(),
            input.len()
        );
    }

    #[test]
    fn automatic_worker_count_keeps_small_files_sequential() {
        assert_eq!(worker_count(1024, None), 1);
        assert_eq!(worker_count(1024, Some(NonZeroUsize::new(3).unwrap())), 3);
        assert_eq!(
            worker_count(20 * 1024 * 1024, Some(NonZeroUsize::new(99).unwrap())),
            64
        );
    }
}
