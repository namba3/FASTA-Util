mod support;

use fasta_util::{LinesInFile, read_lines_from_file};
use std::{
    env,
    fs::{self, File, OpenOptions},
    hint::black_box,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process,
    sync::atomic::{AtomicUsize, Ordering},
    time::{Duration, Instant},
};
use support::{parse_config_for, usage_for};

const LINE_WIDTH: usize = 60;
const RECORD_BASES: usize = 1_000;
const SAMPLE_COUNT: usize = 5;
const WARMUP_ROUNDS: usize = 3;
const BASES: &[u8] = b"ACGTN";

static NEXT_TEMP_FILE_ID: AtomicUsize = AtomicUsize::new(0);

struct TemporaryFasta(PathBuf);

impl TemporaryFasta {
    fn new(sequence_bases: usize, line_ending: &[u8]) -> io::Result<Self> {
        loop {
            let id = NEXT_TEMP_FILE_ID.fetch_add(1, Ordering::Relaxed);
            let path = env::temp_dir().join(format!("fasta-util-bench-{}-{id}.fa", process::id()));
            match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(mut file) => {
                    let temporary = Self(path);
                    let mut writer = BufWriter::new(&mut file);
                    write_fasta(&mut writer, sequence_bases, line_ending)?;
                    writer.flush()?;
                    return Ok(temporary);
                }
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TemporaryFasta {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ScanResult {
    lines: usize,
    bytes: usize,
}

fn write_fasta(
    writer: &mut impl Write,
    sequence_bases: usize,
    line_ending: &[u8],
) -> io::Result<()> {
    let mut remaining = sequence_bases;
    let mut record_number = 0;
    let mut base_index = 0;

    while remaining > 0 {
        record_number += 1;
        let record_length = remaining.min(RECORD_BASES);
        writeln!(writer, ">record_{record_number} length={record_length}")?;
        let mut record_remaining = record_length;
        while record_remaining > 0 {
            let line_length = record_remaining.min(LINE_WIDTH);
            let mut line = Vec::with_capacity(line_length);
            for _ in 0..line_length {
                line.push(BASES[base_index % BASES.len()]);
                base_index += 1;
            }
            writer.write_all(&line)?;
            writer.write_all(line_ending)?;
            record_remaining -= line_length;
        }
        remaining -= record_length;
    }

    Ok(())
}

fn scan_mmap_visitor(path: &Path) -> io::Result<ScanResult> {
    let file = File::open(path)?;
    // SAFETY: The benchmark does not modify the input while the mapping is in use.
    let lines = unsafe { read_lines_from_file(file)? };
    let mut result = ScanResult { lines: 0, bytes: 0 };
    lines.try_for_each_line(|_, line| {
        result.lines += 1;
        result.bytes += line.len();
        Ok::<(), io::Error>(())
    })?;
    Ok(result)
}

fn scan_mmap_iterator(path: &Path) -> io::Result<ScanResult> {
    let file = File::open(path)?;
    // SAFETY: The benchmark does not modify the input while the mapping is in use.
    let lines = unsafe { read_lines_from_file(file)? };
    Ok(scan_lines(lines))
}

fn scan_lines(lines: LinesInFile) -> ScanResult {
    let mut result = ScanResult { lines: 0, bytes: 0 };
    for line in lines {
        result.lines += 1;
        result.bytes += line.as_ref().len();
    }
    result
}

fn scan_buffered(path: &Path) -> io::Result<ScanResult> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let mut line = Vec::with_capacity(LINE_WIDTH + 2);
    let mut result = ScanResult { lines: 0, bytes: 0 };
    loop {
        line.clear();
        let read = reader.read_until(b'\n', &mut line)?;
        if read == 0 {
            break;
        }
        result.lines += 1;
        result.bytes += read;
    }
    Ok(result)
}

type Scanner = fn(&Path) -> io::Result<ScanResult>;

fn measure_sample(
    path: &Path,
    expected: ScanResult,
    sample_duration: Duration,
    scanner: Scanner,
) -> io::Result<f64> {
    let started = Instant::now();
    let mut rounds = 0u64;
    let mut result = ScanResult { lines: 0, bytes: 0 };
    while started.elapsed() < sample_duration {
        result = black_box(scanner(black_box(path))?);
        rounds += 1;
    }
    let elapsed = started.elapsed();
    assert_eq!(result, expected, "benchmark scan result differs");

    let total_bytes = rounds as f64 * expected.bytes as f64;
    Ok(total_bytes / elapsed.as_secs_f64() / (1024.0 * 1024.0))
}

fn benchmark_line_endings(
    name: &str,
    sequence_bases: usize,
    line_ending: &[u8],
    sample_duration: Duration,
) -> io::Result<()> {
    let input = TemporaryFasta::new(sequence_bases, line_ending)?;
    let path = input.path();
    let expected = scan_mmap_visitor(path)?;
    assert_eq!(scan_mmap_iterator(path)?, expected);
    assert_eq!(scan_buffered(path)?, expected);

    let scanners: [(&str, Scanner); 3] = [
        ("mmap borrowed visitor", scan_mmap_visitor),
        ("mmap line iterator", scan_mmap_iterator),
        ("buffered read_until", scan_buffered),
    ];
    let mut samples: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::with_capacity(SAMPLE_COUNT));

    for _ in 0..WARMUP_ROUNDS {
        for &(_, scanner) in &scanners {
            assert_eq!(scanner(path)?, expected);
        }
    }

    for sample_index in 0..SAMPLE_COUNT {
        let mut order = [0, 1, 2];
        let rotation = sample_index % order.len();
        order.rotate_left(rotation);
        for scanner_index in order {
            let throughput =
                measure_sample(path, expected, sample_duration, scanners[scanner_index].1)?;
            samples[scanner_index].push(throughput);
        }
    }

    println!(
        "\n{name}: {} sequence bases, {} FASTA bytes, {} lines",
        sequence_bases, expected.bytes, expected.lines
    );
    for ((scanner_name, _), sample_values) in scanners.into_iter().zip(&mut samples) {
        sample_values.sort_by(f64::total_cmp);
        println!(
            "{scanner_name:>24}: {:>9.2} MiB/s (median of {SAMPLE_COUNT} samples)",
            sample_values[SAMPLE_COUNT / 2]
        );
    }
    Ok(())
}

fn main() -> io::Result<()> {
    let config = match parse_config_for(env::args().skip(1), "fasta_io") {
        Ok(Some(config)) => config,
        Ok(None) => {
            println!("{}", usage_for("fasta_io"));
            return Ok(());
        }
        Err(error) => {
            eprintln!("{error}");
            process::exit(2);
        }
    };

    println!(
        "FASTA warm-cache line scanning: {} sequence bases; {} samples of {} ms after {WARMUP_ROUNDS} warmup rounds",
        config.input_size,
        SAMPLE_COUNT,
        config.sample_duration.as_millis()
    );
    benchmark_line_endings("LF", config.input_size, b"\n", config.sample_duration)?;
    benchmark_line_endings("CRLF", config.input_size, b"\r\n", config.sample_duration)?;
    Ok(())
}
