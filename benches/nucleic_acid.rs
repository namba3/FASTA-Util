use fasta_util::nucleic_acid::{
    NUCLEIC_ACID_SET, is_nucleic_acid_iter, is_nucleic_acid_lut, is_nucleic_acid_match,
};
use std::{
    env,
    hint::black_box,
    process,
    time::{Duration, Instant},
};

const DEFAULT_INPUT_SIZE: usize = 10_000;
const DEFAULT_SAMPLE_MS: u64 = 200;
const SAMPLE_COUNT: usize = 5;
const WARMUP_ROUNDS: usize = 10;
const INVALID_BASES: &[u8] = b"xyz0123?";

struct Config {
    input_size: usize,
    sample_duration: Duration,
}

fn usage() -> &'static str {
    "Usage: cargo bench --bench nucleic_acid -- [--input-size BYTES] [--sample-ms MS]\n\
     Defaults: --input-size 10000 --sample-ms 200"
}

fn parse_config() -> Result<Option<Config>, String> {
    let mut input_size = DEFAULT_INPUT_SIZE;
    let mut sample_ms = DEFAULT_SAMPLE_MS;
    let mut args = env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input-size" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("missing value for {arg}\n{}", usage()))?;
                input_size = value
                    .parse()
                    .map_err(|_| format!("invalid byte count '{value}' for {arg}\n{}", usage()))?;
                if input_size == 0 {
                    return Err(format!("{arg} must be greater than zero\n{}", usage()));
                }
            }
            "--sample-ms" => {
                let value = args
                    .next()
                    .ok_or_else(|| format!("missing value for {arg}\n{}", usage()))?;
                sample_ms = value
                    .parse()
                    .map_err(|_| format!("invalid duration '{value}' for {arg}\n{}", usage()))?;
                if sample_ms == 0 {
                    return Err(format!("{arg} must be greater than zero\n{}", usage()));
                }
            }
            // Cargo appends this flag when it launches a custom benchmark target.
            "--bench" => {}
            "-h" | "--help" => return Ok(None),
            _ => return Err(format!("unknown argument '{arg}'\n{}", usage())),
        }
    }

    Ok(Some(Config {
        input_size,
        sample_duration: Duration::from_millis(sample_ms),
    }))
}

fn count_valid<F>(sequence: &[u8], predicate: F) -> usize
where
    F: Fn(u8) -> bool,
{
    black_box(sequence)
        .iter()
        .copied()
        .filter(|&base| predicate(base))
        .count()
}

fn measure_sample<F>(
    sequence: &[u8],
    expected_count: usize,
    sample_duration: Duration,
    predicate: F,
) -> f64
where
    F: Fn(u8) -> bool + Copy,
{
    let started = Instant::now();
    let mut rounds = 0u64;
    let mut result = 0;
    while started.elapsed() < sample_duration {
        result = black_box(count_valid(black_box(sequence), predicate));
        rounds += 1;
    }
    let elapsed = started.elapsed();

    assert_eq!(result, expected_count, "benchmark result differs");

    let bases_processed = rounds as f64 * sequence.len() as f64;
    elapsed.as_nanos() as f64 / bases_processed
}

fn benchmark_case(
    case_name: &str,
    sequence: &[u8],
    expected_count: usize,
    sample_duration: Duration,
) {
    for _ in 0..WARMUP_ROUNDS {
        black_box(count_valid(black_box(sequence), is_nucleic_acid_match));
        black_box(count_valid(black_box(sequence), is_nucleic_acid_iter));
        black_box(count_valid(black_box(sequence), is_nucleic_acid_lut));
    }

    let mut samples: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::with_capacity(SAMPLE_COUNT));
    for sample_index in 0..SAMPLE_COUNT {
        let mut order = [0, 1, 2];
        let rotation = sample_index % order.len();
        order.rotate_left(rotation);

        for implementation_index in order {
            let ns_per_base = match implementation_index {
                0 => measure_sample(
                    sequence,
                    expected_count,
                    sample_duration,
                    is_nucleic_acid_match,
                ),
                1 => measure_sample(
                    sequence,
                    expected_count,
                    sample_duration,
                    is_nucleic_acid_iter,
                ),
                2 => measure_sample(
                    sequence,
                    expected_count,
                    sample_duration,
                    is_nucleic_acid_lut,
                ),
                _ => unreachable!(),
            };
            samples[implementation_index].push(ns_per_base);
        }
    }

    println!("\n{case_name}:");
    for (name, sample_values) in ["match", "set iteration", "lookup table"]
        .into_iter()
        .zip(&mut samples)
    {
        sample_values.sort_by(f64::total_cmp);
        let median_ns_per_base = sample_values[SAMPLE_COUNT / 2];
        let million_bases_per_second = 1_000.0 / median_ns_per_base;
        println!(
            "{name:>20}: {median_ns_per_base:>8.3} ns/base, {million_bases_per_second:>8.2} Mbase/s (median of {SAMPLE_COUNT} samples)"
        );
    }
}

fn repeated_bytes(bytes: &[u8], input_size: usize) -> Vec<u8> {
    (0..input_size)
        .map(|index| bytes[index % bytes.len()])
        .collect()
}

fn mixed_sequence(valid_percent: u32, mut state: u32, input_size: usize) -> Vec<u8> {
    (0..input_size)
        .map(|_| {
            // Xorshift32 keeps the generated input deterministic without a dependency.
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;

            if state % 100 < valid_percent {
                NUCLEIC_ACID_SET[(state as usize >> 8) % NUCLEIC_ACID_SET.len()]
            } else {
                INVALID_BASES[(state as usize >> 8) % INVALID_BASES.len()]
            }
        })
        .collect()
}

fn main() {
    let config = match parse_config() {
        Ok(Some(config)) => config,
        Ok(None) => {
            println!("{}", usage());
            return;
        }
        Err(error) => {
            eprintln!("{error}");
            process::exit(2);
        }
    };

    println!(
        "Input: {} bytes per pattern; measuring each implementation for {SAMPLE_COUNT} samples of {} ms",
        config.input_size,
        config.sample_duration.as_millis()
    );

    let cases = [
        (
            "all valid",
            repeated_bytes(NUCLEIC_ACID_SET, config.input_size),
        ),
        (
            "all invalid",
            repeated_bytes(INVALID_BASES, config.input_size),
        ),
        (
            "mixed 50% valid",
            mixed_sequence(50, 0x9e37_79b9, config.input_size),
        ),
        (
            "mostly valid (99%)",
            mixed_sequence(99, 0x243f_6a88, config.input_size),
        ),
    ];

    for (case_name, sequence) in cases {
        let expected_count = sequence
            .iter()
            .filter(|&&base| NUCLEIC_ACID_SET.contains(&base))
            .count();
        benchmark_case(case_name, &sequence, expected_count, config.sample_duration);
    }
}
