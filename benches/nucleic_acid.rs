use fasta_util::nucleic_acid::{
    NUCLEIC_ACID_SET, is_nucleic_acid_iter, is_nucleic_acid_lut, is_nucleic_acid_match,
};
use std::{
    hint::black_box,
    time::{Duration, Instant},
};

const INPUT_SIZE: usize = 10_000;
const MEASURE_FOR: Duration = Duration::from_secs(1);
const WARMUP_ROUNDS: usize = 10;
const INVALID_BASES: &[u8] = b"xyz0123?";

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

fn benchmark<F>(name: &str, sequence: &[u8], expected_count: usize, predicate: F)
where
    F: Fn(u8) -> bool,
{
    for _ in 0..WARMUP_ROUNDS {
        black_box(count_valid(black_box(sequence), &predicate));
    }

    let started = Instant::now();
    let mut rounds = 0u64;
    let mut result = 0;
    while started.elapsed() < MEASURE_FOR {
        result = black_box(count_valid(black_box(sequence), &predicate));
        rounds += 1;
    }
    let elapsed = started.elapsed();

    assert_eq!(result, expected_count, "{name} returned a different result");

    let bases_processed = rounds as f64 * sequence.len() as f64;
    let ns_per_base = elapsed.as_nanos() as f64 / bases_processed;
    let million_bases_per_second = bases_processed / elapsed.as_secs_f64() / 1_000_000.0;

    println!(
        "{name:>20}: {ns_per_base:>8.3} ns/base, {million_bases_per_second:>8.2} Mbase/s ({rounds} rounds)"
    );
}

fn repeated_bytes(bytes: &[u8]) -> Vec<u8> {
    (0..INPUT_SIZE)
        .map(|index| bytes[index % bytes.len()])
        .collect()
}

fn mixed_sequence(valid_percent: u32, mut state: u32) -> Vec<u8> {
    (0..INPUT_SIZE)
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
    println!(
        "Input: {INPUT_SIZE} bytes per pattern; measuring each implementation for {} second(s)",
        MEASURE_FOR.as_secs()
    );

    let cases = [
        ("all valid", repeated_bytes(NUCLEIC_ACID_SET)),
        ("all invalid", repeated_bytes(INVALID_BASES)),
        ("mixed 50% valid", mixed_sequence(50, 0x9e37_79b9)),
        ("mostly valid (99%)", mixed_sequence(99, 0x243f_6a88)),
    ];

    for (case_name, sequence) in cases {
        let expected_count = sequence
            .iter()
            .filter(|&&base| NUCLEIC_ACID_SET.contains(&base))
            .count();
        println!("\n{case_name}:");
        benchmark("match", &sequence, expected_count, is_nucleic_acid_match);
        benchmark(
            "set iteration",
            &sequence,
            expected_count,
            is_nucleic_acid_iter,
        );
        benchmark(
            "lookup table",
            &sequence,
            expected_count,
            is_nucleic_acid_lut,
        );
    }
}
