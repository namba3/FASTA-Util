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

fn main() {
    let sequence = (0..INPUT_SIZE)
        .map(|index| NUCLEIC_ACID_SET[index % NUCLEIC_ACID_SET.len()])
        .collect::<Vec<_>>();
    let expected_count = count_valid(&sequence, is_nucleic_acid_match);

    println!(
        "Input: {INPUT_SIZE} bases; measuring each implementation for {} second(s)",
        MEASURE_FOR.as_secs()
    );
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
