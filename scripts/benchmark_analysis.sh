#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s [BASES] [RUNS]\n' "${0##*/}"
    printf 'Defaults: BASES=1000000 RUNS=5; BASES must be at least 64\n'
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if [[ $# -gt 2 ]]; then
    usage >&2
    exit 2
fi

bases="${1:-1000000}"
runs="${2:-5}"
if [[ ! "$bases" =~ ^[1-9][0-9]*$ || ! "$runs" =~ ^[1-9][0-9]*$ || "$bases" -lt 64 ]]; then
    printf 'BASES must be an integer of at least 64 and RUNS must be a positive integer\n' >&2
    usage >&2
    exit 2
fi

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
rustup_bin="$(command -v rustup || true)"
hyperfine="$(command -v hyperfine || true)"
awk_bin="$(command -v awk || true)"
if [[ -z "$rustup_bin" || ! -x "$rustup_bin" ]]; then
    printf 'Required command not found: rustup\n' >&2
    exit 1
fi
if [[ -z "$hyperfine" || ! -x "$hyperfine" ]]; then
    printf 'Required command not found: hyperfine\n' >&2
    exit 1
fi
if [[ -z "$awk_bin" || ! -x "$awk_bin" ]]; then
    printf 'Required command not found: awk\n' >&2
    exit 1
fi

"$rustup_bin" run stable cargo build --release --locked --manifest-path "$repo_root/Cargo.toml" --bin fasta-util

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
binary="$repo_root/target/release/fasta-util"
random_fasta="$work_dir/random.fa"
locate_fasta="$work_dir/homopolymer.fa"

"$rustup_bin" run stable cargo run --quiet --locked --manifest-path "$repo_root/generate_random_data/Cargo.toml" -- "$bases" --seed 42 --line-width 50 > "$random_fasta"

make_stats_input() {
    local lines_per_record="$1"
    local output="$2"
    "$awk_bin" -v lines_per_record="$lines_per_record" '
        /^>/ { next }
        {
            if (line_count % lines_per_record == 0) {
                printf ">record_%d\n", record_count + 1
                record_count++
            }
            print
            line_count++
        }
    ' "$random_fasta" > "$output"
}

stats_short_records="$work_dir/stats-100bp.fa"
stats_long_records="$work_dir/stats-10000bp.fa"
make_stats_input 2 "$stats_short_records"
make_stats_input 200 "$stats_long_records"

"$awk_bin" -v total_bases="$bases" '
    BEGIN {
        record_size = 10000
        record_number = 0
        remaining = total_bases
        while (remaining > 0) {
            record_number++
            printf ">homopolymer_%d\n", record_number
            record_bases = remaining < record_size ? remaining : record_size
            while (record_bases > 0) {
                line_bases = record_bases < 60 ? record_bases : 60
                line = ""
                for (i = 0; i < line_bases; i++) line = line "A"
                print line
                record_bases -= line_bases
                remaining -= line_bases
            }
        }
    }
' > "$locate_fasta"

verify_stats_input() {
    local input="$1"
    local summary="$work_dir/summary.json"
    local stdin_summary="$work_dir/summary-stdin.json"
    "$binary" stats --format json "$input" > "$summary"
    "$binary" stats --format json < "$input" > "$stdin_summary"
    cmp "$summary" "$stdin_summary"

    local expected_records
    expected_records="$("$awk_bin" '/^>/ { count++ } END { print count + 0 }' "$input")"
    "$awk_bin" -v expected_records="$expected_records" -v expected_bases="$bases" '
        /"sequences":/ { records = $2; gsub(/,/, "", records) }
        /"total_len":/ { total = $2; gsub(/,/, "", total) }
        END { exit !(records == expected_records && total == expected_bases) }
    ' "$summary"
}

verify_stats_input "$stats_short_records"
verify_stats_input "$stats_long_records"

motif_for_length() {
    "$awk_bin" -v motif_size="$1" 'BEGIN { for (i = 1; i < motif_size; i++) printf "A"; print "C" }'
}

motif_8="$(motif_for_length 8)"
motif_32="$(motif_for_length 32)"
motif_64="$(motif_for_length 64)"
motif_128="$(motif_for_length 128)"
motif_512="$(motif_for_length 512)"

for motif in "$motif_8" "$motif_32" "$motif_64" "$motif_128" "$motif_512"; do
    output="$work_dir/no-hit.txt"
    "$binary" locate "$locate_fasta" "$motif" > "$output"
    if [[ -s "$output" ]]; then
        printf 'Expected no exact matches for motif %s in the homopolymer fixture\n' "$motif" >&2
        exit 1
    fi
done

frequent_output="$work_dir/frequent-hits.txt"
"$binary" locate "$locate_fasta" A > "$frequent_output"
frequent_count="$(wc -l < "$frequent_output" | tr -d '[:space:]')"
if [[ "$frequent_count" != "$bases" ]]; then
    printf 'Expected %s single-base matches, got %s\n' "$bases" "$frequent_count" >&2
    exit 1
fi

mismatch_output="$work_dir/mismatch-hits.txt"
"$binary" locate "$locate_fasta" "$motif_64" --max-mismatch 1 > "$mismatch_output"
expected_mismatch_hits="$("$awk_bin" -v bases="$bases" -v motif_length=64 '
    BEGIN {
        record_size = 10000
        full_records = int(bases / record_size)
        remainder = bases % record_size
        hits = full_records * (record_size - motif_length + 1)
        if (remainder >= motif_length) hits += remainder - motif_length + 1
        print hits
    }
')"
mismatch_count="$(wc -l < "$mismatch_output" | tr -d '[:space:]')"
if [[ "$mismatch_count" != "$expected_mismatch_hits" ]]; then
    printf 'Expected %s one-mismatch hits, got %s\n' "$expected_mismatch_hits" "$mismatch_count" >&2
    exit 1
fi

quote() {
    printf '%q' "$1"
}

q_binary="$(quote "$binary")"
q_short_records="$(quote "$stats_short_records")"
q_long_records="$(quote "$stats_long_records")"
q_locate_fasta="$(quote "$locate_fasta")"
q_motif_8="$(quote "$motif_8")"
q_motif_32="$(quote "$motif_32")"
q_motif_64="$(quote "$motif_64")"
q_motif_128="$(quote "$motif_128")"
q_motif_512="$(quote "$motif_512")"

printf 'Analysis benchmark input: %s bases, %s runs\n' "$bases" "$runs"
printf 'stats records: %s (about 100 bases each), %s (about 10,000 bases each)\n' \
    "$("$awk_bin" '/^>/ { count++ } END { print count + 0 }' "$stats_short_records")" \
    "$("$awk_bin" '/^>/ { count++ } END { print count + 0 }' "$stats_long_records")"
printf 'locate fixture: %s bases across 10,000-base records; correctness checks passed\n' "$bases"

"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'stats (about 100 bases per record)' "$q_binary stats --format json $q_short_records > /dev/null" \
    -n 'stats (about 10000 bases per record)' "$q_binary stats --format json $q_long_records > /dev/null" \
    -n 'locate (8-base near-match, no hits)' "$q_binary locate $q_locate_fasta $q_motif_8 > /dev/null" \
    -n 'locate (32-base near-match, no hits)' "$q_binary locate $q_locate_fasta $q_motif_32 > /dev/null" \
    -n 'locate (64-base near-match, no hits)' "$q_binary locate $q_locate_fasta $q_motif_64 > /dev/null" \
    -n 'locate (128-base near-match, no hits)' "$q_binary locate $q_locate_fasta $q_motif_128 > /dev/null" \
    -n 'locate (512-base near-match, no hits)' "$q_binary locate $q_locate_fasta $q_motif_512 > /dev/null" \
    -n 'locate (64-base, one mismatch allowed)' "$q_binary locate $q_locate_fasta $q_motif_64 --max-mismatch 1 > /dev/null" \
    -n 'locate (frequent one-base hits)' "$q_binary locate $q_locate_fasta A > /dev/null"
