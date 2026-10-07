#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s [BASES] [RUNS]\n' "${0##*/}"
    printf 'Defaults: BASES=1000000 RUNS=5\n'
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
if [[ ! "$bases" =~ ^[1-9][0-9]*$ || ! "$runs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'BASES and RUNS must be positive decimal integers\n' >&2
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
input="$work_dir/input.fa"
trap 'rm -rf "$work_dir"' EXIT
binary="$repo_root/target/release/fasta-util"
"$rustup_bin" run stable cargo run --quiet --locked --manifest-path "$repo_root/generate_random_data/Cargo.toml" -- "$bases" --seed 42 |
    tr U T |
    "$awk_bin" '
        /^>/ { next }
        line_count % 1000 == 0 { printf ">record_%d\n", record_count + 1; record_count++ }
        { print; line_count++ }
    ' > "$input"

# Keep the seeded input outside the temporary workspace so its size is easy to inspect.
input_size="$(wc -c < "$input" | tr -d ' ')"
printf 'Stable Rust benchmark input: %s bases (%s bytes), seed 42, %s runs\n' \
    "$bases" "$input_size" "$runs"
"$rustup_bin" run stable rustc --version
"$hyperfine" --version

# Check that file and stdin modes produce identical output before timing them.
"$binary" stats "$input" > "$work_dir/stats-file.txt"
"$binary" stats < "$input" > "$work_dir/stats-stdin.txt"
cmp "$work_dir/stats-file.txt" "$work_dir/stats-stdin.txt"

"$binary" filter --min-len 1 "$input" > "$work_dir/filter-file.fa"
"$binary" filter --min-len 1 < "$input" > "$work_dir/filter-stdin.fa"
cmp "$work_dir/filter-file.fa" "$work_dir/filter-stdin.fa"

"$binary" revcomp "$input" > "$work_dir/revcomp-file.fa"
"$binary" revcomp < "$input" > "$work_dir/revcomp-stdin.fa"
cmp "$work_dir/revcomp-file.fa" "$work_dir/revcomp-stdin.fa"

"$binary" filter --min-len 1 "$input" |
    "$binary" revcomp |
    "$binary" stats > "$work_dir/pipeline-stats.txt"
cat "$input" |
    "$binary" filter --min-len 1 |
    "$binary" revcomp |
    "$binary" stats > "$work_dir/pipeline-stdin-stats.txt"
"$binary" stats "$work_dir/revcomp-file.fa" > "$work_dir/expected-pipeline-stats.txt"
cmp "$work_dir/expected-pipeline-stats.txt" "$work_dir/pipeline-stats.txt"
cmp "$work_dir/pipeline-stats.txt" "$work_dir/pipeline-stdin-stats.txt"

quote() {
    printf '%q' "$1"
}

q_binary="$(quote "$binary")"
q_input="$(quote "$input")"
printf '\n%s\n' 'File input and stdin input (1 warmup; output suppressed where applicable):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'stats (file)' "$q_binary stats $q_input > /dev/null" \
    -n 'stats (stdin)' "cat $q_input | $q_binary stats > /dev/null" \
    -n 'filter (file)' "$q_binary filter --min-len 1 $q_input > /dev/null" \
    -n 'filter (stdin)' "cat $q_input | $q_binary filter --min-len 1 > /dev/null" \
    -n 'revcomp (file)' "$q_binary revcomp $q_input > /dev/null" \
    -n 'revcomp (stdin)' "cat $q_input | $q_binary revcomp > /dev/null"

printf '\n%s\n' 'Three-command pipeline (filter | revcomp | stats):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'pipeline (filter reads file)' \
    "$q_binary filter --min-len 1 $q_input | $q_binary revcomp | $q_binary stats > /dev/null" \
    -n 'pipeline (filter reads stdin)' \
    "cat $q_input | $q_binary filter --min-len 1 | $q_binary revcomp | $q_binary stats > /dev/null"
