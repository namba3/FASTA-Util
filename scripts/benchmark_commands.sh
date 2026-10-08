#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s [BASES] [RUNS]\n' "${0##*/}"
    printf 'Defaults: BASES=20000000 RUNS=5\n'
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
if [[ $# -gt 2 ]]; then
    usage >&2
    exit 2
fi

bases="${1:-20000000}"
runs="${2:-5}"
if [[ ! "$bases" =~ ^[1-9][0-9]*$ || ! "$runs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'BASES and RUNS must be positive decimal integers\n' >&2
    usage >&2
    exit 2
fi
if (( bases < 10011000 )); then
    printf 'BASES must be at least 10011000 to benchmark an indexed region in record_000201\n' >&2
    exit 2
fi

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
rustup_bin="$(command -v rustup || true)"
hyperfine="$(command -v hyperfine || true)"
awk_bin="$(command -v awk || true)"
for tool in "$rustup_bin" "$hyperfine" "$awk_bin"; do
    if [[ -z "$tool" || ! -x "$tool" ]]; then
        printf 'Required command not found: %s\n' "${tool:-rustup/hyperfine/awk}" >&2
        exit 1
    fi
done

"$rustup_bin" run stable cargo build --release --locked --manifest-path "$repo_root/Cargo.toml" --bin fasta-util

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
binary="$repo_root/target/release/fasta-util"
nucleotide="$work_dir/nucleotide.fa"
protein="$work_dir/protein.fa"
scan_input="$work_dir/scan.fa"
indexed_input="$work_dir/indexed.fa"
index_input="$work_dir/to-index.fa"

"$rustup_bin" run stable cargo run --quiet --locked --manifest-path "$repo_root/generate_random_data/Cargo.toml" -- "$bases" --seed 42 |
    tr U T |
    "$awk_bin" '
        /^>/ { next }
        line_count % 1000 == 0 { printf ">record_%06d description\n", record_count + 1; record_count++ }
        { print; line_count++ }
    ' > "$nucleotide"

# Use a deterministic protein alphabet with the same total number of symbols.
"$awk_bin" -v bases="$bases" 'BEGIN {
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    printf ">protein_000001 benchmark protein\n"
    for (i = 0; i < bases; i++) {
        printf "%s", substr(alphabet, (i % length(alphabet)) + 1, 1)
        if ((i + 1) % 60 == 0 || i + 1 == bases) printf "\n"
    }
}' > "$protein"

cp "$nucleotide" "$scan_input"
cp "$nucleotide" "$indexed_input"
cp "$nucleotide" "$index_input"
"$binary" index "$indexed_input"

printf 'Stable Rust command benchmark: %s bases, seed 42, %s measured runs per case\n' "$bases" "$runs"
printf 'Input sizes: nucleotide=%s bytes protein=%s bytes; records=%s\n' \
    "$(wc -c < "$nucleotide" | tr -d ' ')" "$(wc -c < "$protein" | tr -d ' ')" \
    "$(awk '/^>/ { count++ } END { print count }' "$nucleotide")"
"$rustup_bin" run stable rustc --version
"$hyperfine" --version

# Verify representative results before timing. File/stdin parity is checked for
# commands that support both forms; get is compared with and without the index.
for subcommand in stats composition filter revcomp format validate; do
    case "$subcommand" in
        filter) args=(--min-len 1) ;;
        format) args=(--width 60) ;;
        *) args=() ;;
    esac
    "$binary" "$subcommand" "${args[@]}" "$nucleotide" > "$work_dir/file.out"
    "$binary" "$subcommand" "${args[@]}" < "$nucleotide" > "$work_dir/stdin.out"
    cmp "$work_dir/file.out" "$work_dir/stdin.out"
done
"$binary" get "$scan_input" record_000201:1001-11000 > "$work_dir/get-scan.out"
"$binary" get "$indexed_input" record_000201:1001-11000 > "$work_dir/get-indexed.out"
cmp "$work_dir/get-scan.out" "$work_dir/get-indexed.out"
"$binary" stats "$protein" --sequence-type protein > "$work_dir/protein-stats.out"
"$binary" composition "$protein" --sequence-type protein > "$work_dir/protein-composition.out"
"$binary" validate "$protein" --sequence-type protein > "$work_dir/protein-validate.out"
"$binary" grep "$nucleotide" record_000001 > "$work_dir/grep.out"
"$binary" locate "$nucleotide" ACGT --max-mismatch 1 > "$work_dir/locate.out"
"$binary" len -i "$nucleotide" > "$work_dir/len.out"

run_analysis() {
    local command="$1"
    local threads="$2"
    local output="$3"
    local thread_args=()
    if [[ "$threads" != "auto" ]]; then
        thread_args=(--threads "$threads")
    fi
    case "$command" in
        len) "$binary" len -i "$nucleotide" "${thread_args[@]}" > "$output" ;;
        stats) "$binary" stats "$nucleotide" "${thread_args[@]}" > "$output" ;;
        composition) "$binary" composition "$nucleotide" "${thread_args[@]}" > "$output" ;;
    esac
}

# Check output parity before comparing worker counts.
for subcommand in len stats composition; do
    run_analysis "$subcommand" 1 "$work_dir/$subcommand-serial.out"
    for threads in auto 2 4; do
        run_analysis "$subcommand" "$threads" "$work_dir/$subcommand-$threads.out"
        cmp "$work_dir/$subcommand-serial.out" "$work_dir/$subcommand-$threads.out"
    done
done

q() { printf '%q' "$1"; }
q_binary="$(q "$binary")"
q_nucleotide="$(q "$nucleotide")"
q_protein="$(q "$protein")"
q_scan="$(q "$scan_input")"
q_indexed="$(q "$indexed_input")"
q_index_input="$(q "$index_input")"
q_cat="$(q "$(command -v cat)")"

printf '\n%s\n' 'All subcommands (1 warmup, output suppressed):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'len (nucleotide)' "$q_binary len -i $q_nucleotide > /dev/null" \
    -n 'validate (nucleotide)' "$q_binary validate $q_nucleotide > /dev/null" \
    -n 'validate (protein)' "$q_binary validate $q_protein --sequence-type protein > /dev/null" \
    -n 'index (write FAI)' "$q_binary index $q_index_input" \
    -n 'stats (nucleotide)' "$q_binary stats $q_nucleotide > /dev/null" \
    -n 'stats (protein)' "$q_binary stats $q_protein --sequence-type protein > /dev/null" \
    -n 'composition (nucleotide)' "$q_binary composition $q_nucleotide > /dev/null" \
    -n 'composition (protein)' "$q_binary composition $q_protein --sequence-type protein > /dev/null" \
    -n 'get (scan 10K)' "$q_binary get $q_scan record_000201:1001-11000 > /dev/null" \
    -n 'get (FAI 10K)' "$q_binary get $q_indexed record_000201:1001-11000 > /dev/null" \
    -n 'filter (nucleotide)' "$q_binary filter --min-len 1 $q_nucleotide > /dev/null" \
    -n 'filter (protein)' "$q_binary filter --min-len 1 --sequence-type protein $q_protein > /dev/null" \
    -n 'revcomp' "$q_binary revcomp $q_nucleotide > /dev/null" \
    -n 'grep' "$q_binary grep $q_nucleotide record_000001 > /dev/null" \
    -n 'locate (mismatch=1)' "$q_binary locate $q_nucleotide ACGT --max-mismatch 1 > /dev/null" \
    -n 'format (width=60)' "$q_binary format --width 60 $q_nucleotide > /dev/null"

printf '\n%s\n' 'Standard input and pipeline (1 warmup, output suppressed):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'stats (stdin)' "$q_cat $q_nucleotide | $q_binary stats > /dev/null" \
    -n 'composition (stdin)' "$q_cat $q_nucleotide | $q_binary composition > /dev/null" \
    -n 'filter (stdin)' "$q_cat $q_nucleotide | $q_binary filter --min-len 1 > /dev/null" \
    -n 'revcomp (stdin)' "$q_cat $q_nucleotide | $q_binary revcomp > /dev/null" \
    -n 'format (stdin)' "$q_cat $q_nucleotide | $q_binary format --width 60 > /dev/null" \
    -n 'validate (stdin)' "$q_cat $q_nucleotide | $q_binary validate > /dev/null" \
    -n 'filter | revcomp | stats' "$q_binary filter --min-len 1 $q_nucleotide | $q_binary revcomp | $q_binary stats > /dev/null"

printf '\n%s\n' 'Analysis worker-count comparison (1 warmup, output suppressed):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'len (1 worker)' "$q_binary len -i $q_nucleotide --threads 1 > /dev/null" \
    -n 'len (2 workers)' "$q_binary len -i $q_nucleotide --threads 2 > /dev/null" \
    -n 'len (4 workers)' "$q_binary len -i $q_nucleotide --threads 4 > /dev/null" \
    -n 'len (auto)' "$q_binary len -i $q_nucleotide > /dev/null" \
    -n 'stats (1 worker)' "$q_binary stats $q_nucleotide --threads 1 > /dev/null" \
    -n 'stats (2 workers)' "$q_binary stats $q_nucleotide --threads 2 > /dev/null" \
    -n 'stats (4 workers)' "$q_binary stats $q_nucleotide --threads 4 > /dev/null" \
    -n 'stats (auto)' "$q_binary stats $q_nucleotide > /dev/null" \
    -n 'composition (1 worker)' "$q_binary composition $q_nucleotide --threads 1 > /dev/null" \
    -n 'composition (2 workers)' "$q_binary composition $q_nucleotide --threads 2 > /dev/null" \
    -n 'composition (4 workers)' "$q_binary composition $q_nucleotide --threads 4 > /dev/null" \
    -n 'composition (auto)' "$q_binary composition $q_nucleotide > /dev/null"
