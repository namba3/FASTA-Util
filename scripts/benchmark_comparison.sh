#!/usr/bin/env bash
set -euo pipefail

usage() {
    printf 'Usage: %s [BASES] [RUNS]\n' "${0##*/}"
    printf 'Defaults: BASES=20000000 RUNS=5 (BASES minimum: 10011000)\n'
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
seqkit="$(command -v seqkit || true)"
awk_bin="$(command -v awk || true)"
for tool in "$rustup_bin" "$hyperfine" "$seqkit" "$awk_bin"; do
    if [[ -z "$tool" || ! -x "$tool" ]]; then
        printf 'Required command not found: %s\n' "${tool:-rustup/hyperfine/seqkit/awk}" >&2
        exit 1
    fi
done

"$rustup_bin" run stable cargo build --release --locked --manifest-path "$repo_root/Cargo.toml" --bin fasta-util

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
binary="$repo_root/target/release/fasta-util"
input="$work_dir/input.fa"
fasta_index_input="$work_dir/fasta-util.fa"
seqkit_index_input="$work_dir/seqkit.fa"

# Alternate record lengths (50 kb and 25 kb) so the length filter selects
# records rather than acting as a pass-through.
"$rustup_bin" run stable cargo run --quiet --locked --manifest-path "$repo_root/generate_random_data/Cargo.toml" -- "$bases" --seed 42 |
    tr U T |
    "$awk_bin" '
        /^>/ { next }
        line_count == 0 { printf ">record_%06d description\n", record_count + 1 }
        { gsub(/[RYSWKMBDHVN]/, "A") }
        { print; line_count++ }
        line_count == (record_count % 2 == 0 ? 1000 : 500) { record_count++; line_count = 0 }
        END { if (line_count > 0) record_count++ }
    ' > "$input"

cp "$input" "$fasta_index_input"
cp "$input" "$seqkit_index_input"
"$binary" index "$fasta_index_input" > /dev/null
"$seqkit" --quiet faidx --update-faidx "$seqkit_index_input" > /dev/null

printf 'Comparative CLI benchmark: %s bases, seed 42, %s measured runs per case\n' "$bases" "$runs"
printf 'fasta-util %s; ' "$($binary --version)"
"$seqkit" version
"$rustup_bin" run stable rustc --version
"$hyperfine" --version

# Compare matching results before timing. The tools format several outputs
# differently, so normalize FASTA widths or locate columns where needed.
fasta_len="$("$binary" len -i "$input")"
seqkit_len="$("$seqkit" --quiet stats --tabular "$input" | "$awk_bin" -F '\t' 'NR == 2 { gsub(/,/, "", $5); print $5 }')"
if [[ "$fasta_len" != "$seqkit_len" ]]; then
    printf 'Length mismatch: fasta-util=%s seqkit=%s\n' "$fasta_len" "$seqkit_len" >&2
    exit 1
fi

"$binary" filter --min-len 40000 "$input" | "$binary" format --width 0 > "$work_dir/filter-fasta-util.fa"
"$seqkit" --quiet seq --min-len 40000 "$input" 2> /dev/null |
    "$seqkit" --quiet seq --line-width 0 2> /dev/null > "$work_dir/filter-seqkit.fa"
cmp "$work_dir/filter-fasta-util.fa" "$work_dir/filter-seqkit.fa"

"$binary" revcomp "$input" 2> /dev/null | "$binary" format --width 0 > "$work_dir/revcomp-fasta-util.fa"
"$seqkit" --quiet seq --seq-type dna --reverse --complement "$input" 2> /dev/null |
    "$seqkit" --quiet seq --line-width 0 2> /dev/null > "$work_dir/revcomp-seqkit.fa"
cmp "$work_dir/revcomp-fasta-util.fa" "$work_dir/revcomp-seqkit.fa"

"$binary" grep "$input" record_000201 | "$binary" format --width 0 > "$work_dir/grep-fasta-util.fa"
"$seqkit" --quiet grep --by-name --use-regexp --pattern record_000201 "$input" 2> /dev/null |
    "$seqkit" --quiet seq --line-width 0 > "$work_dir/grep-seqkit.fa"
cmp "$work_dir/grep-fasta-util.fa" "$work_dir/grep-seqkit.fa"

"$binary" format --width 80 "$input" > "$work_dir/format-fasta-util.fa"
"$seqkit" --quiet seq --line-width 80 "$input" 2> /dev/null > "$work_dir/format-seqkit.fa"
cmp "$work_dir/format-fasta-util.fa" "$work_dir/format-seqkit.fa"

"$binary" get "$fasta_index_input" record_000201:1001-11000 > "$work_dir/get-fasta-util.fa"
"$seqkit" --quiet faidx "$seqkit_index_input" record_000201:1001-11000 > "$work_dir/get-seqkit.fa"
cmp "$work_dir/get-fasta-util.fa" "$work_dir/get-seqkit.fa"

"$binary" locate "$input" ACGA > "$work_dir/locate-fasta-util.tsv"
"$seqkit" --quiet locate --pattern ACGA "$input" > "$work_dir/locate-seqkit.tsv"
LC_ALL=C "$awk_bin" -F '\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' "$work_dir/locate-fasta-util.tsv" |
    LC_ALL=C sort > "$work_dir/locate-fasta-util.sorted"
LC_ALL=C "$awk_bin" -F '\t' 'NR > 1 { print $1 "\t" $5 "\t" $6 "\t" $4 }' "$work_dir/locate-seqkit.tsv" |
    LC_ALL=C sort > "$work_dir/locate-seqkit.sorted"
cmp "$work_dir/locate-fasta-util.sorted" "$work_dir/locate-seqkit.sorted"

q() { printf '%q' "$1"; }
q_binary="$(q "$binary")"
q_seqkit="$(q "$seqkit")"
q_input="$(q "$input")"
q_fasta_index_input="$(q "$fasta_index_input")"
q_seqkit_index_input="$(q "$seqkit_index_input")"

printf '\n%s\n' 'Comparable operations (one warmup, output suppressed):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    -n 'fasta-util len' "$q_binary len -i $q_input > /dev/null" \
    -n 'seqkit stats (sum_len)' "$q_seqkit --quiet stats --tabular $q_input > /dev/null 2>&1" \
    -n 'fasta-util stats' "$q_binary stats $q_input > /dev/null" \
    -n 'seqkit stats --all' "$q_seqkit --quiet stats --all $q_input > /dev/null 2>&1" \
    -n 'fasta-util filter (min-len 40K)' "$q_binary filter --min-len 40000 $q_input > /dev/null" \
    -n 'seqkit seq (min-len 40K)' "$q_seqkit --quiet seq --min-len 40000 $q_input > /dev/null 2>&1" \
    -n 'fasta-util revcomp' "$q_binary revcomp $q_input > /dev/null" \
    -n 'seqkit seq (reverse complement)' "$q_seqkit --quiet seq --seq-type dna --reverse --complement $q_input > /dev/null 2>&1" \
    -n 'fasta-util grep (header)' "$q_binary grep $q_input record_000201 > /dev/null" \
    -n 'seqkit grep (header)' "$q_seqkit --quiet grep --by-name --use-regexp --pattern record_000201 $q_input > /dev/null 2>&1" \
    -n 'fasta-util locate (ACGA)' "$q_binary locate $q_input ACGA > /dev/null" \
    -n 'seqkit locate (ACGA)' "$q_seqkit --quiet locate --pattern ACGA $q_input > /dev/null 2>&1" \
    -n 'fasta-util format (width 80)' "$q_binary format --width 80 $q_input > /dev/null" \
    -n 'seqkit seq (width 80)' "$q_seqkit --quiet seq --line-width 80 $q_input > /dev/null 2>&1" \
    -n 'fasta-util get (FAI, 10K)' "$q_binary get $q_fasta_index_input record_000201:1001-11000 > /dev/null" \
    -n 'seqkit faidx (10K)' "$q_seqkit --quiet faidx $q_seqkit_index_input record_000201:1001-11000 > /dev/null 2>&1"

printf '\n%s\n' 'Index creation (index files removed before every run):'
"$hyperfine" --shell=bash --warmup 1 --runs "$runs" --style basic \
    --prepare "rm -f $q_fasta_index_input.fai $q_seqkit_index_input.fai" \
    -n 'fasta-util index' "$q_binary index $q_fasta_index_input > /dev/null 2>&1" \
    -n 'seqkit faidx index' "$q_seqkit --quiet faidx --update-faidx $q_seqkit_index_input > /dev/null 2>&1"
