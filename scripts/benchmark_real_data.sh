#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
dataset="${1:-$repo_root/dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna}"
binary="$repo_root/target/release/fasta-util"
rustup_bin="$(command -v rustup || true)"
seqkit="$(command -v seqkit || true)"
seqret="$(command -v seqret || true)"
if [[ -z "$seqret" && -x "$repo_root/dataset/emboss/bin/seqret" ]]; then
    seqret="$repo_root/dataset/emboss/bin/seqret"
fi
hyperfine="$(command -v hyperfine || true)"
awk_bin="$(command -v awk || true)"

for tool in "$rustup_bin" "$seqkit" "$hyperfine" "$awk_bin"; do
    if [[ -z "$tool" || ! -x "$tool" ]]; then
        printf 'Required command not found: %s\n' "${tool:-rustup/seqkit/hyperfine/awk}" >&2
        exit 1
    fi
done
if [[ -n "$seqret" && ! -x "$seqret" ]]; then
    seqret=""
fi
if [[ ! -f "$dataset" ]]; then
    printf 'FASTA dataset not found: %s\n' "$dataset" >&2
    exit 1
fi

"$rustup_bin" run stable cargo build --release --manifest-path "$repo_root/Cargo.toml" --bin fasta-util

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
chr1="$work_dir/chr1.fna"
fai="$chr1.fai"
# Keep lowercase soft-masking from the original FASTA in the benchmark input.
"$seqkit" grep -p NC_000001.11 "$dataset" > "$chr1"
LC_ALL=C "$awk_bin" -f "$repo_root/scripts/write_fai.awk" "$chr1" > "$fai"

printf '%s\n' 'Tool versions:'
"$seqkit" version
if [[ -n "$seqret" ]]; then
    "$seqret" -version
else
    printf '%s\n' 'EMBOSS seqret not found; external comparison will be skipped'
fi
"$rustup_bin" run stable rustc --version
printf '\n%s\n' 'Dataset statistics:'
printf 'Source: %s\n' "$dataset"
printf 'Benchmark input: original FASTA, preserving lowercase soft-masking\n'
printf '%s: %s bytes\n' "$dataset" "$(wc -c < "$dataset" | awk '{print $1}')"
printf '%s: %s bytes\n' "$chr1" "$(wc -c < "$chr1" | awk '{print $1}')"
printf '%s: %s bytes\n' "$fai" "$(wc -c < "$fai" | awk '{print $1}')"
"$seqkit" stats "$dataset" "$chr1"
printf 'fasta-util len (assembly): '
"$binary" len -i "$dataset"
printf 'fasta-util len (chr1): '
"$binary" len -i "$chr1"

q() {
    printf '%q' "$1"
}

q_dataset="$(q "$dataset")"
q_chr1="$(q "$chr1")"
q_binary="$(q "$binary")"
q_seqkit="$(q "$seqkit")"
q_fai="$(q "$fai")"

printf '\n%s\n' 'Full assembly length benchmark (1 warmup, 5 measured runs):'
"$hyperfine" --shell=none --warmup 1 --runs 5 --style basic \
    -n 'seqkit stats' "$q_seqkit stats $q_dataset" \
    -n 'fasta-util len' "$q_binary len -i $q_dataset"

printf '\n%s\n' 'chr1 get benchmarks (1 warmup, 5 measured runs):'
while read -r name start length; do
    range_start=$((start + 1))
    range_end=$((start + length))
    fasta_command="$q_binary get $q_chr1 ${range_start}-${range_end} --chars-per-line=60"
    indexed_command="$q_binary get $q_chr1 ${range_start}-${range_end} --fai-index $q_fai --chars-per-line=60"

    "$binary" get "$chr1" "${range_start}-${range_end}" --chars-per-line=60 > "$work_dir/fasta-util.out"
    "$binary" get "$chr1" "${range_start}-${range_end}" --fai-index "$fai" --chars-per-line=60 > "$work_dir/fasta-util-indexed.out"
    if ! cmp -s "$work_dir/fasta-util.out" "$work_dir/fasta-util-indexed.out"; then
        printf 'Output mismatch for get case %s\n' "$name" >&2
        exit 1
    fi
    if [[ -n "$seqret" ]]; then
        seqret_start=$((start + 1))
        seqret_end=$((start + length))
        "$seqret" -sequence "$chr1" -sbegin "$seqret_start" -send "$seqret_end" -auto -stdout > "$work_dir/seqret.out"
        if ! cmp -s "$work_dir/seqret.out" "$work_dir/fasta-util.out"; then
            printf 'seqret output mismatch for get case %s\n' "$name" >&2
            exit 1
        fi
    fi

    printf '\n%s\n' "$name: zero-based offset=$start length=$length (get uses 1-based coordinates)"
    hyperfine_args=(--shell=none --warmup 1 --runs 5 --style basic)
    if [[ -n "$seqret" ]]; then
        seqret_command="$(q "$seqret") -sequence $q_chr1 -sbegin $((start + 1)) -send $((start + length)) -auto -stdout"
        hyperfine_args+=(-n seqret "$seqret_command")
    fi
    hyperfine_args+=(-n fasta-util "$fasta_command" -n 'fasta-util (FAI)' "$indexed_command")
    "$hyperfine" "${hyperfine_args[@]}"
done <<'CASES'
middle_100M 100000000 100000000
middle_100K 100000000 100000
middle_100 100000000 100
start_100M 0 100000000
start_100K 0 100000
start_100 0 100
CASES
