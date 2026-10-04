# FASTA-Util

[日本語](README.md) | [English](README.en.md)

A CLI tool for playing with FASTA files

## About the FASTA format

FASTA is a text format for nucleotide or amino acid sequences. Each record starts with a header line beginning with `>`, followed by one or more sequence lines. Sequence data can wrap across multiple lines.

This tool currently processes nucleotide FASTA files. Amino acid sequences are not supported.

```fasta
>record-1 optional description
ACGTN
UKS-
>record-2
MRY
```

`len` counts sequence symbols, excluding headers and blank lines. `slice` positions also count only sequence symbols, excluding headers and line breaks. For multiple records, the range indexes the sequences concatenated in file order. Headers encountered before the range ends are preserved, so the output may include a header for a record that contributes no symbols to the selected range. This tool accepts uppercase and lowercase `ACGTNUKSYMWRBDHV` symbols and `-` for a gap. `slice` preserves the original letter case. Other symbols are rejected.

## Build

```sh
cargo build --release
```

## Generate Test Data

```sh
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 > test.fna
```

Pass `--seed` to reproduce generated data. Identical output is guaranteed when using the same `rand` version and runtime environment.
Use `--line-width` to change the number of bases per sequence line. The default is 50.

```sh
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 --seed 42 > test.fna
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 --line-width 60 > test.fna
```

## Sub Commands

### len

Count the total length of the sequence

```sh
./target/release/fasta-util len -i test.fna
```

```txt
10000
```

### slice

Cut out a part of the sequence

`--range` uses zero-based Rust range syntax. `2..10` selects positions 2 through 9, while `2..=10` includes position 10. `..10` selects from the beginning through position 9, and `2..` selects from position 2 to the end. The default `..` selects the entire sequence. `--chars-per-line` controls output wrapping and defaults to 60.

When `-o`/`--output` specifies a file, the destination is replaced only after processing succeeds. The input file cannot also be used as the output file.

For a slice from the middle of a large uncompressed FASTA, pass a matching `.fai` index with `--fai-index` to read the selected region directly. Create the index with `samtools faidx`. Rebuild it whenever the FASTA changes. This path validates the sequence symbols in the selected region.

```sh
./target/release/fasta-util slice -i test.fna --range 99..=199
samtools faidx test.fna
./target/release/fasta-util slice -i test.fna --fai-index test.fna.fai --range 100000000..100000100
```

```txt
>TestData 10000 random data
WDCAGVUTRABAKRRNRNHHKTYDNBNTCHMRBRRYHWHKYBHKSBAHVNTCGUMGCMMA
GYMDSVCYRAMWNURRVTCYYCYCWWHTRCAUVSBUVHMHNWTGKGHGATWMHYTWNSUB
SUDKUGDWWTSSYBUCKYUDSAADMMRHMT
```

## Simple tests and benchmarks

Run the nucleotide-check benchmark with:
It measures four input patterns: all valid symbols, all invalid symbols, a 50% valid mix, and a 99% valid mix.
Each implementation is measured for five 200 ms samples, and the median is reported.
The implementation order rotates between samples to reduce order bias.

```sh
cargo bench --bench nucleic_acid
```

Set the input size and duration of each sample with these options. The defaults are 10,000 bytes and 200 ms.

```sh
cargo bench --bench nucleic_acid -- --input-size 100000 --sample-ms 500
```

## Real-data benchmarks

Measurements were taken on 2026-10-04 using the RefSeq GRCh38.p14 FASTA in `dataset/ncbi_dataset`. The full-genome length comparison used `GCF_000001405.40_GRCh38.p14_genomic.fna` (3,339,739,109 bytes), containing 705 records. Slice comparisons used its chromosome 1 record (`NC_000001.11`, 248,956,422 bases).

The benchmark used the source FASTA directly and preserved its lowercase soft-masked sequence. Before timing each slice case, regular `slice` and `.fai` indexed `slice` outputs were verified byte-for-byte. If `seqret` is installed, its output is also verified and included in the comparison.

The regular and `.fai` paths were remeasured by the benchmark script on 2026-10-04 using the same CPU, OS, stable Rust, and hyperfine setup. Their outputs were compared byte-for-byte before five timed runs. `seqret` was unavailable during this remeasurement, so only its values in the table are from the earlier run. Hyperfine reported a statistical outlier for the 100-base slice from the start; treat the small timing difference from the regular path as measurement noise.

Environment: Ubuntu 26.04.1 LTS (WSL2), AMD Ryzen 9 9900X, stable Rust 1.99.0, seqkit 2.10.1, EMBOSS seqret 6.6.0.0, and hyperfine 1.20.0. Each command had one warmup followed by five measured runs. Tables show the mean and standard deviation. Results vary with the machine and file-cache state.

### len

`seqkit stats` computes statistics for all 705 records; this tool counts sequence symbols. Their total lengths matched.

| Command | Result | Time (mean ± standard deviation) |
| --- | ---: | ---: |
| `seqkit stats` | 3,298,430,636 bases, 705 records | 1,618 ± 97 ms |
| `fasta-util len` | 3,298,430,636 bases | 2,646 ± 284 ms |

### slice

Ranges are positions within chromosome 1. `seqret` uses one-based inclusive coordinates; this tool uses zero-based inclusive ranges. Output was wrapped at 60 bases per line.

| Offset | Slice length | seqret (mean ± standard deviation) | fasta-util (regular, mean ± standard deviation) | fasta-util (FAI, mean ± standard deviation) |
| ---: | ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 843 ± 46 ms | 692.4 ± 51.3 ms | 69.5 ± 9.4 ms |
| 100,000,000 | 100,000 | 745 ± 39 ms | 548.1 ± 293.6 ms | 1.4 ± 0.1 ms |
| 100,000,000 | 100 | 731 ± 113 ms | 434.5 ± 60.4 ms | 1.6 ± 0.2 ms |
| 0 | 100,000,000 | 866 ± 43 ms | 384.1 ± 40.5 ms | 62.8 ± 2.1 ms |
| 0 | 100,000 | 769 ± 113 ms | 2.8 ± 0.4 ms | 1.7 ± 0.2 ms |
| 0 | 100 | 618 ± 137 ms | 2.1 ± 0.2 ms | 1.9 ± 0.3 ms |

To reproduce these measurements, install `awk`, `seqkit`, `hyperfine`, and stable Rust, then run this script from the repository root. `seqret` is optional and adds an external-tool comparison when installed. An alternate FASTA path can be supplied as the first argument.

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```
