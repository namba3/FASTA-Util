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

When `-o`/`--output` specifies a file, the destination is replaced only after processing succeeds. The input file cannot also be used as the output file.

```sh
./target/release/fasta-util slice -i test.fna --range 99..=199
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

The benchmark used the source FASTA directly and preserved its lowercase soft-masked sequence. Before timing each slice case, the `seqret` and `fasta-util` outputs were verified byte-for-byte.

Environment: Ubuntu 26.04.1 LTS (WSL2), AMD Ryzen 9 9900X, stable Rust 1.99.0, seqkit 2.10.1, EMBOSS seqret 6.6.0.0, and hyperfine 1.20.0. Each command had one warmup followed by five measured runs. Tables show the mean and standard deviation. Results vary with the machine and file-cache state.

### len

`seqkit stats` computes statistics for all 705 records; this tool counts sequence symbols. Their total lengths matched.

| Command | Result | Time (mean ± standard deviation) |
| --- | ---: | ---: |
| `seqkit stats` | 3,298,430,636 bases, 705 records | 1,618 ± 97 ms |
| `fasta-util len` | 3,298,430,636 bases | 2,646 ± 284 ms |

### slice

Ranges are positions within chromosome 1. `seqret` uses one-based inclusive coordinates; this tool uses zero-based inclusive ranges. Output was wrapped at 60 bases per line.

| Offset | Slice length | seqret (mean ± standard deviation) | fasta-util (mean ± standard deviation) |
| ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 843 ± 46 ms | 643 ± 68 ms |
| 100,000,000 | 100,000 | 745 ± 39 ms | 294 ± 36 ms |
| 100,000,000 | 100 | 731 ± 113 ms | 416 ± 124 ms |
| 0 | 100,000,000 | 866 ± 43 ms | 344 ± 69 ms |
| 0 | 100,000 | 769 ± 113 ms | 2.3 ± 0.3 ms |
| 0 | 100 | 618 ± 137 ms | 1.6 ± 0.2 ms |

To reproduce these measurements, install `seqkit`, `seqret`, `hyperfine`, and stable Rust, then run this script from the repository root. An alternate FASTA path can be supplied as the first argument.

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```
