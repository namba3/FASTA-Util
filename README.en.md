# FASTA-Util

[日本語](README.md) | [English](README.en.md)

A CLI tool for playing with FASTA files

## About the FASTA format

FASTA is a text format for nucleotide or amino acid sequences. Each record starts with a header line beginning with `>`, followed by one or more sequence lines. Sequence data can wrap across multiple lines.

`len` and `get` process nucleotide sequences by default. To process protein sequences, pass `--sequence-type protein`. `revcomp` and `locate` process nucleotide sequences, while `filter`, `stats`, and `composition` detect sequence type automatically.

```fasta
>record-1 optional description
ACGTN
UKS-
>record-2
MRY
```

`len` counts sequence symbols, excluding headers and blank lines. Numeric `get` ranges use 1-based inclusive coordinates across the sequences concatenated in file order. Headers encountered before the range ends are preserved, so the output may include a header for a record that contributes no symbols to the selected range.

Nucleotide sequences accept uppercase and lowercase `ACGTNUKSYMWRBDHV` symbols and `-` for a gap. Protein sequences accept uppercase and lowercase symbols for the 20 standard amino acids, plus `B J O U X Z`, the stop marker `*`, and the gap symbol `-`. `get` preserves the original letter case in both modes. Other symbols are rejected.

Do not modify an input file while running `len` or `get`. While running `get` with an index, do not modify either the FASTA input or its `.fai` index.

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

Commands that produce sequence data write to standard output unless an output file is specified. `stats`, `composition`, `filter`, `revcomp`, `format`, and `validate` read standard input when their input path is omitted. For `get`, `grep`, and `locate`, pass `-` in the input position. `len` reads standard input by default and uses `-i` for file input. `index` requires a file input so it can determine where to write the `.fai` sidecar.

```sh
fasta-util filter --min-len 1000 input.fa |
  fasta-util revcomp |
  fasta-util stats
```

### len

Count the total length of the sequence

```sh
./target/release/fasta-util len -i test.fna
./target/release/fasta-util len --sequence-type protein -i proteins.faa
```

```txt
10000
```

### stats

Show the record count, total/minimum/maximum/mean length, N50, GC percentage, `N` percentage, and sequence type. Add `--each` to show the ID, length, GC percentage, and `N` percentage for every record.

The default `auto` mode selects protein when it finds a protein symbol that is not valid in nucleotide sequences. Otherwise, it infers DNA or RNA from `T` and `U`; if neither appears, the type is reported as ambiguous. Specify `--sequence-type nucleotide` or `--sequence-type protein` to choose explicitly. GC and `N` percentages are `n/a` in text and `null` in JSON for protein sequences.

```sh
./target/release/fasta-util stats genome.fa
./target/release/fasta-util stats genome.fa --each
./target/release/fasta-util stats genome.fa --format json
./target/release/fasta-util stats proteins.faa --sequence-type protein --format json
```

### composition

Show the fraction of each symbol across all sequence data. For nucleotides, the output includes `A`, `C`, `G`, `T` (or `U` for RNA), `N`, any IUPAC ambiguity symbols or gaps found in the input, and `GC` for the combined `G` and `C` fraction. For proteins, it reports the 20 standard amino acids and any extended symbols found. The denominator includes every sequence symbol, including gaps and ambiguity symbols. Sequence type is detected automatically as in `stats`; specify `--sequence-type nucleotide` or `--sequence-type protein` when needed.

```sh
./target/release/fasta-util composition seq.fa
./target/release/fasta-util composition proteins.faa --sequence-type protein
cat seq.fa | ./target/release/fasta-util composition
```

### filter

Select records by sequence length, GC fraction, and `N` fraction. Fractions use values from `0.0` to `1.0`; records must satisfy every supplied condition. GC and `N` fractions use the full sequence length as the denominator, including gaps and ambiguous symbols. Sequence type is detected automatically; a protein-only symbol selects protein mode. For proteins that cannot be distinguished from nucleotide sequences, specify `--sequence-type protein`. GC and `N` conditions are available only for nucleotide sequences. Matching records are written in input order to standard output, or to a file with `-o`/`--output`.

```sh
./target/release/fasta-util filter seq.fa --min-len 1000
./target/release/fasta-util filter seq.fa --max-len 10000
./target/release/fasta-util filter seq.fa --min-gc 0.40 --max-gc 0.60
./target/release/fasta-util filter seq.fa --max-n 0.05
./target/release/fasta-util filter proteins.fa --min-len 100 --sequence-type protein
```

### revcomp

Reverse-complement each nucleotide sequence. For DNA, `A` complements to `T`; for RNA, `A` complements to `U`. IUPAC ambiguity symbols and the `-` gap are supported, and letter case is preserved. A record containing both `T` and `U` is rejected. Output wraps at 60 bases by default; set another width with `--chars-per-line`. Use `-o`/`--output` to write to a file.

```sh
./target/release/fasta-util revcomp seq.fa
./target/release/fasta-util revcomp seq.fa --chars-per-line 80
./target/release/fasta-util revcomp seq.fa --output seq.revcomp.fa
```

### grep

Search record headers for literal text and print each matching FASTA record in input order. Matching is case-sensitive by default. Use `--ignore-case` for ASCII case-insensitive matching and `--invert-match` to select records without a match.

```sh
./target/release/fasta-util grep seq.fa BRCA
./target/release/fasta-util grep seq.fa brca --ignore-case
```

### format

Normalize FASTA line endings to LF and wrap sequences to the requested width. The default width is 60. `--width 0` writes each record's sequence on one line. Use `--uppercase` or `--lowercase` to change sequence letter case, `--remove-gaps` to remove `-`, and `--trim-header` to trim surrounding ASCII whitespace after `>`. Output goes to standard output unless a file is specified.

```sh
./target/release/fasta-util format --width 80 seq.fa > formatted.fa
./target/release/fasta-util format --width 0 seq.fa
./target/release/fasta-util format --uppercase --remove-gaps --trim-header seq.fa
```

### locate

Search nucleotide sequences on both strands and print tab-separated `ID`, start, end, and strand columns. Coordinates are 1-based and inclusive. Motifs accept IUPAC nucleotide symbols, and overlapping matches are reported. Set `--max-mismatch` to allow substitutions. Matches in the input orientation use `+`; matches to the reverse complement use `-`. At ambiguous sequence positions, a match is reported when the possible-base sets intersect at every position.

```sh
./target/release/fasta-util locate genome.fa AATAAA
./target/release/fasta-util locate seq.fa ATGNNNTAA
./target/release/fasta-util locate seq.fa AATAAA --max-mismatch 1
```

### get

Get complete records by ID, record regions using `ID:START-END`, or a global range using `START-END`. Coordinates are 1-based and inclusive. Global ranges index all record sequences concatenated in FASTA order. Multiple IDs and an ID file passed with `--ids` are supported; output follows the records' order in the FASTA file. If a `.fai` sidecar exists next to the input, regions are read directly; otherwise, the file is scanned from the beginning. Pass `--fai-index` to use a specific index. Set output wrapping with `--chars-per-line`.

```sh
./target/release/fasta-util get genome.fa chr1
./target/release/fasta-util get genome.fa chr1 chr3 chrX
./target/release/fasta-util get genome.fa --ids chromosomes.txt
./target/release/fasta-util get genome.fa chr1:1000-2000
./target/release/fasta-util get genome.fa 1000-2000
./target/release/fasta-util get genome.fa 1000-2000 --fai-index genome.fa.fai
```

### validate

Check FASTA record structure, sequence symbols, duplicate record IDs, line endings within each record, and whether sequence line widths can be indexed with `.fai`. Different records may use different line endings. On success, the command reports the record count and sequence type. On failure, it prints compiler-style diagnostics with the file, line, column, and source line, then exits with status 1.

In the default nucleotide mode, sequences containing `T` are reported as DNA and those containing `U` as RNA. Mixing `T` and `U` is an error. If neither occurs, DNA versus RNA is ambiguous. Select `--sequence-type protein` for protein sequences. The gap symbol `-` is accepted.

```sh
./target/release/fasta-util validate seq.fa
./target/release/fasta-util validate proteins.faa --sequence-type protein
```

### index

Stream through a FASTA file once and create a `.fai` index for random access. The output path is the input filename with `.fai` appended. If the input layout is invalid, index creation stops and any existing index is preserved.

```sh
./target/release/fasta-util index genome.fa
./target/release/fasta-util get genome.fa 100000001-100001000
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

Compare FASTA file-scanning paths with this benchmark. It generates LF and CRLF multi-FASTA inputs and measures an mmap borrowed-line visitor, the mmap line iterator, and `BufReader::read_until`. It checks that each path reads the same number of lines and bytes before timing.

```sh
cargo bench --bench fasta_io
cargo bench --bench fasta_io -- --input-size 1000000 --sample-ms 300
```

The default is 10,000 sequence bases. Each method has three warmup scans followed by five 200 ms samples. `--input-size` counts sequence bases. Repeated scans measure data in the OS page cache, not cold-disk throughput.

The CLI benchmark creates a deterministic random multi-FASTA and compares file and standard-input modes for `stats`, `filter`, and `revcomp`. It also measures the complete `filter | revcomp | stats` pipeline. Before timing, it checks that file and stdin modes produce identical output. Stable Rust, `hyperfine`, and `awk` are required. By default, it measures a 1-million-base input with one warmup and five runs per case.

```sh
./scripts/benchmark_pipeline.sh
./scripts/benchmark_pipeline.sh 10000000 7
```

Stdin timings include staging the input in a temporary file; pipeline timings also include process startup. The script creates the input in a temporary directory and removes it on exit.

## Real-data benchmarks

Measurements were taken on 2026-10-04 using the RefSeq GRCh38.p14 FASTA in `dataset/ncbi_dataset`. The full-genome length comparison used `GCF_000001405.40_GRCh38.p14_genomic.fna` (3,339,739,109 bytes), containing 705 records. Slice comparisons used its chromosome 1 record (`NC_000001.11`, 248,956,422 bases).

The benchmark used the source FASTA directly and preserved its lowercase soft-masked sequence. Before timing each case, regular and `.fai` indexed outputs were verified byte-for-byte. If `seqret` is installed, its output is also verified and included in the comparison. The extraction timings are historical values measured with the former `slice` command.

The regular and `.fai` paths were remeasured by the benchmark script on 2026-10-04 using the same CPU, OS, stable Rust, and hyperfine setup. Their outputs were compared byte-for-byte before five timed runs. `seqret` was unavailable during this remeasurement, so only its values in the table are from the earlier run. Hyperfine reported a statistical outlier for the regular path extracting 100,000 bases from the middle; that result has high variance.

Remeasurement environment: Ubuntu 26.04.1 LTS (WSL2), AMD Ryzen 9 9900X, stable Rust 1.99.0, seqkit 2.10.1, and hyperfine 1.20.0. The seqret values in the table are from the earlier run using EMBOSS seqret 6.6.0.0. Each command had one warmup followed by five measured runs. Tables show the mean and standard deviation. Results vary with the machine and file-cache state.

### len

`seqkit stats` computes statistics for all 705 records; this tool counts sequence symbols. Their total lengths matched.

| Command | Result | Time (mean ± standard deviation) |
| --- | ---: | ---: |
| `seqkit stats` | 3,298,430,636 bases, 705 records | 2,987 ± 247 ms |
| `fasta-util len` | 3,298,430,636 bases | 3,097 ± 487 ms |

### get

Ranges are positions within chromosome 1. The former `slice` command used zero-based inclusive ranges. The current `get` command expresses the same regions with 1-based inclusive `START-END` coordinates. Output was wrapped at 60 bases per line.

| Offset | Slice length | seqret (mean ± standard deviation) | fasta-util (regular, mean ± standard deviation) | fasta-util (FAI, mean ± standard deviation) |
| ---: | ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 843 ± 46 ms | 245.8 ± 32.2 ms | 117.9 ± 22.0 ms |
| 100,000,000 | 100,000 | 745 ± 39 ms | 146.3 ± 33.3 ms | 1.7 ± 0.2 ms |
| 100,000,000 | 100 | 731 ± 113 ms | 130.4 ± 22.8 ms | 1.7 ± 0.2 ms |
| 0 | 100,000,000 | 866 ± 43 ms | 117.6 ± 12.4 ms | 70.3 ± 5.2 ms |
| 0 | 100,000 | 769 ± 113 ms | 1.7 ± 0.1 ms | 1.7 ± 0.3 ms |
| 0 | 100 | 618 ± 137 ms | 1.6 ± 0.1 ms | 1.5 ± 0.1 ms |

To reproduce these measurements, install `awk`, `seqkit`, `hyperfine`, and stable Rust, then run this script from the repository root. `seqret` is optional and adds an external-tool comparison when installed. An alternate FASTA path can be supplied as the first argument.

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```

## License

This project is available under either the [MIT License](LICENSE-MIT) or the [Apache License 2.0](LICENSE-APACHE), at your option. The MIT License copyright holder is `namba3 (GitHub: @namba3)`.
