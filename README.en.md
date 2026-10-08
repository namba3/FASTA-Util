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

| Command | Purpose |
| --- | --- |
| [`len`](#len) | Count sequence symbols |
| [`validate`](#validate) | Validate FASTA structure and sequence symbols |
| [`index`](#index) | Create a `.fai` random-access index |
| [`stats`](#stats) | Summarize lengths, N50, GC, and related metrics |
| [`composition`](#composition) | Report sequence symbol frequencies |
| [`get`](#get) | Retrieve records or sequence ranges |
| [`filter`](#filter) | Select records by length, GC, or N fraction |
| [`revcomp`](#revcomp) | Reverse-complement nucleotide sequences |
| [`grep`](#grep) | Search records by header text |
| [`locate`](#locate) | Find motif positions in sequences |
| [`format`](#format) | Reformat line widths, letter case, gaps, and headers |

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
./target/release/fasta-util len -i genome.fa --threads 4
```

File inputs can use multiple workers for length counting. Set the worker count with `--threads N`. When omitted, files smaller than 8 MiB are processed sequentially; larger files use up to 4 workers, limited by available CPUs. Standard input is processed sequentially.

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
./target/release/fasta-util stats genome.fa --threads 4
```

File inputs accept `--threads N` to set the worker count. When omitted, files smaller than 8 MiB are processed sequentially; larger files use up to 4 workers, limited by available CPUs. Standard input is processed sequentially. Statistics split work at record boundaries, so a file containing one record is processed by one worker.

### composition

Show the fraction of each symbol across all sequence data. For nucleotides, the output includes `A`, `C`, `G`, `T` (or `U` for RNA), `N`, any IUPAC ambiguity symbols or gaps found in the input, and `GC` for the combined `G` and `C` fraction. For proteins, it reports the 20 standard amino acids and any extended symbols found. The denominator includes every sequence symbol, including gaps and ambiguity symbols. Sequence type is detected automatically as in `stats`; specify `--sequence-type nucleotide` or `--sequence-type protein` when needed.

```sh
./target/release/fasta-util composition seq.fa
./target/release/fasta-util composition proteins.faa --sequence-type protein
./target/release/fasta-util composition genome.fa --threads 4
cat seq.fa | ./target/release/fasta-util composition
```

File inputs accept `--threads N` to set the worker count. When omitted, files smaller than 8 MiB are processed sequentially; larger files use up to 4 workers, limited by available CPUs. Standard input is processed sequentially.

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

Get complete records by ID, record regions using `ID:START-END`, or a global range using `START-END`. Coordinates are 1-based and inclusive. Global ranges index all record sequences concatenated in FASTA order. Multiple IDs and an ID file passed with `--ids` are supported; output follows the records' order in the FASTA file. If a `.fai` sidecar exists next to the input, it is used automatically; otherwise, the file is scanned from the beginning. Pass `--fai-index` to select an index, or `--no-fai-index` to disable automatic index use. Set output wrapping with `--chars-per-line`.

```sh
./target/release/fasta-util get genome.fa chr1
./target/release/fasta-util get genome.fa chr1 chr3 chrX
./target/release/fasta-util get genome.fa --ids chromosomes.txt
./target/release/fasta-util get genome.fa chr1:1000-2000
./target/release/fasta-util get genome.fa 1000-2000
./target/release/fasta-util get genome.fa 1000-2000 --fai-index genome.fa.fai
./target/release/fasta-util get genome.fa 1000-2000 --no-fai-index
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

## Benchmarks

The procedures and results cover every subcommand, standard input and pipelines, comparisons with SeqKit and seqret, a real genome, and internal Rust implementations in [`docs/benchmark-results.md`](docs/benchmark-results.md). CLI benchmarks use `hyperfine`, `awk`, and stable Rust; the default workload is 20 million bases, with one warmup and five measured runs per case.

```sh
./scripts/benchmark_commands.sh
./scripts/benchmark_commands.sh 50000000 7
./scripts/benchmark_comparison.sh
cargo bench --bench nucleic_acid
cargo bench --bench fasta_io
./scripts/benchmark_real_data.sh
```

Results vary with the CPU, OS, file-cache state, and system load. Process startup time contributes to short command timings.

## License

This project is available under either the [MIT License](LICENSE-MIT) or the [Apache License 2.0](LICENSE-APACHE), at your option. The MIT License copyright holder is `namba3 (GitHub: @namba3)`.
