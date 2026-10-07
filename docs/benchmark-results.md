# Benchmark results

This report records benchmark reruns on 2026-10-08. It includes every `fasta-util` subcommand, comparisons with SeqKit and EMBOSS `seqret` where functions overlap, file and standard-input paths, a multi-command pipeline, the RefSeq GRCh38.p14 dataset, and Rust microbenchmarks. Use the results as measurements of this environment; they are not performance guarantees.

Japanese report: [benchmark-results.ja.md](benchmark-results.ja.md).

## Environment

- Ubuntu 26.04.1 LTS under WSL2
- AMD Ryzen 9 9900X, 11 cores / 22 logical CPUs reported by WSL2
- stable Rust 1.99.0, Cargo 1.99.0
- hyperfine 1.20.0; seqkit 2.10.1
- EMBOSS `seqret` 6.6.0.0 from the Ubuntu 26.04 package bundle in `dataset/emboss/`

## CLI command benchmark

Run from the repository root with stable Rust, `hyperfine`, and `awk` installed:

```sh
./scripts/benchmark_commands.sh
# Optional: set total nucleotide/protein symbols and measured runs
./scripts/benchmark_commands.sh 50000000 7
```

The recorded run used 20,000,000 symbols and five timed runs per case, following one warmup. It generated a seeded nucleotide FASTA (`--seed 42`) with 400 records of 50,000 bases each, plus a deterministic protein FASTA of 20,000,000 symbols. File sizes were 20,410,800 bytes and 20,333,368 bytes respectively. Commands discarded output to `/dev/null`; measured time includes process startup. Standard-input cases also include `cat` and pipe overhead.

`BASES` must be at least 10,011,000 so the benchmark region exists in record 201.

Before measurement the script checked file/stdin output equality for `stats`, `composition`, `filter`, `revcomp`, `format`, and `validate`; checked that `get` returned identical bytes with and without an index; and exercised nucleotide and protein validation/statistics. The `get` region was `record_000201:1001-11000`, which is 10 kb inside a later record. The index was created before timing `get`; the standalone `index` case writes a fresh sidecar on each run.

Times below are mean ± standard deviation from hyperfine. The `get` index result is near the process-startup floor, so treat small absolute differences cautiously.

| Command | Mean ± σ |
| --- | ---: |
| `len` (nucleotide) | 14.7 ± 1.1 ms |
| `validate` (nucleotide) | 103.1 ± 2.4 ms |
| `validate` (protein) | 26.5 ± 1.1 ms |
| `index` (400 records) | 17.6 ± 3.0 ms |
| `stats` (nucleotide) | 117.9 ± 1.0 ms |
| `stats` (protein) | 44.2 ± 3.7 ms |
| `composition` (nucleotide) | 48.1 ± 1.2 ms |
| `composition` (protein) | 45.6 ± 3.0 ms |
| `get` (scan, 10 kb) | 13.7 ± 0.7 ms |
| `get` (FAI, 10 kb) | 1.3 ± 0.3 ms |
| `filter` (nucleotide) | 84.7 ± 4.3 ms |
| `filter` (protein) | 49.8 ± 3.3 ms |
| `revcomp` | 224.9 ± 10.7 ms |
| `grep` | 23.2 ± 5.1 ms |
| `locate` (4-base motif, 1 mismatch allowed) | 668.9 ± 23.2 ms |
| `format` (width 60) | 41.5 ± 3.3 ms |

| Standard-input / pipeline case | Mean ± σ |
| --- | ---: |
| `stats` (stdin) | 183.4 ± 7.9 ms |
| `composition` (stdin) | 91.1 ± 2.7 ms |
| `filter` (stdin) | 115.6 ± 6.6 ms |
| `revcomp` (stdin) | 259.0 ± 9.8 ms |
| `format` (stdin) | 80.2 ± 7.1 ms |
| `validate` (stdin) | 169.3 ± 6.8 ms |
| `filter | revcomp | stats` | 390.6 ± 16.2 ms |

`grep` had higher run-to-run variation than most cases in this sample. The indexed `get` case is below 5 ms, where hyperfine warns that shell startup calibration limits timing accuracy.

## RefSeq GRCh38.p14 dataset

The remeasurement used `dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna` (3,339,739,109 bytes, 705 records, 3,298,430,636 sequence symbols). The chromosome 1 record `NC_000001.11` has 248,956,422 bases. The benchmark preserved lowercase soft masking. The scan path passed `--no-fai-index`, and the indexed path passed `--fai-index` explicitly; their outputs were compared byte-for-byte before timing. Each case had one warmup and five measured runs. `seqkit stats` and `fasta-util len` reported the same total sequence length.

| Command | Result | Mean ± σ |
| --- | ---: | ---: |
| `seqkit stats` (whole assembly) | 705 records, 3,298,430,636 bases | 1.700 ± 0.083 s |
| `fasta-util len` (whole assembly) | 3,298,430,636 bases | 2.133 ± 0.134 s |

`get` results below are for chromosome 1. The offset is zero-based within the input for describing the selected region; the command uses 1-based inclusive coordinates. Output was wrapped at 60 bases per line.

| Offset | Length | `seqret` | Scan | FAI |
| ---: | ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 846.0 ± 56.2 ms | 175.8 ± 7.8 ms | 56.7 ± 2.5 ms |
| 100,000,000 | 100,000 | 693.8 ± 76.1 ms | 98.0 ± 1.8 ms | 1.7 ± 0.3 ms |
| 100,000,000 | 100 | 646.6 ± 37.1 ms | 81.6 ± 2.9 ms | 1.3 ± 0.2 ms |
| 0 | 100,000,000 | 898.7 ± 64.0 ms | 109.1 ± 10.5 ms | 54.0 ± 2.7 ms |
| 0 | 100,000 | 701.6 ± 44.3 ms | 1.6 ± 0.3 ms | 1.4 ± 0.1 ms |
| 0 | 100 | 735.2 ± 44.6 ms | 1.5 ± 0.1 ms | 1.4 ± 0.1 ms |

The `seqret` output matched `fasta-util get` byte-for-byte in all six cases. For 100 bp and 100 kb at the start of the file, process startup dominates both `get` paths, so the measured difference is small. For small regions in the middle, scanning reads sequence data up to the requested start, so seeking through the FAI makes a large difference. A 100 Mb extraction also spends time copying and writing the selected sequence. The previous scan measurements were invalid because that command automatically detected the adjacent `.fai`; this table was rerun with `--no-fai-index`.

Reproduce with the checked-in dataset, or pass another FASTA path as the first argument:

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```

The script requires `seqkit`, `hyperfine`, `awk`, and stable Rust. It uses the source dataset and makes temporary benchmark files under the system temporary directory. It uses `seqret` from `PATH` when available, otherwise it detects the local bundle at `dataset/emboss/bin/seqret`.

The range comparison uses EMBOSS `seqret`'s `-sbegin` and `-send` options, as described in the [official EMBOSS seqret documentation](https://emboss.sourceforge.net/apps/release/6.4/emboss/apps/seqret.html).

## SeqKit comparison

Reproduce the comparable-operation run with:

```sh
./scripts/benchmark_comparison.sh
./scripts/benchmark_comparison.sh 50000000 7
```

This comparison script requires SeqKit, stable Rust, `hyperfine`, and `awk`.

The recorded run used a deterministic 20-million-base FASTA with 533 records of alternating 25 kb and 50 kb lengths. The generated IUPAC ambiguity symbols were mapped to `A`, leaving canonical DNA, so `locate` measured equivalent literal-motif behavior. Before timing, the script compared selected/filter, reverse-complement, header-search, formatting, locate coordinates/strands, and indexed-region outputs. It also checked that `len` and SeqKit's `sum_len` agreed. All timed output was discarded. Each case used one warmup and five runs.

SeqKit's default parallelism is four threads; `fasta-util` processes these inputs in its streaming path. The commands use each tool's normal defaults otherwise. `stats`/`stats --all` calculate overlapping but not identical metric sets, and `len` versus `seqkit stats` compares only the total-length metric. `validate` and `composition` are not included because there is no direct SeqKit command with the same output contract. SeqKit documents the corresponding `stats`, `seq`, `grep`, `locate`, and `faidx` operations in its [usage guide](https://bioinf.shenwei.me/seqkit/usage/).

| Operation | `fasta-util` | SeqKit | Mean ± σ |
| --- | ---: | ---: | ---: |
| Total length | `len` | `stats` (`sum_len` only) | 14.7 ± 0.5 ms / 36.2 ± 2.6 ms |
| Summary statistics | `stats` | `stats --all` | 80.9 ± 1.2 ms / 63.7 ± 2.7 ms |
| Length filter (min 40 kb) | `filter` | `seq --min-len 40000` | 96.7 ± 4.6 ms / 36.0 ± 2.1 ms |
| Reverse complement | `revcomp` | `seq --reverse --complement` | 146.6 ± 2.6 ms / 106.0 ± 6.3 ms |
| Header search | `grep` | `grep --by-name --use-regexp` | 15.0 ± 1.4 ms / 43.4 ± 1.4 ms |
| Motif location (`ACGA`) | `locate` | `locate` | 127.4 ± 8.6 ms / 133.1 ± 7.3 ms |
| Rewrap to width 80 | `format` | `seq --line-width 80` | 35.9 ± 1.2 ms / 32.8 ± 1.7 ms |
| Indexed extraction (10 kb) | `get` | `faidx` | 1.5 ± 0.3 ms / 18.6 ± 1.5 ms |
| Index creation | `index` | `faidx --update-faidx` | 14.2 ± 0.7 ms / 48.3 ± 3.0 ms |

The indexed extraction result for `fasta-util` is below 5 ms, where hyperfine warns that shell startup calibration limits timing accuracy. Results are environment-sensitive and should not be treated as a universal ranking.

## Rust microbenchmarks

These are implementation-focused benchmarks, distinct from full CLI execution. Run them with:

```sh
cargo bench --bench nucleic_acid
cargo bench --bench fasta_io
```

`nucleic_acid` used 10,000 bytes per input pattern and reported the median throughput from five 200 ms samples. The lookup table was fastest for the tested symbol mixes:

| Input pattern | `match` | Set iteration | Lookup table |
| --- | ---: | ---: | ---: |
| All valid | 0.298 ns/base | 3.099 ns/base | 0.177 ns/base |
| All invalid | 0.286 ns/base | 1.821 ns/base | 0.170 ns/base |
| 50% valid | 0.293 ns/base | 3.055 ns/base | 0.172 ns/base |
| 99% valid | 0.292 ns/base | 4.882 ns/base | 0.173 ns/base |

`fasta_io` repeatedly scanned warm-cache LF and CRLF inputs (10,000 sequence bases, 5 samples after 3 warmups, 200 ms per sample). Throughput was:

| Line endings | mmap borrowed visitor | mmap iterator | buffered `read_until` |
| --- | ---: | ---: | ---: |
| LF | 723.57 MiB/s | 706.34 MiB/s | 2,101.32 MiB/s |
| CRLF | 700.68 MiB/s | 719.71 MiB/s | 2,126.95 MiB/s |

This measures repeated warm-cache scanning, not cold storage throughput. Benchmark figures vary with compiler, CPU, operating system, and background load. The commands above can regenerate local results; timings should be compared only under matching conditions.

## Additional benchmark scripts

The focused scripts remain available for investigating specific workloads:

- `./scripts/benchmark_pipeline.sh [BASES] [RUNS]` compares file and stdin modes for `stats`, `filter`, and `revcomp`, plus the complete `filter | revcomp | stats` pipeline.
- `./scripts/benchmark_analysis.sh [BASES] [RUNS]` varies `stats` record count and `locate` motif length, mismatch allowance, and match frequency.
- `./scripts/benchmark_comparison.sh [BASES] [RUNS]` checks and times operations shared with SeqKit.
- `./scripts/benchmark_real_data.sh [FASTA]` compares whole-file `len` with `seqkit stats` and measures `get` with and without `.fai` on chromosome 1. `seqkit` is required; `seqret` is optional.

The script arguments for CLI suites default to 1,000,000 bases and five runs, except `benchmark_commands.sh`, which defaults to 20,000,000 bases. The `cargo bench` programs take their own sampling options; see `-- --help` for the current arguments.
