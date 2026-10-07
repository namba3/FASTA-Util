# Benchmark results

This report records a benchmark rerun on 2026-10-08. It includes every `fasta-util` subcommand, file and standard-input paths, a multi-command pipeline, the RefSeq GRCh38.p14 dataset, and the Rust microbenchmarks. Use the results as measurements of this environment; they are not performance guarantees.

Japanese report: [benchmark-results.ja.md](benchmark-results.ja.md).

## Environment

- Ubuntu 26.04.1 LTS under WSL2
- AMD Ryzen 9 9900X, 11 cores / 22 logical CPUs reported by WSL2
- stable Rust 1.99.0, Cargo 1.99.0
- hyperfine 1.20.0; seqkit 2.10.1
- EMBOSS `seqret` was unavailable

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

The remeasurement used `dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna` (3,339,739,109 bytes, 705 records, 3,298,430,636 sequence symbols). The chromosome 1 record `NC_000001.11` has 248,956,422 bases. The benchmark preserved lowercase soft masking. `get` output with and without `.fai` was compared byte-for-byte before timing. Each case had one warmup and five measured runs. `seqkit stats` and `fasta-util len` reported the same total sequence length.

| Command | Result | Mean ± σ |
| --- | ---: | ---: |
| `seqkit stats` (whole assembly) | 705 records, 3,298,430,636 bases | 1.634 ± 0.056 s |
| `fasta-util len` (whole assembly) | 3,298,430,636 bases | 1.563 ± 0.045 s |

`get` results below are for chromosome 1. The offset is zero-based within the input for describing the selected region; the command uses 1-based inclusive coordinates. Output was wrapped at 60 bases per line.

| Offset | Length | Scan | FAI |
| ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 48.1 ± 0.7 ms | 48.8 ± 1.6 ms |
| 100,000,000 | 100,000 | 1.2 ± 0.1 ms | 1.2 ± 0.2 ms |
| 100,000,000 | 100 | 1.2 ± 0.3 ms | 1.1 ± 0.1 ms |
| 0 | 100,000,000 | 64.8 ± 4.9 ms | 53.6 ± 2.8 ms |
| 0 | 100,000 | 1.5 ± 0.4 ms | 1.4 ± 0.1 ms |
| 0 | 100 | 1.1 ± 0.1 ms | 1.1 ± 0.1 ms |

The small-region timings are dominated by process startup and warm file-cache behavior. The larger extraction from the beginning showed a measurable FAI improvement; the middle 100 Mb extraction was effectively the same in this run. `seqret` was unavailable, so this report does not include an external extraction comparison.

Reproduce with the checked-in dataset, or pass another FASTA path as the first argument:

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```

The script requires `seqkit`, `hyperfine`, `awk`, and stable Rust. It uses the source dataset and makes temporary benchmark files under the system temporary directory.

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
- `./scripts/benchmark_real_data.sh [FASTA]` compares whole-file `len` with `seqkit stats` and measures `get` with and without `.fai` on chromosome 1. `seqkit` is required; `seqret` is optional.

The script arguments for CLI suites default to 1,000,000 bases and five runs, except `benchmark_commands.sh`, which defaults to 20,000,000 bases. The `cargo bench` programs take their own sampling options; see `-- --help` for the current arguments.
