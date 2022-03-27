# FASTA-Util

A CLI tool for playing with FASTA files

## Build

```sh
cargo build --release
```

## Generate Test Data

```sh
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 > test.fna
```

## Count the total length of the sequence

```sh
./target/release/fasta-util len -i test.fna
```

```txt
10000
```

## Cut out a part of the sequence

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

OS: Ubuntu (WSL2)
CPU: AMD Ryzen 9 5900X
Test Data: [Homo sapiens](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/9606/)

### len

Compare with `seqkit stats` command

seqkit:

```sh
seqkit stats ncbi_dataset/data/GCF_000001405.39/chr1.fna
```

```txt
file                                         format  type  num_seqs      sum_len      min_len      avg_len      max_len
ncbi_dataset/data/GCF_000001405.39/chr1.fna  FASTA   DNA          1  248,956,422  248,956,422  248,956,422  248,956,422
```

fasta-util:

```txt
./target/release/fasta-util len -i ncbi_dataset/data/GCF_000001405.39/chr1.fna
```

```txt
248956422
```

| command                                                                          | time (ms) |
| -------------------------------------------------------------------------------- | --------- |
| `seqkit stats ncbi_dataset/data/GCF_000001405.39/chr1.fna`                       | 301.8     |
| `./target/release/fasta-util len -i ncbi_dataset/data/GCF_000001405.39/chr1.fna` | 248.2     |

### slice

Compare with `seqret` command

seqret:

```sh
seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -sbegin 100000000 -send 200000000 -auto -stdout | sha256sum
```

```txt
c570cb67eb05a25922a3fc6f299cdc8bb5763ae3375281505bf15ff0537286cf  -
```

fasta-util:

```sh
./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=199999999 | sha256sum
```

```txt
c570cb67eb05a25922a3fc6f299cdc8bb5763ae3375281505bf15ff0537286cf  -
```

| command                                                                                                                             | time (ms) |
| ----------------------------------------------------------------------------------------------------------------------------------- | --------- |
| `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout`                                                        | 835.9     |
| `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -sbegin 100000000 -send 200000000 -auto -stdout`                      | 747.2     |
| `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -sbegin 200000 -send 300000 -auto -stdout`                            | 610.0     |
| `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60`                              | 366.2     |
| `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=199999999` | 283.1     |
| `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 199999..=299999`      | 1.2       |
