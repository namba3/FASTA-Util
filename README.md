# FASTA-Util

A CLI tool for playing with FASTA files

## Build

```sh
cargo +nightly build --release
```

## Generate Test Data

```sh
cargo +nightly run --manifest-path=generate_random_data/Cargo.toml -- 10000 > test.fna
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

| target     | offset      | slice size  | command                                                                                                                             | time (ms) |
| ---------- | ----------- | ----------- | ----------------------------------------------------------------------------------------------------------------------------------- | --------- |
| seqret     | 100,000,000 | 100,000,000 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -sbegin 100000000 -send 200000000`                      | 665.7     |
| seqret     | 100,000,000 | 100,000     | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -sbegin 100000000 -send 100100000`                      | 563.0     |
| seqret     | 100,000,000 | 100         | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -sbegin 100000000 -send 100000100`                      | 552.7     |
| seqret     | 0           | 100,000,000 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -send 100000000`                                        | 677.1     |
| seqret     | 0           | 100,000     | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -send 100000`                                           | 565.2     |
| seqret     | 0           | 100         | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -send 100`                                              | 558.3     |
| fasta-util | 100,000,000 | 100,000,000 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=199999999` | 251.9     |
| fasta-util | 100,000,000 | 100,000     | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=100099999` | 127.6     |
| fasta-util | 100,000,000 | 100         | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=100000099` | 101.7     |
| fasta-util | 0           | 100,000,000 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range ..=99999999`          | 217.1     |
| fasta-util | 0           | 100,000     | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range ..=99999`             | 0.9       |
| fasta-util | 0           | 100         | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range ..=99`                | 0.7       |
