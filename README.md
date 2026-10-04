# FASTA-Util

[日本語](README.md) | [English](README.en.md)

FASTA ファイルを扱うための CLI ツールです。

## FASTA形式について

FASTAは、塩基配列やアミノ酸配列をテキストで表す形式です。各レコードは`>`で始まるヘッダー行から始まり、その次の行以降に配列を記述します。配列は複数行に折り返して記述できます。

```fasta
>record-1 optional description
ACGTN
UKS-
>record-2
MRY
```

このツールの`len`はヘッダーと空行を除いて配列記号を数えます。`slice`の位置もヘッダーや行区切りを除いた配列上の位置です。対応する配列記号は大文字の`ACGTNUKSYMWRBDHV`とギャップを表す`-`です。小文字やそれ以外の記号は受け付けません。

## ビルド

```sh
cargo build --release
```

## テストデータの生成

```sh
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 > test.fna
```

## サブコマンド

### len

配列の総塩基数を数えます。

```sh
./target/release/fasta-util len -i test.fna
```

```txt
10000
```

### slice

配列の一部を切り出します。

```sh
./target/release/fasta-util slice -i test.fna --range 99..=199
```

```txt
>TestData 10000 random data
WDCAGVUTRABAKRRNRNHHKTYDNBNTCHMRBRRYHWHKYBHKSBAHVNTCGUMGCMMA
GYMDSVCYRAMWNURRVTCYYCYCWWHTRCAUVSBUVHMHNWTGKGHGATWMHYTWNSUB
SUDKUGDWWTSSYBUCKYUDSAADMMRHMT
```

## 簡単なテストとベンチマーク

塩基判定方式のベンチマークを実行するには、次のコマンドを使います。
有効な塩基のみ、無効な文字のみ、有効率50%の混在、有効率99%の混在の4パターンを計測します。
各方式を200ミリ秒ずつ5回計測し、中央値を表示します。
計測サンプルごとに方式の実行順を入れ替え、順序による偏りを抑えます。

```sh
cargo bench --bench nucleic_acid
```

入力サイズと各サンプルの計測時間は、次のように変更できます。標準では10,000バイト、200ミリ秒です。

```sh
cargo bench --bench nucleic_acid -- --input-size 100000 --sample-ms 500
```

OS: Ubuntu (WSL2)

CPU: AMD Ryzen 9 5900X

テストデータ: [Homo sapiens](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/9606/)

### len

`seqkit stats` コマンドと比較します。

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

| コマンド | 時間 (ms) |
| --- | ---: |
| `seqkit stats ncbi_dataset/data/GCF_000001405.39/chr1.fna` | 301.8 |
| `./target/release/fasta-util len -i ncbi_dataset/data/GCF_000001405.39/chr1.fna` | 248.2 |

### slice

`seqret` コマンドと比較します。

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

| 対象 | オフセット | 切り出し長 | コマンド | 時間 (ms) |
| --- | ---: | ---: | --- | ---: |
| seqret | 100,000,000 | 100,000,000 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -sbegin 100000000 -send 200000000` | 665.7 |
| seqret | 100,000,000 | 100,000 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -sbegin 100000000 -send 100100000` | 563.0 |
| seqret | 100,000,000 | 100 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -sbegin 100000000 -send 100000100` | 552.7 |
| seqret | 0 | 100,000,000 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -send 100000000` | 677.1 |
| seqret | 0 | 100,000 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -send 100000` | 565.2 |
| seqret | 0 | 100 | `seqret -sequence ncbi_dataset/data/GCF_000001405.39/chr1.fna -auto -stdout -send 100` | 558.3 |
| fasta-util | 100,000,000 | 100,000,000 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=199999999` | 251.9 |
| fasta-util | 100,000,000 | 100,000 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=100099999` | 127.6 |
| fasta-util | 100,000,000 | 100 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range 99999999..=100000099` | 101.7 |
| fasta-util | 0 | 100,000,000 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range ..=99999999` | 217.1 |
| fasta-util | 0 | 100,000 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range ..=99999` | 0.9 |
| fasta-util | 0 | 100 | `./target/release/fasta-util slice -i ncbi_dataset/data/GCF_000001405.39/chr1.fna --chars-per-line=60 --range ..=99` | 0.7 |
