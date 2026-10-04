# FASTA-Util

[日本語](README.md) | [English](README.en.md)

FASTA ファイルを扱うための CLI ツールです。

## FASTA形式について

FASTAは、塩基配列やアミノ酸配列をテキストで表す形式です。各レコードは`>`で始まるヘッダー行から始まり、その次の行以降に配列を記述します。配列は複数行に折り返して記述できます。

このツールが処理するのは核酸配列のFASTAです。アミノ酸配列には対応していません。

```fasta
>record-1 optional description
ACGTN
UKS-
>record-2
MRY
```

このツールの`len`はヘッダーと空行を除いて配列記号を数えます。`slice`の位置もヘッダーや行区切りを除いた配列上の位置です。複数レコードの場合、配列をファイル順に連結した位置で範囲を指定します。出力には範囲の終端までに現れたヘッダーを残すため、選択範囲の配列記号がないレコードのヘッダーも含まれることがあります。対応する配列記号は大文字・小文字の`ACGTNUKSYMWRBDHV`とギャップを表す`-`です。`slice`は配列記号の大文字・小文字を維持します。それ以外の記号は受け付けません。

## ビルド

```sh
cargo build --release
```

## テストデータの生成

```sh
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 > test.fna
```

同じ乱数シードでデータを再生成するには、`--seed`を指定します。出力の一致は、同じ`rand`バージョンと実行環境で保証されます。
`--line-width`で配列行の長さを変更できます。既定値は50です。

```sh
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 --seed 42 > test.fna
cargo run --manifest-path=generate_random_data/Cargo.toml -- 10000 --line-width 60 > test.fna
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

`--range`は0始まりのRust範囲記法です。`2..10`は位置2〜9、`2..=10`は位置2〜10を選択します。`..10`は先頭から位置9まで、`2..`は位置2から末尾までを選択します。既定値の`..`は配列全体です。`--chars-per-line`は出力時の折り返し幅で、既定値は60です。

`-o`/`--output`で出力ファイルを指定すると、処理が成功した場合にだけ出力先を置き換えます。入力ファイルと同じファイルは出力先に指定できません。

大きな非圧縮FASTAの中ほどから切り出す場合は、対応する`.fai`インデックスを作成して`--fai-index`に渡すと、選択範囲を直接読み込めます。インデックス作成には`samtools faidx`を使えます。FASTAを変更した場合はインデックスを作り直してください。この経路では選択範囲の配列記号を検証します。

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

## 実データのベンチマーク

2026-10-04 に `dataset/ncbi_dataset` 内の RefSeq GRCh38.p14 FASTA を使って測定しました。全ゲノムの長さ計測には705レコードを含む `GCF_000001405.40_GRCh38.p14_genomic.fna`（3,339,739,109 bytes）を使い、slice にはその中の chr1（`NC_000001.11`、248,956,422塩基）を使いました。

元データに含まれる小文字のソフトマスク配列を維持したまま測定しました。slice の各ケースでは、計測前に通常の`slice`と`.fai`を使った`slice`の出力がバイト単位で一致することを確認しています。`seqret`がインストールされている場合は、その出力との一致も確認して比較に加えます。

通常経路と`.fai`経路は2026-10-04に同じCPU・OS・stable Rust・hyperfine環境でベンチマークスクリプトから再測定し、出力一致を確認してから各5回計測しました。`seqret`がなかったため、seqretの値のみ前回測定の結果です。先頭から100塩基を切り出すFAI測定では外れ値の警告が出たため、通常経路との差は測定誤差として扱ってください。

測定環境は Ubuntu 26.04.1 LTS（WSL2）、AMD Ryzen 9 9900X、stable Rust 1.99.0、seqkit 2.10.1、EMBOSS seqret 6.6.0.0、hyperfine 1.20.0 です。各コマンドを1回ウォームアップした後、5回測定した平均と標準偏差を示します。結果は実行環境やファイルキャッシュの状態で変わります。

### len

`seqkit stats` は705レコードの統計を計算し、本ツールは配列記号の総数を計算します。両方の合計値は一致しました。

| コマンド | 結果 | 時間（平均 ± 標準偏差） |
| --- | ---: | ---: |
| `seqkit stats` | 3,298,430,636塩基、705レコード | 1,618 ± 97 ms |
| `fasta-util len` | 3,298,430,636塩基 | 2,646 ± 284 ms |

### slice

切り出し範囲は chr1 内の位置です。`seqret` は1始まりの両端包含、本ツールは0始まりの両端包含で指定し、出力は60塩基ごとに折り返しました。

| オフセット | 切り出し長 | seqret（平均 ± 標準偏差） | fasta-util（通常、平均 ± 標準偏差） | fasta-util（FAI、平均 ± 標準偏差） |
| ---: | ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 843 ± 46 ms | 692.4 ± 51.3 ms | 69.5 ± 9.4 ms |
| 100,000,000 | 100,000 | 745 ± 39 ms | 548.1 ± 293.6 ms | 1.4 ± 0.1 ms |
| 100,000,000 | 100 | 731 ± 113 ms | 434.5 ± 60.4 ms | 1.6 ± 0.2 ms |
| 0 | 100,000,000 | 866 ± 43 ms | 384.1 ± 40.5 ms | 62.8 ± 2.1 ms |
| 0 | 100,000 | 769 ± 113 ms | 2.8 ± 0.4 ms | 1.7 ± 0.2 ms |
| 0 | 100 | 618 ± 137 ms | 2.1 ± 0.2 ms | 1.9 ± 0.3 ms |

同じ条件で再測定するには、`awk`、`seqkit`、`hyperfine` と stable Rust を用意し、リポジトリのルートで次を実行します。`seqret`は任意で、インストールされている場合は外部ツール比較も行います。データセットのパスは引数でも指定できます。

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```
