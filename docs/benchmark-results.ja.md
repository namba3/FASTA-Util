# ベンチマーク結果

2026-10-08に再測定した結果を記録します。全サブコマンド、標準入力とパイプ、RefSeq GRCh38.p14データセット、Rust内のマイクロベンチマークを含みます。以下の値はこの環境での測定結果で、性能を保証するものではありません。

English: [benchmark-results.md](benchmark-results.md).

## 測定環境

- Ubuntu 26.04.1 LTS（WSL2）
- AMD Ryzen 9 9900X（WSL2から11コア・22論理CPUを認識）
- stable Rust 1.99.0、Cargo 1.99.0
- hyperfine 1.20.0、seqkit 2.10.1
- EMBOSS `seqret`は未導入

## CLIサブコマンド

stable Rust、`hyperfine`、`awk`を用意し、リポジトリルートから実行します。

```sh
./scripts/benchmark_commands.sh
# 塩基数と計測回数は変更できます
./scripts/benchmark_commands.sh 50000000 7
```

記録した測定では2000万記号、各ケース5回の計測を行い、その前に1回ウォームアップしました。`--seed 42`で生成した核酸FASTAは、5万塩基ずつの400レコードで、ファイルサイズは20,410,800 bytesでした。タンパク質FASTAは決定的な配列20,000,000記号、20,333,368 bytesです。出力は`/dev/null`へ捨てていますが、プロセス起動時間は測定に含まれます。標準入力ケースには`cat`とパイプの時間も含まれます。

計測前に、`stats`・`composition`・`filter`・`revcomp`・`format`・`validate`のファイル入力と標準入力で出力が一致することを確認しました。`get`ではFAIの有無による出力一致も確認し、核酸・タンパク質の検証と統計コマンドを実行しました。対象範囲は後半レコードの`record_000201:1001-11000`です。`get`のFAIは計測前に作成し、独立した`index`ケースでは毎回FAIを書き込みます。

以下はhyperfineによる平均 ± 標準偏差です。FAIを使う`get`はプロセス起動時間に近いため、小さな時間差は参考値として扱ってください。

| コマンド | 平均 ± 標準偏差 |
| --- | ---: |
| `len`（核酸） | 14.7 ± 1.1 ms |
| `validate`（核酸） | 103.1 ± 2.4 ms |
| `validate`（タンパク質） | 26.5 ± 1.1 ms |
| `index`（400レコード） | 17.6 ± 3.0 ms |
| `stats`（核酸） | 117.9 ± 1.0 ms |
| `stats`（タンパク質） | 44.2 ± 3.7 ms |
| `composition`（核酸） | 48.1 ± 1.2 ms |
| `composition`（タンパク質） | 45.6 ± 3.0 ms |
| `get`（FAIなし、10 kb） | 13.7 ± 0.7 ms |
| `get`（FAIあり、10 kb） | 1.3 ± 0.3 ms |
| `filter`（核酸） | 84.7 ± 4.3 ms |
| `filter`（タンパク質） | 49.8 ± 3.3 ms |
| `revcomp` | 224.9 ± 10.7 ms |
| `grep` | 23.2 ± 5.1 ms |
| `locate`（4記号モチーフ、不一致1を許容） | 668.9 ± 23.2 ms |
| `format`（幅60） | 41.5 ± 3.3 ms |

| 標準入力・パイプ | 平均 ± 標準偏差 |
| --- | ---: |
| `stats`（stdin） | 183.4 ± 7.9 ms |
| `composition`（stdin） | 91.1 ± 2.7 ms |
| `filter`（stdin） | 115.6 ± 6.6 ms |
| `revcomp`（stdin） | 259.0 ± 9.8 ms |
| `format`（stdin） | 80.2 ± 7.1 ms |
| `validate`（stdin） | 169.3 ± 6.8 ms |
| `filter | revcomp | stats` | 390.6 ± 16.2 ms |

この測定では`grep`のばらつきが他の多くのケースより大きくなりました。FAIを使う`get`は5 ms未満で、hyperfineからシェル起動時間の校正精度について警告が出ています。

`BASES`は少なくとも10,011,000を指定してください。測定対象の領域が201番目のレコード内に必要です。

## RefSeq GRCh38.p14データセット

`dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna`を再測定に使用しました。ファイルは3,339,739,109 bytes、705レコードで、配列記号の合計は3,298,430,636です。第1染色体レコード`NC_000001.11`は248,956,422塩基です。小文字のソフトマスク配列は維持しました。計測前に`get`の通常経路とFAI経路の出力がバイト単位で一致することを確認しました。各ケースは1回ウォームアップした後、5回計測しています。`seqkit stats`と`fasta-util len`の総配列長は一致しました。

| コマンド | 結果 | 平均 ± 標準偏差 |
| --- | ---: | ---: |
| `seqkit stats`（全ゲノム） | 705レコード、3,298,430,636塩基 | 1.634 ± 0.056 s |
| `fasta-util len`（全ゲノム） | 3,298,430,636塩基 | 1.563 ± 0.045 s |

以下は第1染色体からの`get`の結果です。オフセットは範囲を説明するための0始まりの値で、コマンドでは1始まり・両端を含む座標を使います。出力は1行60塩基で折り返しています。

| オフセット | 長さ | FAIなし | FAIあり |
| ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 48.1 ± 0.7 ms | 48.8 ± 1.6 ms |
| 100,000,000 | 100,000 | 1.2 ± 0.1 ms | 1.2 ± 0.2 ms |
| 100,000,000 | 100 | 1.2 ± 0.3 ms | 1.1 ± 0.1 ms |
| 0 | 100,000,000 | 64.8 ± 4.9 ms | 53.6 ± 2.8 ms |
| 0 | 100,000 | 1.5 ± 0.4 ms | 1.4 ± 0.1 ms |
| 0 | 100 | 1.1 ± 0.1 ms | 1.1 ± 0.1 ms |

短い範囲の計測はプロセス起動と温まったファイルキャッシュの影響が大きくなります。先頭から100 Mbを切り出すケースではFAIによる改善が確認できましたが、中央からの100 Mb切り出しでは差がほぼありませんでした。`seqret`がなかったため、外部ツールとの比較は含みません。

同梱データセットを使うか、別のFASTAを第1引数に渡して再測定できます。

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```

このスクリプトには`seqkit`、`hyperfine`、`awk`、stable Rustが必要です。元データを読み、作業用ファイルはOSの一時ディレクトリに作ります。

## Rustマイクロベンチマーク

これらはCLI全体ではなく、内部処理方式を測ります。

```sh
cargo bench --bench nucleic_acid
cargo bench --bench fasta_io
```

`nucleic_acid`は各パターン1万bytesで、200 msの測定を5回行った中央値です。今回の記号パターンではlookup tableが最も速い結果でした。

| 入力パターン | `match` | Set iteration | Lookup table |
| --- | ---: | ---: | ---: |
| すべて有効 | 0.298 ns/base | 3.099 ns/base | 0.177 ns/base |
| すべて無効 | 0.286 ns/base | 1.821 ns/base | 0.170 ns/base |
| 有効50% | 0.293 ns/base | 3.055 ns/base | 0.172 ns/base |
| 有効99% | 0.292 ns/base | 4.882 ns/base | 0.173 ns/base |

`fasta_io`はLF・CRLFの入力をページキャッシュ上で繰り返し走査しました。入力は塩基1万個、各方式3回ウォームアップ後に200 msの測定を5回行っています。

| 改行 | mmap borrowed visitor | mmap iterator | buffered `read_until` |
| --- | ---: | ---: | ---: |
| LF | 723.57 MiB/s | 706.34 MiB/s | 2,101.32 MiB/s |
| CRLF | 700.68 MiB/s | 719.71 MiB/s | 2,126.95 MiB/s |

この結果は温まったページキャッシュ上の反復走査を示し、ストレージから初回に読む速度ではありません。コンパイラ、CPU、OS、バックグラウンド負荷で値は変わります。数値を比較する際は同じ条件で再測定してください。

## 用途別ベンチマークスクリプト

個別の条件を調べるスクリプトも利用できます。

- `./scripts/benchmark_pipeline.sh [BASES] [RUNS]`は`stats`・`filter`・`revcomp`のファイル入力とstdin、および`filter | revcomp | stats`全体を比較します。
- `./scripts/benchmark_analysis.sh [BASES] [RUNS]`は`stats`のレコード数と`locate`のモチーフ長・不一致許容数・一致頻度を変えて測ります。
- `./scripts/benchmark_real_data.sh [FASTA]`は全体の`len`と`seqkit stats`を比較し、第1染色体のFAI有無による`get`を測ります。`seqkit`が必要で、`seqret`は任意です。

CLIスクリプトの既定値は`benchmark_commands.sh`が2000万塩基、その他が100万塩基で、計測回数は5回です。`cargo bench`のサンプル設定はそれぞれの`-- --help`を参照してください。
