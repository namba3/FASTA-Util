# ベンチマーク結果

2026-10-08に再測定した結果を記録します。全サブコマンド、SeqKit・EMBOSS `seqret`との機能比較、標準入力とパイプ、RefSeq GRCh38.p14データセット、Rust内のマイクロベンチマークを含みます。以下の値はこの環境での測定結果で、性能を保証するものではありません。

English: [benchmark-results.md](benchmark-results.md).

## 測定環境

- Ubuntu 26.04.1 LTS（WSL2）
- AMD Ryzen 9 9900X（WSL2から11コア・22論理CPUを認識）
- stable Rust 1.99.0、Cargo 1.99.0
- hyperfine 1.20.0、seqkit 2.10.1
- EMBOSS `seqret` 6.6.0.0（Ubuntu 26.04向けの`dataset/emboss/`内ローカルパッケージ）

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

`dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna`を再測定に使用しました。ファイルは3,339,739,109 bytes、705レコードで、配列記号の合計は3,298,430,636です。第1染色体レコード`NC_000001.11`は248,956,422塩基です。小文字のソフトマスク配列は維持しました。隣接FAIがあっても走査する条件には`--no-fai-index`を指定し、FAI経路は`--fai-index`で明示しました。計測前に両経路の出力がバイト単位で一致することを確認しました。各ケースは1回ウォームアップした後、5回計測しています。`seqkit stats`と`fasta-util len`の総配列長は一致しました。

| コマンド | 結果 | 平均 ± 標準偏差 |
| --- | ---: | ---: |
| `seqkit stats`（全ゲノム） | 705レコード、3,298,430,636塩基 | 1.700 ± 0.083 s |
| `fasta-util len`（全ゲノム） | 3,298,430,636塩基 | 2.133 ± 0.134 s |

以下は第1染色体からの`get`の結果です。オフセットは範囲を説明するための0始まりの値で、コマンドでは1始まり・両端を含む座標を使います。出力は1行60塩基で折り返しています。

| オフセット | 長さ | `seqret` | FAIなし | FAIあり |
| ---: | ---: | ---: | ---: | ---: |
| 100,000,000 | 100,000,000 | 846.0 ± 56.2 ms | 175.8 ± 7.8 ms | 56.7 ± 2.5 ms |
| 100,000,000 | 100,000 | 693.8 ± 76.1 ms | 98.0 ± 1.8 ms | 1.7 ± 0.3 ms |
| 100,000,000 | 100 | 646.6 ± 37.1 ms | 81.6 ± 2.9 ms | 1.3 ± 0.2 ms |
| 0 | 100,000,000 | 898.7 ± 64.0 ms | 109.1 ± 10.5 ms | 54.0 ± 2.7 ms |
| 0 | 100,000 | 701.6 ± 44.3 ms | 1.6 ± 0.3 ms | 1.4 ± 0.1 ms |
| 0 | 100 | 735.2 ± 44.6 ms | 1.5 ± 0.1 ms | 1.4 ± 0.1 ms |

6ケースすべてで`seqret`と`fasta-util get`の出力がバイト単位で一致しました。先頭からの100 bp・100 kbでは走査してもFAIを使ってもプロセス起動時間が支配的で、測定差は小さくなります。中間位置の小範囲では走査が開始位置までの配列を読むため、FAIによるseekの効果が大きく出ました。100 Mb切り出しでは出力・配列コピーのコストも加わります。前回の計測ではFAIなし側のコマンドも隣接`.fai`を自動検出していたため、走査値として無効でした。この表は`--no-fai-index`を指定して再計測した値です。

同梱データセットを使うか、別のFASTAを第1引数に渡して再測定できます。

```sh
./scripts/benchmark_real_data.sh
./scripts/benchmark_real_data.sh path/to/genomic.fna
```

このスクリプトには`seqkit`、`hyperfine`、`awk`、stable Rustが必要です。元データを読み、作業用ファイルはOSの一時ディレクトリに作ります。`seqret`はPATH上の実行ファイルを使い、見つからない場合は`dataset/emboss/bin/seqret`を探します。

範囲切り出し比較では、[EMBOSS公式seqret説明](https://emboss.sourceforge.net/apps/release/6.4/emboss/apps/seqret.html)にある`-sbegin`・`-send`を使います。

## SeqKitとの比較

機能が重なるコマンドの比較は次のスクリプトで再現できます。

```sh
./scripts/benchmark_comparison.sh
./scripts/benchmark_comparison.sh 50000000 7
```

この比較スクリプトにはSeqKit、stable Rust、`hyperfine`、`awk`が必要です。

測定には、25 kbと50 kbのレコードを交互に含む、2000万塩基・533レコードの決定的なFASTAを使いました。モチーフ検索の条件を揃えるため、生成時にIUPAC曖昧記号を`A`へ変換し、通常のDNA文字だけにしています。計測前に長さフィルター・逆相補・ヘッダー検索・整形・モチーフ位置と鎖向き・FAIを使う領域取得の出力一致を確認し、`len`とSeqKitの`sum_len`も照合しました。各ケースは1回ウォームアップし、5回測定して出力を破棄します。

SeqKitの既定スレッド数は4で、各ツールの既定設定で測定しています。`stats`と`stats --all`の統計項目は一部異なり、`len`と`seqkit stats`は総長だけを比較します。`validate`と`composition`には同じ出力仕様の直接対応コマンドがないため比較していません。SeqKit公式の[使用方法](https://bioinf.shenwei.me/seqkit/usage/)では、対応する`stats`・`seq`・`grep`・`locate`・`faidx`の機能を説明しています。

今回の実装改善では、`stats`と`filter`の塩基ごとの重複チェックを減らし、長さだけで選別する場合はGC/Nを数えないようにしました。逆相補では変換表を使い、`format`は変換なしの入力をまとめてコピーします。同じ環境での改善前後の測定では、`filter`は106.8 msから39.2 ms、`revcomp`は149.7 msから60.9 ms、`format`は54.5 msから28.2 msになりました。`stats`は111.6 msから88.0 msに改善しましたが、今回のSeqKit比較ではまだ遅い結果です。`stats`と`stats --all`の項目は完全一致せず、SeqKitは4スレッドを使います。負荷により測定値は変動します。

| 処理 | `fasta-util` | SeqKit | 平均 ± 標準偏差 |
| --- | ---: | ---: | ---: |
| 総配列長 | `len` | `stats`（`sum_len`のみ） | 15.3 ± 1.6 ms / 31.1 ± 2.5 ms |
| 統計 | `stats` | `stats --all` | 88.0 ± 6.0 ms / 73.6 ± 5.4 ms |
| 長さフィルター（最小40 kb） | `filter` | `seq --min-len 40000` | 39.2 ± 3.3 ms / 37.2 ± 1.1 ms |
| 逆相補 | `revcomp` | `seq --reverse --complement` | 60.9 ± 1.1 ms / 102.7 ± 4.7 ms |
| ヘッダー検索 | `grep` | `grep --by-name --use-regexp` | 19.7 ± 1.0 ms / 48.1 ± 0.5 ms |
| モチーフ位置（`ACGA`） | `locate` | `locate` | 113.7 ± 3.0 ms / 143.4 ± 3.2 ms |
| 幅80に整形 | `format` | `seq --line-width 80` | 28.2 ± 2.3 ms / 40.7 ± 4.9 ms |
| FAIを使う領域取得（10 kb） | `get` | `faidx` | 2.3 ± 0.3 ms / 22.2 ± 1.6 ms |
| インデックス作成 | `index` | `faidx --update-faidx` | 16.1 ± 1.5 ms / 50.6 ± 2.4 ms |

`fasta-util`のFAI取得は5 ms未満で、hyperfineはシェル起動時間の校正精度に警告を出しています。測定値はこの環境における参考値です。

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
- `./scripts/benchmark_comparison.sh [BASES] [RUNS]`はSeqKitと機能が重なる処理を照合して測定します。
- `./scripts/benchmark_real_data.sh [FASTA]`は全体の`len`と`seqkit stats`を比較し、第1染色体のFAI有無による`get`を測ります。`seqkit`が必要で、`seqret`は任意です。

CLIスクリプトの既定値は`benchmark_commands.sh`と`benchmark_comparison.sh`が2000万塩基、その他が100万塩基で、計測回数は5回です。`cargo bench`のサンプル設定はそれぞれの`-- --help`を参照してください。
