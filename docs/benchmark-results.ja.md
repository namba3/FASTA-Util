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

計測前に、`stats`・`composition`・`filter`・`revcomp`・`format`・`validate`のファイル入力と標準入力で出力が一致することを確認しました。`len`・`stats`・`composition`は1 worker、2 worker、4 worker、自動設定で出力一致も確認しました。`get`ではFAIの有無による出力一致も確認し、核酸・タンパク質の検証と統計コマンドを実行しました。対象範囲は後半レコードの`record_000201:1001-11000`です。`get`のFAIは計測前に作成し、独立した`index`ケースでは毎回FAIを書き込みます。

以下はhyperfineによる平均 ± 標準偏差です。FAIを使う`get`はプロセス起動時間に近いため、小さな時間差は参考値として扱ってください。

| コマンド | 平均 ± 標準偏差 |
| --- | ---: |
| `len`（核酸） | 12.9 ± 1.4 ms |
| `validate`（核酸） | 107.3 ± 4.0 ms |
| `validate`（タンパク質） | 29.8 ± 2.4 ms |
| `index`（400レコード） | 13.8 ± 0.2 ms |
| `stats`（核酸） | 51.1 ± 3.5 ms |
| `stats`（タンパク質） | 41.3 ± 3.6 ms |
| `composition`（核酸） | 24.3 ± 3.5 ms |
| `composition`（タンパク質） | 21.8 ± 0.6 ms |
| `get`（FAIなし、10 kb） | 11.4 ± 0.6 ms |
| `get`（FAIあり、10 kb） | 1.9 ± 0.2 ms |
| `filter`（核酸） | 42.1 ± 4.0 ms |
| `filter`（タンパク質） | 40.4 ± 1.4 ms |
| `revcomp` | 62.4 ± 4.2 ms |
| `grep` | 14.5 ± 0.9 ms |
| `locate`（4記号モチーフ、不一致1を許容） | 649.8 ± 15.4 ms |
| `format`（幅60） | 21.0 ± 2.1 ms |

| 標準入力・パイプ | 平均 ± 標準偏差 |
| --- | ---: |
| `stats`（stdin） | 183.8 ± 4.0 ms |
| `composition`（stdin） | 91.6 ± 1.2 ms |
| `filter`（stdin） | 59.4 ± 4.2 ms |
| `revcomp`（stdin） | 78.9 ± 7.8 ms |
| `format`（stdin） | 39.4 ± 5.3 ms |
| `validate`（stdin） | 167.1 ± 8.6 ms |
| `filter | revcomp | stats` | 275.0 ± 10.3 ms |

この測定では`grep`のばらつきが他の多くのケースより大きくなりました。FAIを使う`get`は5 ms未満で、hyperfineからシェル起動時間の校正精度について警告が出ています。

## 並列worker数の比較

`benchmark_commands.sh`はファイル入力の`len`・`stats`・`composition`について、1・2・4 workerと自動設定を同じ入力で測ります。自動設定は8 MiB未満なら逐次処理し、それ以上では利用可能なCPU数に応じて最大4 workerを使います。`stats`はrecord境界で分割するため、recordが400件あるこの入力では複数workerを使えます。

以下は上記と同じ2000万塩基・400レコードの核酸FASTAを使い、各条件を1回ウォームアップ後に5回測定した平均 ± 標準偏差です。出力は破棄しています。

| コマンド | 1 worker | 2 workers | 4 workers | 自動設定 |
| --- | ---: | ---: | ---: | ---: |
| `len` | 19.6 ± 2.4 ms | 14.8 ± 1.5 ms | 12.2 ± 0.9 ms | 12.4 ± 1.1 ms |
| `stats` | 157.1 ± 3.5 ms | 81.7 ± 3.2 ms | 49.6 ± 4.3 ms | 52.3 ± 5.2 ms |
| `composition` | 50.5 ± 2.4 ms | 34.3 ± 5.0 ms | 23.9 ± 1.1 ms | 21.1 ± 0.8 ms |

この入力と環境では、1 workerに対し4 workerで`len`は約1.6倍、`stats`は約3.2倍、`composition`は約2.1倍でした。自動設定も4 workerに近い結果です。これはページキャッシュが温まった状態の同一環境での値で、CPU負荷やFASTAのrecord数・長さの偏りで変わります。特に`stats`はrecord単位で分割するため、単一recordの巨大FASTAでは並列化されません。

`BASES`は少なくとも10,011,000を指定してください。測定対象の領域が201番目のレコード内に必要です。

## RefSeq GRCh38.p14データセット

`dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna`を再測定に使用しました。ファイルは3,339,739,109 bytes、705レコードで、配列記号の合計は3,298,430,636です。第1染色体レコード`NC_000001.11`は248,956,422塩基です。小文字のソフトマスク配列は維持しました。隣接FAIがあっても走査する条件には`--no-fai-index`を指定し、FAI経路は`--fai-index`で明示しました。計測前に両経路の出力がバイト単位で一致することを確認しました。`get`は1回ウォームアップ後に5回、全ゲノムの`len`・`stats`・`composition`は1回ウォームアップ後に3回計測しました。全ゲノム3コマンドは1 workerと自動設定の出力も一致しました。

| コマンド | 結果 | 平均 ± 標準偏差 |
| --- | ---: | ---: |
| `seqkit stats`（全ゲノム） | 705レコード、3,298,430,636塩基 | 1.617 ± 0.021 s |
| `fasta-util len`（全ゲノム、自動設定） | 3,298,430,636塩基 | 1.046 ± 0.020 s |

### 全ゲノムでのworker数比較

下表は3.3 GBの同じFASTAを、1 worker・2 worker・4 worker・自動設定で計測した平均 ± 標準偏差です。hyperfineのウォームアップは1回、計測は各3回です。

| コマンド | 1 worker | 2 workers | 4 workers | 自動設定 |
| --- | ---: | ---: | ---: | ---: |
| `len` | 1.504 ± 0.013 s | 1.211 ± 0.050 s | 1.054 ± 0.063 s | 1.046 ± 0.020 s |
| `stats` | 18.067 ± 0.527 s | 10.300 ± 0.116 s | 6.377 ± 0.031 s | 6.512 ± 0.115 s |
| `composition` | 7.089 ± 0.226 s | 4.454 ± 0.116 s | 2.783 ± 0.062 s | 2.811 ± 0.075 s |

4 workerでは1 workerより`len`が約1.4倍、`stats`が約2.8倍、`composition`が約2.5倍速くなりました。`stats`はrecord境界で分割し、このデータセットの705 recordsをworker間で処理します。計測は温まったページキャッシュ上で行いました。`composition`の4 worker測定にはhyperfineが外れ値の警告を出しているため、差の小さい値は参考値として扱ってください。

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

今回の測定では、`fasta-util`の`len`・`stats`・`composition`は既定で自動worker設定を使います。`stats`とSeqKitの`stats --all`は項目が完全一致せず、SeqKitは4スレッドを使います。測定値は負荷により変動します。

| 処理 | `fasta-util` | SeqKit | 平均 ± 標準偏差 |
| --- | ---: | ---: | ---: |
| 総配列長 | `len` | `stats`（`sum_len`のみ） | 12.0 ± 0.9 ms / 30.4 ± 1.9 ms |
| 統計 | `stats` | `stats --all` | 35.0 ± 2.2 ms / 69.0 ± 6.7 ms |
| 長さフィルター（最小40 kb） | `filter` | `seq --min-len 40000` | 41.6 ± 3.9 ms / 36.1 ± 2.1 ms |
| 逆相補 | `revcomp` | `seq --reverse --complement` | 60.7 ± 1.1 ms / 96.7 ± 4.9 ms |
| ヘッダー検索 | `grep` | `grep --by-name --use-regexp` | 16.8 ± 2.2 ms / 43.0 ± 4.3 ms |
| モチーフ位置（`ACGA`） | `locate` | `locate` | 113.3 ± 3.8 ms / 141.9 ± 10.0 ms |
| 幅80に整形 | `format` | `seq --line-width 80` | 20.0 ± 2.1 ms / 33.6 ± 1.6 ms |
| FAIを使う領域取得（10 kb） | `get` | `faidx` | 1.8 ± 0.4 ms / 25.2 ± 1.9 ms |
| インデックス作成 | `index` | `faidx --update-faidx` | 14.7 ± 0.7 ms / 52.1 ± 3.6 ms |

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
