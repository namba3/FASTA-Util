# FASTA-Util

[日本語](README.md) | [English](README.en.md)

FASTA ファイルを扱うための CLI ツールです。

## FASTA形式について

FASTAは、塩基配列やアミノ酸配列をテキストで表す形式です。各レコードは`>`で始まるヘッダー行から始まり、その次の行以降に配列を記述します。配列は複数行に折り返して記述できます。

`len`と`get`は、既定で核酸配列を処理します。タンパク質配列を処理する場合は、`--sequence-type protein`を指定してください。`revcomp`と`locate`は核酸配列を処理し、`filter`、`stats`、`composition`は配列種別を自動判定します。

```fasta
>record-1 optional description
ACGTN
UKS-
>record-2
MRY
```

このツールの`len`はヘッダーと空行を除いて配列記号を数えます。`get`の数値範囲は1始まり・両端を含み、複数レコードの配列をファイル順に連結した位置で範囲を指定します。出力には範囲の終端までに現れたヘッダーを残すため、選択範囲の配列記号がないレコードのヘッダーも含まれることがあります。

核酸配列では、大文字・小文字の`ACGTNUKSYMWRBDHV`とギャップを表す`-`を受け付けます。タンパク質配列では、大文字・小文字の標準20アミノ酸記号に加え、`B J O U X Z`、終止記号`*`、ギャップ`-`を受け付けます。`get`は元の大文字・小文字を維持し、許可されていない記号はエラーになります。

入力ファイルを使って`len`または`get`を実行している間は、そのファイルを変更しないでください。`.fai`を使った`get`の実行中も、FASTA入力とインデックスの両方を変更しないでください。

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

| コマンド | 機能 |
| --- | --- |
| [`len`](#len) | 配列記号の総数を数える |
| [`validate`](#validate) | FASTA構造と配列記号を検証する |
| [`index`](#index) | ランダムアクセス用の`.fai`を作成する |
| [`stats`](#stats) | 長さ・N50・GC率などを集計する |
| [`composition`](#composition) | 配列記号ごとの割合を集計する |
| [`get`](#get) | レコードまたは指定範囲を取得する |
| [`filter`](#filter) | 長さ・GC率・N率でレコードを選ぶ |
| [`revcomp`](#revcomp) | 核酸配列を逆相補鎖に変換する |
| [`grep`](#grep) | ヘッダー文字列でレコードを検索する |
| [`locate`](#locate) | モチーフ配列の位置を検索する |
| [`format`](#format) | 改行・文字ケース・ギャップなどを整形する |

配列を入力に取るコマンドは、出力先を指定しなければ結果を標準出力に書きます。`stats`・`composition`・`filter`・`revcomp`・`format`・`validate`は入力ファイルを省略すると標準入力を読みます。`get`・`grep`・`locate`では入力位置に`-`を指定してください。`len`は既定で標準入力を読み、ファイル入力には`-i`を使います。`index`は`.fai`の保存先を決めるためファイル入力が必要です。

```sh
fasta-util filter --min-len 1000 input.fa |
  fasta-util revcomp |
  fasta-util stats
```

### len

配列の総記号数を数えます。

```sh
./target/release/fasta-util len -i test.fna
./target/release/fasta-util len --sequence-type protein -i proteins.faa
```

```txt
10000
```

### stats

FASTA全体のレコード数、合計・最小・最大・平均長、N50、GC率、`N`率、配列種別を表示します。`--each`を付けると、レコードごとのID・長さ・GC率・`N`率を表示します。

既定の`auto`では、核酸配列にないタンパク質記号があればタンパク質として判定します。それ以外は`T`と`U`からDNA/RNAを推定し、どちらもない場合は判別不能と表示します。`--sequence-type nucleotide`または`--sequence-type protein`で明示指定できます。タンパク質のGC率と`N`率はテキストでは`n/a`、JSONでは`null`です。

```sh
./target/release/fasta-util stats genome.fa
./target/release/fasta-util stats genome.fa --each
./target/release/fasta-util stats genome.fa --format json
./target/release/fasta-util stats proteins.faa --sequence-type protein --format json
```

### composition

配列全体に含まれる記号ごとの割合を表示します。核酸配列では`A`・`C`・`G`・`T`（RNAでは`U`）・`N`と、入力に含まれるIUPAC曖昧記号やギャップを表示し、`GC`には`G`と`C`の合計割合を表示します。タンパク質では20種類の標準アミノ酸と、入力に含まれる拡張記号の割合を表示します。割合の分母はギャップや曖昧記号を含む配列記号の総数です。配列種別は`stats`と同様に自動判定し、必要に応じて`--sequence-type nucleotide`または`--sequence-type protein`を指定できます。

```sh
./target/release/fasta-util composition seq.fa
./target/release/fasta-util composition proteins.faa --sequence-type protein
cat seq.fa | ./target/release/fasta-util composition
```

### filter

配列ごとの長さ、GC割合、`N`割合でレコードを選びます。割合は`0.0`から`1.0`で指定し、指定した条件はすべて満たす必要があります。GCと`N`の割合は、ギャップや曖昧な記号を含む配列全体の長さを分母にします。配列種別は自動判定し、タンパク質だけにある記号が含まれればタンパク質として扱います。核酸配列と区別できないタンパク質配列は`--sequence-type protein`で明示できます。GCと`N`の条件は核酸配列でのみ使えます。選択したレコードは入力順で標準出力に出力し、`-o`/`--output`でファイルにも保存できます。

```sh
./target/release/fasta-util filter seq.fa --min-len 1000
./target/release/fasta-util filter seq.fa --max-len 10000
./target/release/fasta-util filter seq.fa --min-gc 0.40 --max-gc 0.60
./target/release/fasta-util filter seq.fa --max-n 0.05
./target/release/fasta-util filter proteins.fa --min-len 100 --sequence-type protein
```

### revcomp

各レコードの核酸配列を逆相補鎖に変換します。DNAでは`A`を`T`、RNAでは`A`を`U`に対応させます。IUPAC曖昧塩基とギャップ`-`にも対応し、元の大文字・小文字を維持します。同じレコード内に`T`と`U`が混在する配列はエラーになります。出力は既定で60塩基ごとに折り返し、幅は`--chars-per-line`で変更できます。`-o`/`--output`を指定するとファイルに保存します。

```sh
./target/release/fasta-util revcomp seq.fa
./target/release/fasta-util revcomp seq.fa --chars-per-line 80
./target/release/fasta-util revcomp seq.fa --output seq.revcomp.fa
```

### grep

FASTAヘッダーに含まれる文字列でレコードを検索し、一致したレコード全体を入力順に出力します。既定では大文字・小文字を区別します。`--ignore-case`で区別をなくし、`--invert-match`で一致しないレコードを選べます。

```sh
./target/release/fasta-util grep seq.fa BRCA
./target/release/fasta-util grep seq.fa brca --ignore-case
```

### format

FASTAの改行をLFに統一し、配列行を指定幅で折り返します。既定の幅は60です。`--width 0`は各レコードの配列を1行にまとめます。`--uppercase`または`--lowercase`で配列記号の大文字・小文字を変更でき、`--remove-gaps`で`-`を除去します。`--trim-header`は`>`の後のヘッダー文字列の前後にあるASCII空白を除去します。出力先を指定しない場合は標準出力に書きます。

```sh
./target/release/fasta-util format --width 80 seq.fa > formatted.fa
./target/release/fasta-util format --width 0 seq.fa
./target/release/fasta-util format --uppercase --remove-gaps --trim-header seq.fa
```

### locate

核酸配列内のモチーフを両鎖から検索し、タブ区切りで`ID・開始・終了・strand`を出力します。座標は1始まり・両端包含です。検索語はIUPAC塩基記号に対応し、重複する位置も報告します。`--max-mismatch`で許容する不一致数を指定できます。逆相補配列と一致する位置は`-`、入力配列の向きで一致する位置は`+`です。曖昧塩基を含む配列では、各位置の塩基集合がモチーフの集合と重なる場合に一致とみなします。

```sh
./target/release/fasta-util locate genome.fa AATAAA
./target/release/fasta-util locate seq.fa ATGNNNTAA
./target/release/fasta-util locate seq.fa AATAAA --max-mismatch 1
```

### get

指定したIDのレコード全体、`ID:開始-終了`形式のレコード内領域、または`開始-終了`形式の全体範囲を取得します。座標は1始まりで両端を含みます。全体範囲は全レコードの配列をファイル順に連結した位置です。複数IDや`--ids`ファイルに対応し、出力はFASTA内の出現順です。入力の隣に`.fai`があれば自動で使い、なければ先頭から検索します。`--fai-index`でインデックスを明示でき、`--no-fai-index`で自動利用を無効にできます。配列の折り返し幅は`--chars-per-line`で指定できます。

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

FASTAのレコード構造、配列記号、レコードIDの重複、レコード内の改行形式、`.fai`で扱える行幅かを確認します。成功時はレコード数と配列種別を表示し、失敗時はファイル名・行・列と該当行を含む診断を表示して、終了コード1を返します。レコードごとにLFとCRLFが異なるファイルは受け付けます。

既定の核酸モードでは、`T`を含む配列をDNA、`U`を含む配列をRNAと表示します。`T`と`U`の両方がある入力はエラーです。どちらも含まない配列はDNA/RNAを判別できないため、その旨を表示します。タンパク質配列は`--sequence-type protein`で指定します。ギャップ`-`は許可されます。

```sh
./target/release/fasta-util validate seq.fa
./target/release/fasta-util validate proteins.faa --sequence-type protein
```

### index

FASTAを1回ストリーミングして、ランダムアクセス用の`.fai`インデックスを作成します。出力先は入力ファイル名に`.fai`を追加したパスです。入力形式に問題がある場合は作成を中止し、既存のインデックスを保持します。

```sh
./target/release/fasta-util index genome.fa
./target/release/fasta-util get genome.fa 100000001-100001000 --fai-index genome.fa.fai
```

## ベンチマーク

全サブコマンド、標準入力とパイプ、SeqKit・seqretとの機能比較、実ゲノム、Rust内の処理方式を測定する手順と結果を[`docs/benchmark-results.ja.md`](docs/benchmark-results.ja.md)にまとめています。CLI測定は`hyperfine`・`awk`・stable Rustを使い、既定で2000万塩基、各ケース1回ウォームアップ・5回計測です。

```sh
./scripts/benchmark_commands.sh
./scripts/benchmark_commands.sh 50000000 7
./scripts/benchmark_comparison.sh
cargo bench --bench nucleic_acid
cargo bench --bench fasta_io
./scripts/benchmark_real_data.sh
```

ベンチマーク結果はCPU、OS、ファイルキャッシュや実行負荷によって変わります。小さい処理ではプロセス起動時間が結果に含まれます。

## ライセンス

このプロジェクトは [MIT License](LICENSE-MIT) または [Apache License 2.0](LICENSE-APACHE) のいずれかの条件で利用できます。MIT License の著作権者表記は `namba3 (GitHub: @namba3)` です。
