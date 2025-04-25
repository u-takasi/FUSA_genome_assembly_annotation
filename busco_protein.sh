#!/bin/bash
#SBATCH --job-name=busco_protein_analysis
#SBATCH --output=busco_protein_%A_%a.log
#SBATCH --error=busco_protein_%A_%a.log
#SBATCH --mem=144G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

set -e
trap 'echo "Error on line $LINENO"; rm -rf $BUSCO_DOWNLOAD_PATH' ERR EXIT

# TSVファイルのパス
tsv_file="sample_data.tsv"

# SLURM_ARRAY_TASK_IDに応じてTSVからデータを取得
read line < <(awk -F'	' -v line=$SLURM_ARRAY_TASK_ID 'NR==line+1 {print $0}' $tsv_file | sed 's/\r//g' | sed 's/[[:space:]]*$//')

# TSVから値を変数に割り当て
IFS=$'	' read -r directory barcode name species isolate <<< "$line"

# デバッグ出力
echo "Raw line: $line"
echo "Directory: $directory"
echo "Barcode: $barcode"
echo "Name: $name"
echo "Species: $species"
echo "Isolate: $isolate"



# 作業ディレクトリ設定
dir=`pwd`
output_name="$name"
output_dir="$dir/$output_name"

# スペースをアンダースコアに置き換えてから、ファイルパスを組み立てる
species=$(echo "$species" | sed 's/ /_/g')

echo "Species: $species"
echo "Isolate: $isolate"
echo "Assembled genome path: $output_dir/funannotate/funannotate/annotate_results/${species}_${isolate}.scaffolds.fa"

# アセンブル済みタンパク質ファイルのパス
protein_file="$output_dir/funannotate/funannotate/annotate_results/${species}_${isolate}.proteins.fa"

# 出力先ディレクトリがない場合は作成
busco_output_dir="$output_dir/busco_protein"
mkdir -p $busco_output_dir  # busco_proteinディレクトリを作成

# 個別のダウンロードディレクトリを `/tmp/` に作成
BUSCO_DOWNLOAD_PATH="/tmp/busco_downloads_protein_$SLURM_ARRAY_TASK_ID"
mkdir -p $BUSCO_DOWNLOAD_PATH

# BUSCOの実行 (タンパク質ファイルの場合は -m proteins を指定)
singularity exec /usr/local/biotools/b/busco\:5.8.2--pyhdfd78af_0 \
    busco -i "$protein_file" \
    -o "busco_output_protein_${species}_${isolate}" \
    --out_path "$busco_output_dir" \
    -m proteins \
    --download_path $BUSCO_DOWNLOAD_PATH \
    --auto-lineage-euk \
    -f --cpu 4

# 成功時のメッセージ
echo "BUSCO protein analysis completed successfully for $species_$isolate"

# 解析後にダウンロードディレクトリを削除
rm -rf $BUSCO_DOWNLOAD_PATH
echo "Temporary BUSCO download directory removed: $BUSCO_DOWNLOAD_PATH"
