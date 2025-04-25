#!/bin/bash
#SBATCH --job-name=funannotate_compare
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=medium
#SBATCH --output=funannotate_compare_%j.log
#SBATCH --error=funannotate_compare_%j.log

# 現在のディレクトリ
dir=`pwd`

# funannotate_compareディレクトリ作成
mkdir -p $dir/funannotate_compare/input
mkdir -p $dir/funannotate_compare/output

# .gbkファイルをコピー
cp $dir/*/funannotate/funannotate/annotate_results/*.gbk $dir/funannotate_compare/input/

# Funannotate環境変数を設定
export SINGULARITYENV_FUNANNOTATE_DB=/lustre10/home/umeyama/funannotate_db/
export SINGULARITYENV_GENEMARK_PATH=/lustre10/home/umeyama/gmes_linux_64_4/
export SINGULARITYENV_PATH=/venv/bin:/home/umeyama/.local/bin:/home/umeyama/gmes_linux_64_4:/home/umeyama/eggnog-mapper:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# 既存の環境変数設定の後に追加
export SINGULARITYENV_OMP_NUM_THREADS=8
export SINGULARITYENV_NUMEXPR_MAX_THREADS=8
export SINGULARITYENV_MKL_NUM_THREADS=8

# スクリプトの最初の方に追加
if [ ! -d "/lustre10/home/umeyama/funannotate_db" ]; then
  echo "Error: Database directory does not exist on host system"
  exit 1
fi

cd $dir/funannotate_compare/

singularity exec /lustre10/home/umeyama/funannotate_latest.sif bash -c "echo \$FUNANNOTATE_DB && ls -la \$FUNANNOTATE_DB"

# compareの実行
# compareの実行部分を修正
singularity exec \
  --bind /lustre10:/lustre10 \
  /lustre10/home/umeyama/funannotate_latest.sif funannotate compare \
  -i $dir/funannotate_compare/input/*.gbk \
  -o funannotate_compare_FUSA_output/ \
  --go_fdr 0.01 \
  --cpus 8
