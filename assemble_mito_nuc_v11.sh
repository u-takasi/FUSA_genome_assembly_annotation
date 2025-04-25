#!/bin/bash
#SBATCH --job-name=genome_analysis
#SBATCH --output=assemble_annotation_%A_%a.log
#SBATCH --error=assemble_annotation_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=320G
#SBATCH --array=1-22
#SBATCH -p medium

set -e
trap 'echo "Error on line $LINENO"' ERR

# Singularityコマンドを関数として定義
filtlong() {
    singularity exec /usr/local/biotools/f/filtlong:0.2.1--hdcf5f25_4 filtlong "$@"
}
flye() {
    singularity exec /usr/local/biotools/f/flye:2.9.5--py39hdf45acc_1 flye "$@"
}
get_organelle() {
    singularity exec /usr/local/biotools/g/getorganelle:1.7.7.1--pyhdfd78af_0 get_organelle_from_assembly.py "$@"
}
minimap2() {
    singularity exec /usr/local/biotools/m/minimap2:2.28--he4a0461_3 minimap2 "$@"
}
racon() {
    singularity exec /usr/local/biotools/r/racon:1.5.0--hdcf5f25_5 racon "$@"
}
samtools() {
    singularity exec /usr/local/biotools/s/samtools:1.9--h91753b0_8 samtools "$@"
}

lrge() {
    singularity exec /usr/local/biotools/l/lrge\:0.1.3--h9f13da3_1 lrge "$@"
}    

# TSVファイルのパス
tsv_file="sample_data.tsv"

# カレントディレクトリを取得
dir=$(pwd)

# TSVファイルをチェック
if [ ! -f "$tsv_file" ]; then
    echo "Error: TSV file $tsv_file not found!"
    exit 1
fi

echo "------------------------------------------------"
echo "Job Information:"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Working directory: $dir"
echo "TSV file: $tsv_file"
echo "------------------------------------------------"

# ヘッダー行をスキップして、SLURM_ARRAY_TASK_IDに応じたデータ行を取得
# +1 するのはヘッダー行をスキップするため
row_num=$((SLURM_ARRAY_TASK_ID + 1))
line=$(sed -n "${row_num}p" "$tsv_file" | sed 's/\r//g' | sed 's/[[:space:]]*$//')

echo "Processing row $row_num from TSV file"
echo "Raw data line: $line"

# 正確に6つのフィールドに分割
# format: directory barcode name species species2 isolate 
read -r directory barcode name species species2 isolate <<< $(echo "$line" | awk '{print $1, $2, $3, $4, $5, $6}')

# 種名を結合
species="$species $species2"

echo "Parsed data:"
echo "Directory: $directory"
echo "Barcode: $barcode"
echo "Name: $name"
echo "Species: $species"
echo "Isolate: $isolate"
echo "------------------------------------------------"

# 入力ファイルパスの設定
# TSVのディレクトリ情報はパスとして扱わず、ベースディレクトリに固定のinput_fastqパスを使用
input_dir="$dir/input_fastq"
echo "Input directory: $input_dir"

output_name="$name"
output_dir="$dir/$output_name"

echo "Output name: $output_name"
echo "Output directory: $output_dir"

# 出力ディレクトリの作成
#mkdir -p "$output_dir"

# 入力ファイルのパスを構築 - 複数の可能なパターンをチェック
fastq_patterns=(
    "$input_dir/$barcode.fastq"
    "$input_dir/$barcode.fq"
    "$input_dir/${barcode}.fastq.gz"
    "$input_dir/${barcode}.fq.gz"
    "$input_dir/$name.fastq"
    "$input_dir/$name.fq"
    "$input_dir/${name}.fastq.gz"
    "$input_dir/${name}.fq.gz"
    "$input_dir/barcode*.fastq"
    "$input_dir/${barcode}*.fastq"
    "$input_dir/*${name}*.fastq"
    "$input_dir/*${barcode}*.fastq"
)

echo "Looking for input FASTQ files:"
for pattern in "${fastq_patterns[@]}"; do
    echo "  Checking $pattern"
done

# 入力ファイルを検索
input_file=""
for pattern in "${fastq_patterns[@]}"; do
    if [ -f "$pattern" ]; then
        input_file="$pattern"
        echo "Found input file: $input_file"
        break
    fi
done

# 入力ファイルが見つからない場合、ディレクトリを調査して終了
if [ -z "$input_file" ]; then
    echo "Error: No input file found!"
    exit 1
fi


# 動的に生成した変数を使用
dir=`pwd`
input_dir="$dir/input_fastq"
output_name="$name"
output_dir="$dir/$output_name"

# 出力用ディレクトリの作成
mkdir -p $output_dir

# barcodeをnameに変換
echo "Processing directory: $input_dir"
echo "Output name: $output_name"

# フィルタリングされたFASTQの出力パス
filtered_fastq="$output_dir/${output_name}_filtered.fastq.gz"

# FiltLongでフィルタリング
filtlong --min_length 1000 --keep_percent 95 "$input_file" | gzip > "$filtered_fastq"
echo "Filtering completed. Output: $filtered_fastq"

# Flyeでアセンブリ
flye --nano-hq "$filtered_fastq" -o $output_dir/flye_$output_name/ -t 12

# ミトコンドリアゲノムの特定
get_organelle -F fungus_mt -g $output_dir/flye_$output_name/assembly_graph.gfa -o $output_dir/mtDNA_$output_name --expected-max-size 100000


# 出力ディレクトリ内のすべてのファスタファイルを取得
output_files=($output_dir/mtDNA_$output_name/fungus_mt.*.path_sequence.fasta)

largest_size=0
organelle_fasta=""

if [[ ${#output_files[@]} -gt 0 ]]; then
    # 各ファイルのサイズをチェック
    for file in "${output_files[@]}"; do
        if [[ -f "$file" ]]; then
            # ファイルサイズを取得（バイト単位）
            current_size=$(stat -c %s "$file")
            
            # より大きいサイズを見つけた場合、更新
            if (( current_size > largest_size )); then
                largest_size=$current_size
                organelle_fasta="$file"
            fi
        fi
    done
    
    if [[ -n "$organelle_fasta" ]]; then
        echo "Selected largest file: $organelle_fasta ($(numfmt --to=iec --suffix=B $largest_size))"
    else
        echo "Error: No valid files found."
        echo "Contents of $output_dir/mtDNA_$output_name:"
        ls -l $output_dir/mtDNA_$output_name
        exit 1
    fi
else
    echo "Error: No output files found."
    echo "Contents of $output_dir/mtDNA_$output_name:"
    ls -l $output_dir/mtDNA_$output_name
    exit 1
fi

# デバッグ出力
echo "Using organelle fasta: $organelle_fasta"

# minimap2でのリードマッピングに使用
minimap2 -t 12 -ax map-ont "$organelle_fasta" "$filtered_fastq" | \
  samtools sort -@ 12 -o $output_dir/mapped_$output_name.bam
samtools view -b -F  4 -@ 12 $output_dir/mapped_$output_name.bam > $output_dir/sorted_mapped_$output_name.bam
samtools fastq -@ 12 $output_dir/sorted_mapped_$output_name.bam > $output_dir/mt_mapped_$output_name.fastq

# Mapped Readsの処理
reads=$output_dir/mt_mapped_$output_name.fastq
threads=12

genome_size=$(lrge -Q 1200 -t "$threads" "$reads")
#genome_size=$($HOME/autocycler-helper-scripts/genome_size_raven.sh "$reads" "$threads")
#genome_size="100000"
echo "LRGE estimated genome size :  $genome_size"

# Step 1: subsample the long-read set into multiple files
# subsampled_reads が存在する場合に削除
if [ -d "$output_dir/subsampled_reads" ]; then
  echo "Removing existing directory: $output_dir/subsampled_reads"
  rm -rf "$output_dir/subsampled_reads"
fi
autocycler subsample --reads "$reads" --out_dir $output_dir/subsampled_reads --genome_size "$genome_size" --min_read_depth 20

# Step 2: assemble each subsampled file
# assemblies が存在する場合に削除                                                                                                                                                                                                                                                                                                                                  
if [ -d "$output_dir/assemblies" ]; then
  echo "Removing existing directory: $output_dir/assemblies"
  rm -rf "$output_dir/assemblies"
fi
mkdir -p $output_dir/assemblies

#除外する　necat nextdenovo raven miniasm

for assembler in canu flye ; do
    for i in 01 02 03 04; do
        $HOME/autocycler-helper-scripts/"$assembler".sh $output_dir/subsampled_reads/sample_"$i".fastq $output_dir/assemblies/"$assembler"_"$i" "$threads" "$genome_size"
    done
done

# Optional step: remove the subsampled reads to save space
#rm $output_dir/subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
autocycler compress -i $output_dir/assemblies -a $output_dir/autocycler_out

# Step 4: cluster the input contigs into putative genomic sequences
autocycler cluster -a $output_dir/autocycler_out --max_contigs 200

# Steps 5 and 6: trim and resolve each QC-pass cluster
for c in $output_dir/autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c"
    autocycler resolve -c "$c"
done

# Step 7: combine resolved clusters into a final assembly
autocycler combine -a $output_dir/autocycler_out -i $output_dir/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa


# autocyclerの出力ファイルを指定
assembly_fasta="$output_dir/autocycler_out/consensus_assembly.fasta"

# ファイルが存在するか確認
if [[ ! -f "$assembly_fasta" ]]; then
    echo "Error: Consensus assembly file not found at $assembly_fasta"
    exit 1
fi


# ポリッシング処理 (autocycler output)
for i in {1..2}; do
    echo "Running minimap2 and racon (iteration $i)..."
    minimap2 -t 12 -x map-ont \
        "$assembly_fasta" \
        "$filtered_fastq" > \
        "$output_dir/minimap2_$i_$output_name.paf" 2>> "$output_dir/minimap2_$i_$output_name.log"

    if [[ ! -s $output_dir/minimap2_$i_$output_name.paf ]]; then
        echo "Error: minimap2 PAF file not generated in iteration $i."
        exit 1
    fi

    racon -t 12 \
        "$filtered_fastq" \
        "$output_dir/minimap2_$i_$output_name.paf" \
        "$assembly_fasta" > \
        "$output_dir/racon_${i}_$output_name.fasta" 2>> "$output_dir/racon_${i}_$output_name.log"

    if [[ ! -s $output_dir/racon_${i}_$output_name.fasta ]]; then
        echo "Error: racon output not generated in iteration $i."
        exit 1
    fi

    # raconの結果を次の入力として指定
    mv "$output_dir/racon_${i}_$output_name.fasta" "$assembly_fasta"
done

# Medakaによるポリッシング
echo "Running medaka consensus..."
medaka_mapped_dir="$output_dir/mapped_$output_name"
mkdir -p "$medaka_mapped_dir"

singularity exec /usr/local/biotools/m/medaka:2.0.1--py39hf77f13f_0 bash -c "
    export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH &&
    medaka_consensus \
        -i $output_dir/mapped_$output_name.bam \
        -d $output_dir/autocycler_out/consensus_assembly.fasta \
        -o $medaka_mapped_dir \
        -m r1041_e82_400bps_sup_v5.0.0 \
        -t 12
"

if [[ ! -d "$medaka_mapped_dir" || -z "$(ls -A $medaka_mapped_dir)" ]]; then
    echo "Error: Medaka output not generated for mapped reads."
    exit 1
fi

# Unmapped Readsの処理
echo "Processing unmapped reads..."
samtools view -b -f 4 "$output_dir/mapped_$output_name.bam" > \
    "$output_dir/unmapped_$output_name.bam"

if [[ ! -s $output_dir/unmapped_$output_name.bam ]]; then
    echo "Error: Unmapped BAM file not generated."
    exit 1
fi

samtools fastq -@ 12 "$output_dir/unmapped_$output_name.bam" > "$output_dir/unmapped_$output_name.fastq"

if [[ ! -s $output_dir/unmapped_$output_name.fastq ]]; then
    echo "Error: Unmapped FASTQ file not generated."
    exit 1
fi

# Flyeでアセンブリ
echo "Running Flye on unmapped reads..."
flye --nano-hq "$output_dir/unmapped_$output_name.fastq" \
    -o "$output_dir/flye_unmapped_$output_name/" \
    -t 12

if [[ $? -ne 0 ]]; then
    echo "Error: Flye assembly on unmapped reads failed."
    exit 1
fi

if [[ ! -s "$output_dir/flye_unmapped_$output_name/assembly.fasta" ]]; then
    echo "Error: Flye output assembly not generated."
    exit 1
fi

# ポリッシング処理 (unmapped reads)
for i in {1..2}; do
    echo "Running minimap2 and racon on unmapped reads (iteration $i)..."
    minimap2 -t 12 -x map-ont \
        "$output_dir/flye_unmapped_$output_name/assembly.fasta" \
        "$output_dir/unmapped_$output_name.fastq" > \
        "$output_dir/minimap2_unmapped_$i_$output_name.paf" 2>> "$output_dir/minimap2_unmapped_$i_$output_name.log"

    if [[ ! -s $output_dir/minimap2_unmapped_$i_$output_name.paf ]]; then
        echo "Error: minimap2 PAF file on unmapped reads not generated in iteration $i."
        exit 1
    fi

    racon -t 12 \
        "$output_dir/unmapped_$output_name.fastq" \
        "$output_dir/minimap2_unmapped_$i_$output_name.paf" \
        "$output_dir/flye_unmapped_$output_name/assembly.fasta" > \
        "$output_dir/racon_unmapped_${i}_$output_name.fasta" 2>> "$output_dir/racon_unmapped_${i}_$output_name.log"

    if [[ ! -s $output_dir/racon_unmapped_${i}_$output_name.fasta ]]; then
        echo "Error: racon output on unmapped reads not generated in iteration $i."
        exit 1
    fi

    mv "$output_dir/racon_unmapped_${i}_$output_name.fasta" \
        "$output_dir/flye_unmapped_$output_name/assembly.fasta"
done

# Medakaによるポリッシング (unmapped reads)
echo "Running medaka consensus on unmapped reads..."
unmapped_medaka_dir="$output_dir/unmapped_medaka_$output_name"
mkdir -p "$unmapped_medaka_dir"

singularity exec /usr/local/biotools/m/medaka:2.0.1--py39hf77f13f_0 bash -c "
    export LD_LIBRARY_PATH=/usr/local/lib:\$LD_LIBRARY_PATH &&
    medaka_consensus \
        -i $output_dir/unmapped_$output_name.fastq \
        -d $output_dir/flye_unmapped_$output_name/assembly.fasta \
        -o $unmapped_medaka_dir \
        -m r1041_e82_400bps_sup_v5.0.0 \
        -t 12
"

if [[ ! -d "$unmapped_medaka_dir" || -z "$(ls -A $unmapped_medaka_dir)" ]]; then
    echo "Error: Medaka output not generated for unmapped reads."
    exit 1
fi

# `longstitch` 用のディレクトリ作成                                                                                                                                                                                                                                                                                                                                       
longstitch_dir="$output_dir/longstitch"
mkdir -p $longstitch_dir

# 必要なファイルをコピー                                                                                                                                                                                                                                                                                                                                                  
cp "$output_dir/unmapped_medaka_$output_name/consensus.fasta" "$longstitch_dir/"
cp "$output_dir/unmapped_$output_name.fastq" "$longstitch_dir/"

# `longstitch` ディレクトリに移動                                                                                                                                                                                                                                                                                                                                         
cd $longstitch_dir

# consensus.fasta を consensus.fa にリネーム                                                                                                                                                                                                                                                                                                                              
mv consensus.fasta consensus.fa

# consensus.fa の総塩基数を計算                                                                                                                                                                                                                                                                                                                                           
total_length=$(awk '!/^>/ {total += length($0)} END {print total}' consensus.fa)
echo "Total length of consensus.fa (bp): $total_length"

# Mb（メガベース）単位に変換し、小数点以下を切り捨てて整数にする                                                                                                                                                                                                                                                                                                          
G_value=$(( total_length / 1000000 ))
echo "Computed G value (in Mb): $G_value"

# `longstitch` の実行                                                                                                                                                                                                                                                                                                                                                     
echo "Running longstitch with G=$G_value in $longstitch_dir"
singularity exec /usr/local/biotools/l/longstitch\:1.0.5--hdfd78af_0 \
    longstitch tigmint-ntLink-arks draft=consensus reads=unmapped_$output_name G=$G_value

# 成功メッセージ                                                                                                                                                                                                                                                                                                                                                          
echo "Longstitch analysis completed successfully for $output_name"

cd $dir


# ディレクトリ構造の設定
dir=`pwd`
input_dir="$dir/input_fastq"
output_name="$name"
output_dir="$dir/$output_name"
threads=12
echo "Directory structure prepared: $output_dir"

# MITOSによるアノテーション
mkdir -p "$output_dir/mitos"
echo "Starting MITOS annotation..."
singularity exec /usr/local/biotools/m/mitos:2.1.9--pyhdfd78af_0 runmitos.py \
  -i "$output_dir/mapped_$output_name/consensus.fasta" \
  -o "$output_dir/mitos" \
  --refdir "$HOME" \
  --refseqver refseq89f \
  -c 04
echo "MITOS annotation completed."


output_mitos_dir="$output_dir/mitos/"
bed_file="${output_mitos_dir}/result.bed"
fasta_file="$output_dir/mapped_$output_name/consensus.fasta"
output_file="${output_mitos_dir}/final_mtDNA.fasta"


# 一時ファイル作成
filtered_bed=$(mktemp)

# result.bedからrrnLかつscore=0.0の行を抽出
awk '$4 ~ /rrnL/ && $5 == "0.0"' "$bed_file" > "$filtered_bed"

# 抽出結果の確認
if [ ! -s "$filtered_bed" ]; then
    echo "No rrnL entries with score 0.0 found in BED file."
    rm "$filtered_bed"
    exit 1
fi

# rrnL行を1行ずつ処理
while read -r line; do
    # BED情報を取得
    seq_id=$(echo "$line" | cut -f1)
    start=$(echo "$line" | cut -f2)
    end=$(echo "$line" | cut -f3)
    strand=$(echo "$line" | cut -f6)

    if [ "$strand" == "-" ]; then
        # -方向の場合の処理
        restart_pos=$((end + 1))
        seqkit restart -i "$restart_pos" -o restarted.fasta "$fasta_file"
        seqkit seq -t dna -r -p -o "$output_file" restarted.fasta
        rm restarted.fasta
        echo "Processed strand '-' with restart position $restart_pos"

    elif [ "$strand" == "+" ]; then
        # +方向の場合の処理
        restart_pos=$((start + 1))
        seqkit restart -i "$restart_pos" -o "$output_file" "$fasta_file"
        echo "Processed strand '+' with restart position $restart_pos"
    else
        echo "Unknown strand direction in line: $line"
        continue
    fi
done < "$filtered_bed"

# ヘッダー名を変更
sed -i "1s/^>.*/>${name}_mtDNA/" "$output_file"

# 一時ファイルの削除
rm "$filtered_bed"

# 実行完了メッセージ
if [ -f "$output_file" ]; then
    echo "Restarted sequence saved to $output_file"
else
    echo "Failed to process rrnL sequence."
fi


# MITOSによるアノテーション
mkdir -p "$output_dir/mitos2"
echo "Starting MITOS annotation..."
cp $output_file "$output_dir/mitos2"
singularity exec /usr/local/biotools/m/mitos:2.1.9--pyhdfd78af_0 runmitos.py \
  -i "$output_dir/mitos2/final_mtDNA.fasta" \
  -o "$output_dir/mitos2" \
  --refdir "$HOME" \
  --refseqver refseq89f \
  -c 04
echo "2nd MITOS annotation completed."

# Singularityコマンドを関数として定義
funannotate() {
    singularity exec ~/funannotate_latest.sif funannotate "$@"
}

# TSVファイルのパス
tsv_file="sample_data.tsv"

# SLURM_ARRAY_TASK_IDに応じてTSVからデータを取得
read line < <(awk -F'\t' -v line=$SLURM_ARRAY_TASK_ID 'NR==line+1 {print $0}' $tsv_file | sed 's/\r//g' | sed 's/[[:space:]]*$//')

# TSVから値を変数に割り当て
IFS=$'\t' read -r directory barcode name species isolate <<< "$line"
species=${species// /_}

# デバッグ出力
echo "Raw line: $line"
echo "Directory: $directory"
echo "Barcode: $barcode"
echo "Name: $name"
echo "Species: $species"
echo "Isolate: $isolate"

# 動的に生成した変数を使用
dir=`pwd`
output_name="$name"
output_dir="$dir/$output_name"
thread=12

# 出力ディレクトリ作成
mkdir -p $dir/$name/funannotate

# Funannotate環境変数を設定
export SINGULARITYENV_FUNANNOTATE_DB=/home/umeyama/funannotate_db/
export SINGULARITYENV_GENEMARK_PATH=/home/umeyama/gmes_linux_64_4/
export SINGULARITYENV_PATH=/venv/bin:/home/umeyama/.local/bin:/home/umeyama/gmes_linux_64_4:/home/umeyama/eggnog-mapper:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# ファイルの名前変換と初期処理
funannotate clean \
    -i $dir/$name/longstitch/consensus.k32.w100.tigmint-ntLink.longstitch-scaffolds.fa \
    -o $dir/$name/funannotate/$name.clean.fasta
echo "Step: funannotate clean completed."

funannotate sort \
    -i $dir/$name/funannotate/$name.clean.fasta \
    -o $dir/$name/funannotate/$name.sort.fasta \
    -b contig \
    --minlen 500
echo "Step: funannotate sort completed."

funannotate mask \
    -i $dir/$name/funannotate/$name.sort.fasta \
    -o $dir/$name/funannotate/$name.softmasked.fasta
echo "Step: funannotate mask completed."

# predict
funannotate predict \
    -i $dir/$name/funannotate/$name.softmasked.fasta \
    -o $dir/$name/funannotate/funannotate \
    -s "$species" \
    --isolate "$isolate" \
    --cpus $thread \
    --optimize_augustus \
    --min_training_models 180 \
    --name $name
echo "Step: funannotate predict completed."


gbk=$dir/$name/funannotate/funannotate/predict_results/${species}_$isolate.gbk
protein=$dir/$name/funannotate/funannotate/predict_results/${species}_$isolate.proteins.fa


# interproscan
mkdir -p $dir/$name/funannotate/iprscan
cd $dir/$name/funannotate/iprscan
~/my_interproscan/interproscan-5.72-103.0/interproscan.sh \
    -i $protein \
    -o $dir/$name/funannotate/iprscan/iprscan.xml \
    -goterms \
    -f XML \
    -T $dir/$name/funannotate/iprscan/tmp/ \
    -cpu $thread \
    -dp
echo "Step: funannotate iprscan completed."

# antismash
mkdir -p $dir/$name/funannotate/antismash
cd $dir/$name/funannotate/antismash
export SINGULARITYENV_LANG=C.UTF-8
singularity exec /home/umeyama/standalone_latest.sif antismash \
	    --output-dir $dir/$name/funannotate/antismash/ \
	    --taxon fungi \
	    --genefinding-tool glimmerhmm \
	    --cpus $thread \
	    $gbk
echo "Step: antismash completed."


# eggnog
export EGGNOG_DATA_DIR=/home/umeyama/eggnog-mapper-data
mkdir -p $dir/$name/funannotate/eggnog
/home/umeyama/eggnog-mapper/emapper.py \
    -i $protein \
    --output_dir $dir/$name/funannotate/eggnog \
    --output $name \
    --cpu $thread
echo "Step: eggnog completed."


# phobius
mkdir -p $dir/$name/funannotate/phobius
cd $dir/$name/funannotate/phobius
perl ~/phobius/phobius.pl -short $protein > $dir/$name/funannotate/phobius/phobius.results.txt
echo "Step: phobius completed."


# SignalP6の環境を有効化
cd ~/signalp6_slow_sequential/
source signalp_env/bin/activate

# signalp6
mkdir -p $dir/$name/funannotate/signalp6
cd $dir/$name/funannotate/signalp6
signalp6 \
    --fastafile $protein \
    --output_dir $dir/$name/funannotate/signalp6 \
    --organism eukarya \
    --mode slow-sequential
echo "Step: signalp6 completed."

# オプション：環境を非アクティブ化（必要な場合）
deactivate

#iprscan.xmlのコピー
mkdir -p $dir/$name/funannotate/funannotate/annotate_misc/
cp $dir/$name/funannotate/iprscan/iprscan.xml $dir/$name/funannotate/funannotate/annotate_misc/iprscan.xml

# annotate
cd $dir/$name/funannotate
funannotate annotate \
	    -i $dir/$name/funannotate/funannotate/predict_results/ \
	    --isolate "$isolate" \
	    --cpus $thread \
	    --eggnog $dir/$name/funannotate/eggnog/$name.emapper.annotations \
	    --iprscan $dir/$name/funannotate/iprscan/iprscan.xml \
	    --antismash $dir/$name/funannotate/antismash/${species}_$isolate.gbk \
	    --phobius $dir/$name/funannotate/phobius/phobius.results.txt \
	    --signalp $dir/$name/funannotate/signalp6/prediction_results.txt \
	    --sbt ~/template.sbt \
	    --force
echo "Step: funannotate annotate completed."
