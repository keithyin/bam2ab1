#!/bin/bash

# 配置根目录（请根据实际情况修改或通过参数传入）
ROOT_DIR="/mnt/e/ab1-debug/20260420-plasmid/Plasmid"

# 创建带时间戳的日志文件
LOG_FILE="bam2ab1_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== 开始批量转换任务: $(date) ==="
echo "根目录: $ROOT_DIR"

process_folder() {
    local folder=$1
    local base_name=$(basename "$folder" _plasmid)
    
    local bam_file="$folder/${base_name}.assembly.bam"
    local ref_file="$folder/${base_name}_plassembler_plasmids.fasta"

    if [[ -f "$bam_file" && -f "$ref_file" ]]; then
        echo "处理文件夹: $base_name"
        echo "BAM文件: $bam_file"
        echo "参考文件: $ref_file"
        
        # 如需添加 --chunkSize 和 --ovlpSize，可在此行追加
        if bam2ab1 --bam "$bam_file" --ref "$ref_file" --insIdent N --chunkSize 1000 --ovlpSize 50; then
            echo "转换成功: $base_name"
            return 0
        else
            echo "转换失败: $base_name" >&2
            return 1
        fi
    else
        echo "文件夹 $base_name 中缺少所需文件（${base_name}.assembly.bam 或 ${base_name}_plassembler_plasmids.fasta)" >&2
        return 2
    fi
}

# 遍历根目录下的一级子文件夹
find "$ROOT_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r folder; do
    process_folder "$folder"
    echo "----------------------------------------"
done

echo "=== 批量转换任务完成: $(date) ==="
echo "详细日志请查看: $LOG_FILE"