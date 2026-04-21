
#!/bin/bash

# 配置目录路径
BAM_DIR="/mnt/e/ab1-debug/20260420-Consensus/ab1"
REF_DIR="/mnt/e/ab1-debug/20260420-Consensus/ab1"

# 创建带时间戳的日志文件
LOG_FILE="bam2ab1_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== 开始批量转换任务: $(date) ==="
echo "BAM目录: $BAM_DIR"
echo "参考序列目录: $REF_DIR"

process_file() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .sort.bam)
    local ref_file="$REF_DIR/${base_name}.consensus.fasta"

    if [[ -f "$ref_file" ]]; then
        echo "处理文件: $base_name"
        echo "BAM文件: $bam_file"
        echo "参考文件: $ref_file"
        
        # --chunkSize 500 --ovlpSize 50
        if bam2ab1 --bam "$bam_file" --ref "$ref_file" --insIdent N ; then
            echo "转换成功: $base_name"
            return 0
        else
            echo "转换失败: $base_name" >&2
            return 1
        fi
    else
        echo "未找到匹配的参考文件: $base_name" >&2
        return 2
    fi
}

# ./target/release/bam2ab1 --bam test-data/channel_340564.sort.bam --ref test-data/channel_340564.consensus.fasta --chunkSize 500 --ovlpSize 50 --insIdent N


# 主处理循环
find "$BAM_DIR" -name "*.sort.bam" | while read -r bam; do
    process_file "$bam"
    echo "----------------------------------------"
done

echo "=== 批量转换任务完成: $(date) ==="
echo "详细日志请查看: $LOG_FILE"
