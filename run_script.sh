
#!/bin/bash

# 配置目录路径
BAM_DIR="/data-slow/qingke-deliver-data/20250818_240601Y0012_Run0003/Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Bam"
REF_DIR="/data-slow/qingke-deliver-data/20250818_240601Y0012_Run0003/Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Sequences"

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
        
        if bam2ab1 --bam "$bam_file" --ref "$ref_file"; then
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



# 主处理循环
find "$BAM_DIR" -name "*.bam" | while read -r bam; do
    process_file "$bam"
    echo "----------------------------------------"
done

echo "=== 批量转换任务完成: $(date) ==="
echo "详细日志请查看: $LOG_FILE"
