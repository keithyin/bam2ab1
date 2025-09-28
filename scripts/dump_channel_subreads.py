import pysam
import sys
from tqdm import tqdm


def extract_records_by_ch(bam_file, output_bam, ch_value):
    """
    从BAM文件中提取CH标签等于指定值的reads并保存到新BAM文件
    :param bam_file: 输入BAM文件路径
    :param output_bam: 输出BAM文件路径
    :param ch_value: 要匹配的CH标签值
    """
    # 打开输入BAM文件
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40) as input_bam:
        # 创建输出BAM文件
        with pysam.AlignmentFile(output_bam, "wb", template=input_bam, check_sq=False, threads=40) as output_bam:
            # 遍历输入BAM文件中的每个read
            for read in tqdm(input_bam, desc=f"reading {bam_file}"):
                # 检查read是否包含CH标签且值匹配
                if read.has_tag('ch') and read.get_tag('ch') == ch_value:
                    # 将匹配的read写入输出BAM文件
                    output_bam.write(read)


if __name__ == "__main__":

    extract_records_by_ch("/data-slow/md0-backup/data/ccs_data/data2025Q1/S_aureus_2h/20250207_Sync_Y0003_08_H01_Run0001_called.adapter.bam",
                          "test-data/channel_340564_subreads.bam", 340564)
    # print(f"Extracted records with CH={ch_value} to {output_bam}")
