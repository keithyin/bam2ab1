import pysam


def find_first_high_rq_read(bam_file):
    """
    从BAM文件中查找并打印第一个RQ≥0.999的序列
    :param bam_file: BAM文件路径
    """
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for read in bam:
            # 获取RQ标签值，不同测序平台可能使用不同标签名
            rq = read.get_tag('rq') if read.has_tag('rq') else None

            # 如果找到RQ≥0.999的read则打印并退出
            if rq is not None and rq >= 0.999:
                print(f"Found read with rq={rq}:")
                print(f"channel:{read.get_tag("ch")}")
                print(read.to_string())
                print(f"Sequence: {read.query_sequence}")
                return

    print("No read found with RQ≥0.999")


if __name__ == "__main__":

    find_first_high_rq_read(
        "/data-slow/md0-backup/data/ccs_data/data2025Q1/S_aureus_2h/20250207_Sync_Y0003_08_H01_Run0001_called.smicing.smc_all_reads.bam")
