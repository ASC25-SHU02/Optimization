import pandas as pd

def load_and_filter(file_path):
    """加载 TSV 文件并筛选出 passed=true 且 pval<=1e-6 的行"""
    df = pd.read_csv(file_path, sep='\t')
    filtered_df = df[(df['passed'] == True) & (df['pval'] <= 1e-6)]
    return filtered_df

def find_intersection(file1, file2, file3, output_file):
    """找到三个文件中满足条件的交集并输出到文件"""
    # 加载并筛选数据
    df1 = load_and_filter(file1)
    df2 = load_and_filter(file2)
    df3 = load_and_filter(file3)

    # 找到交集（基于 ref 和 pos 列）
    intersection = pd.merge(df1, df2, on=['ref', 'pos', 'strand'], how='inner')
    intersection = pd.merge(intersection, df3, on=['ref', 'pos', 'strand'], how='inner')

    # 计算三个 ur 列的平均值
    intersection['ur_avg'] = intersection[['ur_x', 'ur_y', 'ur']].mean(axis=1)

    intersection = intersection[['ref', 'pos', 'strand', 'ur_avg']]

    # 保存结果到输出文件
    intersection.to_csv(output_file, sep='\t', index=False)
    print(f"交集已保存到 {output_file}，共 {len(intersection)} 行。")

# 使用示例
file1 = '../../process/SRR23538290/SRR23538290.filtered.tsv'  # 第一个文件路径
file2 = '../../process/SRR23538291/SRR23538291.filtered.tsv'  # 第二个文件路径
file3 = '../../process/SRR23538292/SRR23538292.filtered.tsv'  # 第三个文件路径
output_file = 'detected.tsv'  # 输出文件路径

find_intersection(file1, file2, file3, output_file)
