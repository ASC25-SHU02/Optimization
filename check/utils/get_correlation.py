import pandas as pd

def match_and_extract(file1, file2, output_file):
    """
    匹配两个文件中前三列一致的行，并提取第一个文件的 ur_avg 和第二个文件的 ur。
    """
    # 读取两个 TSV 文件
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t')
    df1.rename(columns={'ratio': 'ur', 'chromosome': 'ref', 'position': 'pos'}, inplace=True)
    df1['ref'] = df1['ref'].astype(str)
    df2['ref'] = df2['ref'].astype(str)
    df1['pos'] = df1['pos'].astype(str)  # 将 pos 列转换为字符串
    df2['pos'] = df2['pos'].astype(str)  # 将 pos 列转换为字符串
    # 合并两个数据框，基于前三列（ref, pos, strand）
    merged_df = pd.merge(df1, df2, on=['ref', 'pos', 'strand'], how='inner')

    # 提取所需的列（ur_avg 和 ur）
    result_df = merged_df[['ur_avg', 'ur']]

    # 保存结果到输出文件
    result_df.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"结果已保存到 {output_file}，共 {len(result_df)} 行。")

# 使用示例
file1 = 'true.tsv'  # 第一个文件路径
file2 = 'detected.tsv'  # 第二个文件路径
output_file = 'correlation.tsv'  # 输出文件路径

match_and_extract(file1, file2, output_file)
