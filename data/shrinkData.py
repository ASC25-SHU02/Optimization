import random

def process_file(input_file, output_file):
    # 读取文件并过滤以@开头的行
    with open(input_file, 'r') as file:
        lines = [line for line in file if not line.startswith('@')]

    # 随机选取10000行，保持原有顺序
    if len(lines) > 10000:
        selected_indices = sorted(random.sample(range(len(lines)), 10000))
        selected_lines = [lines[i] for i in selected_indices]
    else:
        selected_lines = lines

    # 将选中的行写入输出文件
    with open(output_file, 'w') as file:
        file.writelines(selected_lines)

# 使用示例
if __name__ == "__main__":
    input_file = 'input.tmp'  # 输入文件名
    output_file = 'Shrinked.in'  # 输出文件名
    process_file(input_file, output_file)
