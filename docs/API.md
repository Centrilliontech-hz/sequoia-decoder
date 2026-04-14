# Sequoia Decoder API 文档

本文档描述了 `sequoia-decoder` crate 的内部模块和关键结构体。

## 模块: `optimized_processor`

核心处理模块，负责所有的解码逻辑。

### `struct OptimizedProcessor`
管理解码状态和不可变资源（如编码库、索引）。

#### 方法
*   `new(...)`: 初始化处理器，构建删除邻域索引 (Deletion Neighborhood Index)。
    *   **输入**: 编码库 (Codebook), 分隔符配置, 反向互补设置, 最大编辑距离。
    *   **输出**: `OptimizedProcessor` 实例。
*   `process_paired_fastq_streaming(...)`: 处理文件的主要入口点。
    *   **输入**: R1 路径, R2 路径, 输出路径, 配置参数。
    *   **行为**: 启动生产者/消费者线程池进行并行处理。
*   `find_best_match_dual_window(...)`: 处理变长序列的高级匹配函数。
    *   **输入**: 原始序列切片, 模式标志 (Mode flags)。
    *   **输出**: 最佳 `MatchResult`。
    *   **逻辑**: 同时尝试 16bp 和 15bp 窗口，通过仲裁逻辑解决垃圾碱基干扰。
*   `find_best_match_optimized(...)`: 使用索引+Levenshtein的底层匹配函数。
    *   **输入**: 查询字符串。
    *   **输出**: `MatchResult` (最佳编码, 距离, 候选列表)。

### `struct MatchResult`
*   `best_code`: `Option<String>` - 最佳匹配的编码 ID。
*   `distance`: `usize` - 编辑距离。
*   `candidates`: `Vec<String>` - 具有相同最佳距离的所有候选编码列表（用于歧义检测）。

## 模块: `shape_qc`

处理图像生成。

### `struct ShapeQC`
*   `new(codebook: &[(String, String)])`: 初始化 QC 模块，建立编码库映射。
*   `generate_images(output_base: &str)`: 读取统计数据（或内部状态）并生成 PNG 热图。

### `fn gaussian_blur_f32`
*   **输入**: 扁平化的 `Vec<f32>`, 宽度, 高度, Sigma。
*   **输出**: 模糊处理后的 `Vec<f32>`。
*   **实现**: 使用 `rayon` 进行并行化计算，加速图像生成。

## 模块: `dna_utils`

### `fn build_regex_pattern`
构建用于提取 R1 结构的正则表达式。
*   **输入**: 分隔符配置 (如 "DD,CCC,DD"), 模糊标志, UMI 位置标志。
*   **输出**: 正则表达式字符串。

### `fn rev_comp`
计算 DNA 序列的反向互补序列。
*   **输入**: `&str` (如 "ATCG")。
*   **输出**: `String` (如 "CGAT")。

## 模块: `stats`

### `struct ProcessingStats`
保存处理过程中的计数器，用于生成摘要报告。
*   `regex_matched`, `regex_failed`: 正则匹配成功/失败数。
*   `group1_matched`, `group2_matched`: 各组解码成功数。
*   `ambiguous_matched`: 歧义匹配数。
*   `final_output`: 最终输出序列数。
*   `sequences_per_second`: 处理速度 (序列/秒)。

## 配置 (CLI)
*   `-i, --input`: 输入的 R1 和 R2 FASTQ 文件路径 (Gzipped, 必需, 例如 `-i R1.fq.gz R2.fq.gz`)。
*   `-c, --codebook`: 编码库 CSV 文件路径 (可选, 默认使用内置编码库)。
*   `-o, --output-dir`: 输出目录 (建议使用此选项，会自动管理所有输出文件名)。
*   `-T, --output-table`: 指定输出 CSV 表格路径 (覆盖默认)。
*   `--output-r2`: 指定输出 R2 FASTQ 路径 (覆盖默认)。
*   `--json-summary`: 指定输出 JSON 统计路径 (覆盖默认)。
*   `-j, --threads`: 工作线程数 (默认: CPU 核心数)。
*   `-s, --sep`: 分隔符模式 (默认: "HH,GGG,HH")。
*   `--rc`: 反向互补模式 (0=关闭, 1=开启, -1=自动, 默认: -1)。
*   `--max-distance`: 最大编辑距离 (默认: 3)。
*   `--shape-qc`: 生成 Shape QC 热图 (默认: false)。
*   `--readable-debug-example`: 生成可读的调试示例文件 (ambiguous.example.txt, fail.example.txt)。
*   `--dump-fail`: 输出解码失败的原始序列到 FASTQ 文件。
*   `--ambiguous-mode`: 歧义输出模式 (all/unique/ambiguous, 默认: all)。
