# Sequoia Decoder - 高性能 DNA 序列解码器

Sequoia Decoder 是一个高度优化的多线程工具，用于处理和解码双端 DNA 测序数据 (FASTQ)。它专为每秒处理数百万条 Reads 而设计，能够根据提供的编码库提取特定的条形码 ("Zipcodes") 和 UMI。

## 功能特性

*   **高性能**: 使用 Rust 编写，利用多线程、SIMD 指令集和高效的流式 I/O (支持 GZIP)。
*   **双端处理**: 同时处理 R1 (结构信息) 和 R2 (生物学序列)。
*   **便捷易用**: 
    *   **内置编码库**: 默认内置标准 Codebook，无需额外指定。
    *   **智能默认值**: 优化的默认参数配置。
*   **稳健的匹配算法**:
    *   **双窗口搜索 (Dual Window Search)**: 自动测试 16bp 和 15bp 窗口，以处理测序伪影。
    *   **混合索引**: 使用基于哈希的删除邻域索引进行 O(1) 查找，结合 SIMD Levenshtein 验证。
    *   **歧义消解**: 严格执行“全局最小距离”优先级 (精确 > 距离1 > 距离2)。
*   **灵活配置**:
    *   可自定义分隔符模式 (支持模糊匹配)。
    *   自动检测反向互补模式。
*   **调试与质控**:
    *   **Shape QC**: 生成空间热图以可视化解码质量。
    *   **详细日志**: 支持生成人类可读的调试文件 (`ambiguous.example.txt`, `fail.example.txt`)。
    *   **故障转储**: 支持导出解码失败的原始序列以便分析。

## 安装

```bash
cd /pipeline/tools/simple_decoder
cargo build --release
# 二进制文件位于: target/release/sequoia-decoder
```

## 使用方法

### 基本用法

最简单的调用方式只需要指定输入文件（R1 和 R2）：

```bash
./target/release/sequoia-decoder -i input_R1.fastq.gz input_R2.fastq.gz -o output_dir
```

### 完整参数示例

```bash
./target/release/sequoia-decoder \
    -i input_R1.fastq.gz input_R2.fastq.gz \
    -o output_dir \
    -j 16 \
    -s "GGG" \
    --shape-qc 1 \
    --fuzzy-sep 1 \
    --strict-anchor 1 \
    --output-mode merge \
    --strategy min-total-edit
```

### 参数说明

| 参数 | 缩写 | 描述 | 默认值 |
|------|------|------|--------|
| `--input <R1> <R2>` | `-i` | 输入的双端 FASTQ 文件路径 (必需) | - |
| `--output-dir <DIR>` | `-o` | 输出目录 (必需) | - |
| `--codebook <FILE>` | `-c` | 编码库 CSV 文件路径 (可选) | 使用内置库 |
| `--threads <NUM>` | `-j` | 线程数 | CPU 核心数 |
| `--sep <PATTERN>` | `-s` | 分隔符序列 (例如: GGG 或 CCC) | "GGG" |
| `--rc <MODE>` | - | 反向互补 (-1=自动, 0=关, 1=开) | -1 |
| `--fuzzy-sep <0/1>` | - | 允许 Separator 错配 (0=禁用, 1=启用) | 1 |
| `--strict-anchor <0/1>` | - | 严格 Anchor 模式 (0=禁用, 1=启用) | 1 |
| `--shape-qc <0/1>` | - | 生成 QC 热图图片 (0=禁用, 1=启用) | 1 |
| `--output-mode <MODE>` | - | 输出模式 (all, good, warning, merge) | merge |
| `--strategy <MODE>` | - | 解码策略 (greedy, min-total-edit) | min-total-edit |
| `--readable-debug-example <0/1>` | - | 生成可读的调试示例文件 (0=禁用, 1=启用) | 0 |
| `--dump-fail` | - | 输出解码失败的序列到 FASTQ | false |

## 输出文件

在输出目录中，您将看到：

*   **`decoded_zipcodes.csv.gz`**: 主要解码结果表。
    *   格式: `SeqID, UMI, X-Pos, X-Edit, X-Distance, Y-Pos, Y-Edit, Y-Distance`
*   **`decoded_r2.fastq.gz`**: 对应的 R2 序列 (ID 已包含解码信息，详见 DESIGN 文档)。
*   **`stats.json`**: 详细的运行统计信息（包含解码率、歧义率等）。
*   **`shape_qc_*.png`**: (如果启用 `--shape-qc`) 空间分布热图。
*   **`ambiguous.example.txt`**: (如果启用 `--readable-debug-example`) 歧义匹配的可视化示例。
*   **`fail.example.txt`**: (如果启用 `--readable-debug-example`) 解码失败的可视化示例。
*   **`fail.R1.fastq.gz` / `fail.R2.fastq.gz`**: (如果启用 `--dump-fail`) 解码失败的原始序列。

## 运行测试

本项目包含多层级的测试体系，从单元测试到基于模拟数据的全流程验证。

### 1. 模拟数据与准确率评估 (推荐)

使用 `test_generator.py` 工具可以生成带有 Ground Truth 的模拟数据，并自动评估解码器的准确率、鲁棒性以及索引偏移情况。

```bash
# 生成模拟数据并运行评估 (以 HD 模式为例)
python3 tests/scripts/test_generator.py --mode HD --n 500

# 运行所有支持的模式 (HD, HDC, HDCv3, HDC-TCR, HDCv3-TCR)
python3 tests/scripts/test_all_anchor.py
```

评估结果将保存至 `tests/reports/` 目录下，包含各场景下的解码率、X/Y 端独立准确率、平均偏移距离 (Mean)、Q90 偏移等详细指标。

### 2. 集成测试 (自动)

使用 Rust 内置的测试框架运行针对真实资产数据的集成测试：

```bash
cargo test --test integration_tests -- --nocapture
```

### 3. 单元测试

运行核心算法和工具函数的单元测试：

```bash
cargo test
```

更多详细信息请参阅 [测试指南](docs/TESTING.md)。
