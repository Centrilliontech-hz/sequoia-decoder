# Sequoia Decoder 测试指南

## 1. 单元测试 (Unit Tests)
本项目包含针对核心逻辑的单元测试。

### 运行测试
```bash
cargo test
```

### 关键测试用例
*   **`test_extract_groups_from_sequence`** (`src/dna_utils.rs`):
    *   验证正则表达式的生成和分组提取逻辑。
    *   测试各种分隔符配置。
*   **`test_unified_matching_priority`** (`src/optimized_processor.rs`):
    *   **核心算法验证**。
    *   验证优先级逻辑：精确匹配 (16bp) > 精确匹配 (15bp 前缀) > 1 Mismatch。
    *   验证歧义检测逻辑 (Ambiguity Detection)。
*   **`test_trim_trailing_mismatches`** (`src/optimized_processor.rs`):
    *   验证在生成输出 CSV 时，是否正确截去了末端的错配碱基。

## 2. 集成测试 (自动)
本项目包含完善的 Rust 集成测试套件，位于 `tests/integration_tests.rs`。

### 运行测试
```bash
# 运行所有集成测试
cargo test --test integration_tests -- --nocapture
```

该测试会自动：
1.  构建 `release` 版本的二进制文件。
2.  遍历 `tests/assets` 下的所有测试数据集。
3.  使用 `--head 200000` 参数运行解码器，验证功能和性能。
4.  检查输出文件是否存在且非空。

## 3. 模拟仿真与准确率评估 (Simulation & Accuracy)

为了深入评估解码器在各种复杂场景下的表现，我们开发了基于 Python 的模拟测试框架 `test_generator.py`。该框架能够模拟真实的测序错误，并根据 Ground Truth 自动计算各项性能指标。

### 3.1 核心特性
*   **Mixin 架构**: 测试脚本采用模块化设计，将 DNA 工具、锚点模式 (HD/HDC/TCR) 以及干扰场景 (Scenario) 分离，便于扩展和维护。
*   **全场景模拟**:
    *   `normal`: 理想情况。
    *   `ambiguity`: 模拟编码库中极其相似的序列，测试歧义消解能力。
    *   `interference`: 在关键位置引入随机碱基干扰。
    *   `shift`: 模拟序列整体偏移（如插入/缺失导致的错位）。
    *   `sep_error`: 在分隔符 (Separator) 位置引入突变。
    *   `adaptor_error`: 模拟接头序列污染。
    *   `structural_error`: 破坏序列的整体结构。
*   **高级突变模型**: `mutate_dna` 支持替换 (Sub)、插入 (Ins) 和缺失 (Del)，并能通过参数控制比例，确保所有突变至少产生 1 个编辑距离。

### 3.2 运行评估
```bash
python3 tests/scripts/test_generator.py --mode <MODE> --n <READS_PER_SCENARIO>
```
*   `--mode`: 支持 `HD`, `HDC`, `HDCv3`, `HDC-TCR`, `HDCv3-TCR`。
*   `--n`: 每个场景生成的 Read 数量（推荐 500+）。

### 3.3 统计指标说明
评估报告提供 X/Y 两端独立的统计数据，以便精准定位问题：
*   **Decode %**: 成功被解码器识别并输出的比例。
*   **X-Corr % / Y-Corr %**: 在已解码的 Read 中，X/Y 索引完全正确的比例。
*   **X-Mean / Y-Mean**: 误判索引与真实索引之间的平均偏移距离。
*   **X-Q90 / Y-Q90**: 索引偏移距离的 90 分位数。如果该值较小而 Mean 较大，说明存在极个别严重的误判，而非系统性偏移。
*   **X-Max / Y-Max**: 最大偏移距离。

### 3.4 自动化批处理
使用 `test_all_anchor.py` 可以一键运行所有锚点模式的测试，并生成一份汇总报告 `full_test_report.md`。

## 4. 手动验证命令示例
如果您需要手动运行特定数据集进行调试：

```bash
./target/release/sequoia-decoder \
  -i tests/assets/2026-01-28T05-19-17Z/QGMB0122-1_L2_UDI502.R1.fastq.gz \
     tests/assets/2026-01-28T05-19-17Z/QGMB0122-1_L2_UDI502.R2.fastq.gz \
  -o output/manual_test \
  -j 16 \
  --readable-debug-example
```

## 4. 验证标准
1.  **速度**: 在 32 线程机器上，处理 100 万条 reads 约 30 秒。
2.  **统计**:
    *   正则匹配率 (Regex Match Rate) > 95%
    *   解码率 (Decoding Rate) > 70%
    *   歧义率 (Ambiguity Rate) < 30% (具体取决于文库质量)
3.  **输出**:
    *   `decoded.csv.gz` 应为有效的 GZIP 文件。
    *   `stats.json` 应包含有效的 JSON 数据。
    *   `*.png` 图片文件应被正确生成且可查看。

## 5. 性能分析 (Performance Profiling)
如需检查性能瓶颈，请运行工具并观察标准输出中的 `Detailed Performance Breakdown` 部分：
```text
Detailed Performance Breakdown:
  - Core Processing (Threads): ...
  - Shape QC Generation: ...
  - Output Table Merging: ...
```
*   **Shape QC Generation** 应非常快 (约 10s 或更少)。
*   **Output Table Merging** 应几乎瞬间完成 (< 1s)。
