# 更新日志

本项目的所有重要更改都将记录在此文件中。

## [0.4.0] - 2026-02-13

### 新增
- **UMI 提取逻辑**: 实现了基于 Code 1 位置的相对定位 UMI 提取逻辑。增加了空间检查（过滤 < 8bp 的 UMI，为 8bp UMI 补 N）。
- **统计信息**: 在 `stats.json` 中新增了 `umi_filtered`（UMI 过滤数）、`umi_with_n`（含 N 的 UMI 数）和 `umi_padded`（补 N 的 UMI 数）字段。
- **R2 输出格式**: 重构了 R2 序列 ID 以包含解码元数据，格式为 `@{RealID} {OurInfo} {Annot}`。
- **CLI 增强**: 将所有输出路径参数统一为 `--output-dir`。移除了已弃用的 `--output-table`、`--output-r2` 和 `--json-summary` 参数。

### 变更
- **默认配置**: 更新了默认参数：`strategy` 默认为 `min-total-edit`，`output-mode` 默认为 `merge`，`shape-qc` 和 `fuzzy-sep` 默认启用。
- **文档**: 更新了 `README.md` 和 `docs/DESIGN.md` 以反映新的 UMI 提取逻辑和参数变更。

## [0.3.0] - 2026-02-12

### 新增
- **解码策略**: 新增 `min-total-edit` 策略，优先选择总编辑距离 (G1+Sep+G2) 最小的候选结果。
- **输出模式**: 新增 `Merge` 模式，同时输出 Good 和 Ambiguous 的 Reads (要求 Total Edit <= 3)。
- **Anchor 控制**: 新增 `--anchor` 参数用于指定搜索范围，以及 `--strict-anchor` 用于控制是否允许回退搜索。
- **测试框架**: 添加了全面的测试数据生成器 (`tests/scripts/test_generator.py`) 和性能对比脚本。

### 变更
- **歧义解决**: 改进了歧义处理逻辑，在 `min-total-edit` 模式下，如果多个最佳候选的 Code 范围 <= 4，则取中值输出。
- **Codebook**: 增加了对内置 Codebook 的支持，不再强制要求提供外部 CSV 文件。

### 修复
- **匹配选择**: 修复了当多个候选结果编辑距离相同时，未能优先选择最长匹配的问题。

## [0.2.0] - 2026-02-05

### 新增
- **调试输出**: 增加了生成可读调试示例文件 (`ambiguous.example.txt`, `fail.example.txt`) 和导出失败序列 (`fail.R1/R2.fastq.gz`) 的选项。
- **统计信息**: 增强了 `stats.json` 输出，提供更详细的匹配结果分类统计。

### 变更
- **系统架构**: 将并行处理重构为“生产者-消费者-写入器”流水线架构，显著提升了扩展性和吞吐量。

## [0.1.0] - 2026-02-04

### 新增
- **初始版本**: 发布高性能 Rust DNA 序列解码器。
- **双窗口搜索**: 实现了 16bp/15bp 双窗口搜索机制，提高 Anchor 检测的鲁棒性。
- **Shape QC**: 支持生成空间分布热图 (Shape QC) 以进行质量控制。
- **模糊匹配**: 集成了基于 SymSpell 的模糊匹配算法和 SIMD 加速的 Levenshtein 距离计算。
- **并行计算**: 利用 `rayon` 实现高效的数据并行处理。
