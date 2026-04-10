# Sequoia Decoder 技术设计文档

## 1. 概述
Sequoia Decoder 是一个高性能、多线程的工具，旨在从双端 FASTQ 文件中解码 DNA 序列。它结合了锚点搜索 (Anchor Search) 和模糊字符串匹配（Levenshtein 距离）来识别序列中的特定“邮编”（探针）。该系统针对速度和准确性进行了优化，提供详细的质量控制 (QC) 指标和灵活的输出过滤功能。

## 2. 系统架构

解码器采用基于 `crossbeam-channel` 库的 **生产者-消费者** 架构，以实现高效的并行处理。

### 流水线阶段
1.  **读取线程 (生产者)**:
    *   高效读取双端 FASTQ 文件 (R1/R2)。
    *   将序列打包（默认批次大小：50,000）以减少通道开销。
    *   将批次发送到工作队列。

2.  **工作线程 (消费者)**:
    *   从工作队列中获取批次。
    *   执行核心解码逻辑（锚点搜索、编辑距离计算、Warning 分类）。
    *   更新线程局部统计信息以避免锁竞争。
    *   将处理结果（解码数据、统计信息、QC 更新）发送给写入线程。

3.  **写入线程**:
    *   汇总所有工作线程的结果。
    *   将解码序列写入 `decoded_zipcodes.csv.gz`。
    *   将模棱两可的匹配写入 `ambiguous.csv`（如果启用）。
    *   更新全局 ShapeQC 热图。
    *   将线程局部统计信息合并到全局统计信息中。

## 3. 解码逻辑

每个 Read 对的解码过程涉及以下几个步骤。可以通过 `--strategy` 参数选择解码策略：`greedy` (默认) 或 `min-total-edit`。

### 3.1 锚点搜索 (Anchor Search)
工具**不使用**正则表达式进行初始分割，而是直接搜索预定义的分隔符 (Separator) 序列作为锚点。
*   **目的**: 定位 Separator 的位置，从而推断 Code1, Code2 和 UMI 的边界。
*   **机制**:
    *   **区域搜索**: 首先在预期的位置窗口内（例如 23-30bp）搜索 Separator。
    *   **严格模式**: 默认情况下，如果未启用 `fuzzy-sep`，Separator 必须完全匹配（零错误）。
    *   **模糊模式**: 如果启用 `fuzzy-sep`，允许 Separator 存在一定的编辑距离偏差。
    *   **回退 (Fallback)**: 如果在预期窗口内未找到锚点，且未启用严格模式，则尝试在全序列范围内搜索锚点。

### 3.2 UMI 提取逻辑 (UMI Extraction)
UMI (Unique Molecular Identifier) 的提取基于相对于 **Code 1** (Sequence 1) 的位置，而非固定的绝对位置。

*   **前提**: UMI 始终位于 Code 1 的外侧（Away from Separator）。
*   **提取规则**:
    1.  确定 Code 1 的起始位置 `G1_Start`（通过 Separator 位置和 Code 1 长度反推）。
    2.  **空间检查**: 检查 `G1_Start` 之前（或之后，视方向而定）是否有足够的空间容纳 UMI。
        *   **不足 8bp**: 空间不足，视为解码失败 (UMI Filtered)。
        *   **8bp**: 提取这 8bp，并在远离 Code 1 的一侧补一个 `N`，凑成 9bp。
        *   **>= 9bp**: 提取紧邻 Code 1 的 9bp。
*   **方向性**:
    *   **Standard Mode (非 RC)**: UMI 在 Code 1 左侧。提取 `Seq[G1_Start-9 .. G1_Start]`。
    *   **RC Mode**: 经过反向互补后，逻辑与 Standard Mode 一致（UMI 在 Code 1 左侧）。
*   **统计**: 统计因 UMI 空间不足被过滤的 Read 数 (`UMI Filtered`) 以及输出中包含 `N` 的 UMI 数。

### 3.3 解码策略

#### 3.2.1 贪婪策略 (Greedy Strategy - 默认)
*   **逻辑**: 优先寻找与 Separator 匹配度最高的锚点。
*   **流程**:
    1.  在搜索范围内找到所有潜在的 Separator 候选者。
    2.  按 Separator 质量（错误数越少越好）排序。
    3.  依次尝试每个候选者进行解码 (Match X/Y)。
    4.  一旦找到一个能成功解码 X 和 Y 的组合（编辑距离在允许范围内），即停止并返回结果。
*   **特点**: 偏向于 Separator 的正确性，可能错过 Separator 稍差但 X/Y 匹配更好的组合。

#### 3.2.2 最小总编辑距离策略 (Min Total Edit Strategy)
*   **逻辑**: 搜索链上最小的总编辑距离 (Total Edit Distance)。
*   **公式**: `Total Edit = Dist(X) + Dist(Sep) + Dist(Y)`。
*   **流程**:
    1.  在指定的 Anchor 范围，搜索所有潜在 Separator。
    2.  对每个候选者进行解码，计算 Total Edit。
    3.  对所有结果按 Total Edit 进行排序，找到最小的组合。
    4.  **歧义解决**: 如果存在多个最小 Edit 组合，且 Code 的极差 (Max - Min) 不超过 4，则取中位 Value 输出；否则视为解码失败。
*   **特点**: 综合考虑整体链的质量，能容忍 Separator 的错误以换取更好的 Payload 匹配。

### 3.4 模糊匹配 (Levenshtein 距离)
根据锚点位置提取出的 Code1 和 Code2 候选序列，将与预加载的 **Codebook**（编码本）进行匹配。
*   **算法**: Levenshtein Distance (编辑距离)。
*   **最大距离**: 固定为 3。
*   **机制**:
    *   计算候选序列与 Codebook 中所有条目之间的编辑距离。
    *   识别“最佳匹配”（距离最小的条目）。

### 3.5 Warning 分类
探针根据匹配质量进行分类。单个探针可能同时触发多个 Warning。

1.  **Ambiguous (歧义/模棱两可)**:
    *   **定义**: 
        *   **Greedy**: 对于给定的 Code，多个 Codebook 条目共享*相同*的最佳编辑距离。
        *   **Min Total Edit**: 存在多个最小 Total Edit 的组合，且无法通过中位逻辑解决。
    *   **规则**: 仅当 `Total Edit <= 3` 时计入 Ambiguous。
    *   **处理**: 解码器记录所有候选者。探针被标记为 Ambiguous。

2.  **Low Quality (低质量)**:
    *   **定义**: 总编辑距离过高。
    *   **规则**: `Total Edit > 3` (包含 X, Y 和 Separator 的编辑距离)。
    *   **含义**: 序列与预期 Code 偏差较大，提示测序错误或低质量 Read。

3.  **Separator Error (分隔符错误)**:
    *   **定义**: 识别为 Separator 的序列片段与预期配置不匹配（存在编辑距离）。
    *   **作用**: 主要作为 QC 指标 (Quality Control)。在 `Good` 和 `Merge` 输出模式下，只要 Total Edit 达标，允许存在 Sep Error。

### 3.6 分类逻辑 (Good vs. Bad)
*   **Good 探针**: `Total Edit <= 3` 且无歧义 (Not Ambiguous)。(注意：在 Good 模式下，Sep Error 不作为过滤条件，只要 Total Edit 满足即可)
*   **Bad/Warning 探针**: 上述条件不满足。

### 3.7 R2 ID 处理与输出格式
在输出 R2 FASTQ (`decoded_r2.fastq.gz`) 时，ID 会经过特殊处理以保留关键信息并去除冗余：
1.  **解析**: 读取原始 R1/R2 的 Header，以 `/` 为界截断（舍弃 `/1`, `/2` 及其后内容），同时提取空格后的注释（Annot）。
2.  **重组**: 
    *   新 ID 格式: `@{RealID} {OurInfo} {Annot}`
    *   `OurInfo`: 包含解码元数据的字符串，格式为 `[Code1:Dist1:Start1:Len1:Code2:Dist2:Start2:Len2:UMI]`。
3.  **示例**:
    *   Input: `@SEQ001/1 1:N:0:1`
    *   Output: `@SEQ001 [101:0:0:15:202:0:20:15:NNNNNNNNN] 1:N:0:1`

## 4. 输出模式

解码器通过 `--output-mode` 参数支持四种输出过滤模式：

1.  **All (模式 0)**:
    *   输出所有处理过的序列，无论质量或 Warning 如何。

2.  **Good (模式 1)**:
    *   **仅**输出被分类为 "Good" 的序列。
    *   规则: `!is_ambiguous && !is_low_quality`。
    *   过滤掉 Ambiguous 和 Low Quality 的 Reads。

3.  **Warning (模式 2)**:
    *   **仅**输出至少有一个 Warning 的序列。
    *   用于调试和分析有问题的 Reads。

4.  **Merge (模式 3, 默认)**:
    *   输出 Good 或 Ambiguous 的序列。
    *   规则: `!is_low_quality` (即 `Total Edit <= 3`)。
    *   允许输出存在歧义但质量尚可的 Reads。

## 5. 质量控制 (ShapeQC)

工具生成 "ShapeQC" 报告，在虚拟芯片/流通池上可视化探针质量的空间分布。

### 5.1 矩阵结构
*   **网格**: 一个 2D 矩阵（例如 2048x2048），代表从 Code1 和 Code2 派生的物理或逻辑坐标。
*   **映射**: `Code1` -> X 坐标, `Code2` -> Y 坐标。

### 5.2 每个坐标的指标
对于每个 (X, Y) 位置，系统跟踪：
*   **Good Count**: "Good" 探针（干净匹配）的数量。
*   **Ambiguous Count**: 存在歧义的探针数量。
*   **Low Quality Count**: 高编辑距离的探针数量。
*   **Error Sep Count**: 存在分隔符错误的探针数量。

*注意：由于一个探针可能存在多个 Warning，这些计数的总和可能会超过该坐标的总 Read 数。*

## 6. 统计与监控 (`stats.json`)

`stats.json` 文件提供解码运行的综合摘要。

### 关键指标
*   `total_sequences`: 处理的总输入 Read 对数。
*   `decoded_sequences`: 找到*某种*匹配的总序列数（即使质量很差）。
*   `final_output`: 写入输出文件的序列数（取决于输出模式）。

### Warning 指标
*   `ambiguous_sequences`: 标记为 Ambiguous 的 Read 数量。
*   `low_quality_sequences`: 标记为 Low Quality (`Total Edit > 3`) 的 Read 数量。
*   `error_sep_sequences`: 标记为 Separator Error 的 Read 数量。
*   `ambiguous_matched`: 模棱两可事件的总数（如果一个 Read 中的 Code1 和 Code2 都模棱两可，则此值加 2）。
*   `umi_filtered`: 因 UMI 不足 8bp 而被过滤的 Read 数量。
*   `umi_with_n`: 输出的 UMI 中包含 N 的数量 (包含原生 N 和补的 N)。
*   `umi_padded`: 因 UMI 长度为 8bp 而补 N 的数量。

### 距离分布
*   `group1_dist_counts`: Code1 的编辑距离分布 (0, 1, 2, 3...)。
*   `group2_dist_counts`: Code2 的编辑距离分布。

## 7. 文件格式

### 7.1 解码输出 (`decoded_zipcodes.csv.gz`)
包含解码结果的 CSV 文件。
*   **列**:
    1.  `SeqID`: 序列标识符。
    2.  `UMI`: 唯一分子标识符。
    3.  `Code1`: 解码的 Code1 字符串。
    4.  `Score1`: (保留/常量) 通常为 0。
    5.  `Dist1`: Code1 的编辑距离。
    6.  `Code2`: 解码的 Code2 字符串。
    7.  `Score2`: (保留/常量) 通常为 0。
    8.  `Dist2`: Code2 的编辑距离。

### 7.2 Ambiguous 调试输出 (`ambiguous.csv`)
如果检测到歧义则生成。
*   **列**:
    1.  `SeqID`: 序列标识符。
    2.  `MatchGroup`: "match1" 或 "match2"。
    3.  `SeqType`: "seq1" or "seq2"。
    4.  `Distance`: 编辑距离。
    5.  `Candidates`: 以 `|` 分隔的候选 Code 列表。
    6.  `OriginalSequence`: 原始序列片段。
    7.  `StartPos`: Read 中的起始位置。

## 8. 实现细节

*   **语言**: Rust (2021 Edition)。
*   **依赖**:
    *   `clap`: 命令行参数解析。
    *   `regex`: 模式匹配 (注意：核心解码路径已迁移到直接的 Anchor 搜索，不再依赖 Regex 进行分割)。
    *   `flate2`: Gzip 压缩/解压。
    *   `serde`/`serde_json`: 统计信息的 JSON 序列化。
    *   `rayon`: (可选) 数据并行。
    *   `crossbeam`: 并发原语 (Channel)。
*   **性能优化**:
    *   `FxHasher`: 内部 HashMap 使用更快的哈希算法。
    *   线程局部统计累加，以最小化原子操作。
    *   读写缓冲 IO。

---
*文档最后更新: 2026-02-11*
