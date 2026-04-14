use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingStats {
    // 执行信息
    pub command: String,
    pub input_file: String,
    pub output_file: String,
    pub codebook_file: String,
    pub pattern: String,
    pub separator_config: String,
    pub reverse_complement_enabled: bool,
    pub max_edit_distance: usize,
    pub optimized_mode: bool,
    pub threads: usize,
    pub batch_size: usize,

    // 处理统计
    pub total_sequences: usize,
    pub regex_matched: usize,
    pub regex_failed: usize,
    pub group1_matched: usize,
    pub group2_matched: usize,
    pub final_output: usize,
    pub decoded_sequences: usize, // 成功解码的总数 (Unique + Ambiguous)
    pub ambiguous_sequences: usize, // 存在歧义的 Read 数量
    pub low_quality_sequences: usize, // 错误过多的 Read 数量 (G1+G2 dist > 3)
    pub error_sep_sequences: usize, // Separator 错误的 Read 数量
    pub ambiguous_matched: usize, // 模棱两可的匹配事件数量 (输出到 ambiguous.csv 的行数)
    pub ambiguous_ratio: f64,     // ambiguous_sequences / decoded_sequences
    pub umi_filtered: usize,      // 因 UMI 不足 8bp 而被过滤的 Read 数量
    pub umi_with_n: usize,        // 输出的 UMI 中包含 N 的数量
    pub umi_padded: usize,        // 因 UMI 长度为 8bp 而补 N 的数量

    // 距离分布统计
    pub group1_dist_counts: HashMap<usize, usize>,
    pub group2_dist_counts: HashMap<usize, usize>,

    // 性能统计 (Profiling)
    pub perf_total_levenshtein_calcs: usize,
    pub perf_heuristic_skips: usize,
    pub perf_sep_search_ns: u64,
    pub perf_match_search_ns: u64,

    // 性能统计
    pub processing_time_seconds: f64,
    pub sequences_per_second: f64,
}

impl ProcessingStats {
    pub fn new(
        command: String,
        input_file: String,
        output_file: String,
        codebook_file: String,
        pattern: String,
        separator_config: String,
        reverse_complement_enabled: bool,
        max_edit_distance: usize,
        optimized_mode: bool,
        threads: usize,
        batch_size: usize,
    ) -> Self {
        Self {
            command,
            input_file,
            output_file,
            codebook_file,
            pattern,
            separator_config,
            reverse_complement_enabled,
            max_edit_distance,
            optimized_mode,
            threads,
            batch_size,
            total_sequences: 0,
            regex_matched: 0,
            regex_failed: 0,
            group1_matched: 0,
            group2_matched: 0,
            final_output: 0,
            decoded_sequences: 0,
            ambiguous_sequences: 0,
            low_quality_sequences: 0,
            error_sep_sequences: 0,
            ambiguous_matched: 0,
            ambiguous_ratio: 0.0,
            umi_filtered: 0,
            umi_with_n: 0,
            umi_padded: 0,
            group1_dist_counts: HashMap::new(),
            group2_dist_counts: HashMap::new(),
            perf_total_levenshtein_calcs: 0,
            perf_heuristic_skips: 0,
            perf_sep_search_ns: 0,
            perf_match_search_ns: 0,
            processing_time_seconds: 0.0,
            sequences_per_second: 0.0,
        }
    }

    pub fn merge(&mut self, other: &ProcessingStats) {
        self.total_sequences += other.total_sequences;
        self.regex_matched += other.regex_matched;
        self.regex_failed += other.regex_failed;
        self.group1_matched += other.group1_matched;
        self.group2_matched += other.group2_matched;
        self.final_output += other.final_output;
        self.decoded_sequences += other.decoded_sequences;
        self.ambiguous_sequences += other.ambiguous_sequences;
        self.low_quality_sequences += other.low_quality_sequences;
        self.error_sep_sequences += other.error_sep_sequences;
        self.ambiguous_matched += other.ambiguous_matched;
        self.umi_filtered += other.umi_filtered;
        self.umi_with_n += other.umi_with_n;
        self.perf_total_levenshtein_calcs += other.perf_total_levenshtein_calcs;
        self.perf_heuristic_skips += other.perf_heuristic_skips;
        self.perf_sep_search_ns += other.perf_sep_search_ns;
        self.perf_match_search_ns += other.perf_match_search_ns;

        for (k, v) in &other.group1_dist_counts {
            *self.group1_dist_counts.entry(*k).or_insert(0) += v;
        }
        for (k, v) in &other.group2_dist_counts {
            *self.group2_dist_counts.entry(*k).or_insert(0) += v;
        }
    }
}
