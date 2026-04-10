use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use clap::Parser;
use serde::{Deserialize, Serialize};

mod debug_formatter;
mod dna_utils;
mod lookup_table;
mod optimized_processor;
mod shape_qc;
mod stats;

use dna_utils::{build_regex_pattern, should_enable_reverse_complement};
use optimized_processor::*;
use stats::ProcessingStats;
#[derive(clap::ValueEnum, Clone, Debug)]
enum OutputMode {
    All,
    Good,
    Warning,
    Merge,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入 R1 和 R2 FASTQ.GZ 文件路径 (例如: -i R1.fq.gz R2.fq.gz)
    #[arg(short, long, num_args = 2, required = true)]
    input: Vec<String>,

    /// 编码库文件路径 (默认: 使用内置编码库)
    #[arg(short, long)]
    codebook: Option<String>,

    /// 输出目录 (如果设置，将覆盖 output_table, output_r2, json_summary 为默认路径)
    #[arg(short, long)]
    output_dir: String,

    /// 线程数
    #[arg(short = 'j', long, default_value_t = num_cpus::get())]
    threads: usize,

    /// 批处理大小
    #[arg(short, long, default_value_t = 50000)]
    batch_size: usize,

    /// 进度显示间隔（序列数）
    #[arg(long, default_value_t = 10000)]
    progress_interval: usize,

    /// 分隔符序列 (例如: GGG 或 CCC)
    #[arg(short, long, default_value = "GGG")]
    sep: String,

    /// 反向互补设置：0=不启用, 1=启用, -1=自动决定 (默认: -1)
    #[arg(long, default_value_t = -1, allow_hyphen_values = true)]
    rc: i32,

    /// 只处理前 K 个序列对 (可选)
    #[arg(long)]
    head: Option<usize>,

    /// 是否允许 Separator 错配 (0=禁用, 1=启用) (默认: 1)
    #[arg(long, default_value_t = 1)]
    fuzzy_sep: u8,

    /// Anchor 搜索范围设置
    /// 格式支持:
    /// 1. 单个数字: 在此位置 +/- 2bp 范围内搜索 (例如: --anchor 20)
    /// 2. 两个数字: 指定起始和结束范围 [start, end) (例如: --anchor 18 25)
    /// 3. 预设模式名称: HD, HDC, HDCv3, HDC-TCR, HDCv3-TCR
    /// 默认: 根据 sep 自动决定 (包含 GGG -> HD, 包含 CCC -> HDC-TCR)
    #[arg(long, num_args = 1..=2)]
    anchor: Option<Vec<String>>,

    /// 是否生成 Shape QC 图片 (0=禁用, 1=启用) (默认: 1)
    #[arg(long, default_value_t = 1)]
    shape_qc: u8,

    /// 强制覆盖输出文件 (默认: false)
    #[arg(long, default_value_t = false)]
    force: bool,

    /// Output Mode (默认: merge)
    #[arg(long, value_enum, default_value_t = OutputMode::Merge)]
    output_mode: OutputMode,

    /// 生成可读的调试示例文件 (0=禁用, 1=启用) (默认: 0)
    #[arg(long, default_value_t = 0)]
    readable_debug_example: u8,

    /// 输出解码失败的序列 (fail.R1.fastq.gz, fail.R2.fastq.gz)
    #[arg(long, default_value_t = false)]
    dump_fail: bool,

    /// 启用严格的 Anchor 模式 (0=禁用, 1=启用) (默认: 1)
    /// 如果禁用此选项 (设置为 0)，当 Anchor 搜索失败时，将回退到全序列正则表达式搜索。
    /// 默认情况下，强制使用严格模式 (失败直接报错)。
    #[arg(long, default_value_t = 1)]
    strict_anchor: u8,

    /// Decoding Strategy (greedy, min-total-edit)
    #[arg(long, value_enum, default_value_t = DecodingStrategy::MinTotalEdit)]
    strategy: DecodingStrategy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessedSequence {
    pub seqid_r1: String,
    pub seq1: String,
    pub code1: Option<String>,
    pub seq2: String,
    pub code2: Option<String>,
    pub seq3: String,
    // R2 信息
    pub seqid_r2: String,
    pub seq_r2: String,
    pub qual_r2: String,
}

#[derive(Debug, Clone)]
pub struct SequenceMatch {
    pub group1: String,
    pub group2: String,
    pub group3: String,
}

const DEFAULT_CODEBOOK_CSV: &str = include_str!("codebook.csv");

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let input_r1 = &args.input[0];
    let input_r2 = &args.input[1];

    println!("开始处理DNA序列 (双端模式)...");
    println!("输入 R1: {}", input_r1);
    println!("输入 R2: {}", input_r2);

    if args.strict_anchor == 0 && args.anchor.is_none() {
        return Err("--strict-anchor 0 requires --anchor to be specified.".into());
    }

    // 处理输出路径逻辑
    let output_dir = &args.output_dir;
    let path = std::path::Path::new(output_dir);
    if !path.exists() {
        std::fs::create_dir_all(output_dir)?;
    } else if !args.force && path.read_dir()?.next().is_some() {
        // 如果目录存在且不为空，且未指定 --force，则警告 (但在批处理场景下，用户可能希望覆盖)
    }

    let output_table = path.join("decoded_zipcodes.csv.gz")
        .to_string_lossy()
        .to_string();
    let output_r2 = path.join("R2.fastq.gz")
        .to_string_lossy()
        .to_string();
    let json_summary = path.join("stats.json")
        .to_string_lossy()
        .to_string();

    println!("输出 解码表: {}", output_table);
    println!("输出 R2: {}", output_r2);
    println!("输出 Summary: {}", json_summary);
    println!("编码库: {}", args.codebook.as_deref().unwrap_or("Built-in"));
    println!("线程数: {}", args.threads);

    let start_time = Instant::now();

    // 设置线程池
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // 加载编码库
    println!("加载编码库...");
    let codebook = load_codebook(args.codebook.as_deref())?;
    println!("编码库加载完成，共{}个编码", codebook.len());

    // 构建正则表达式模式
    let enable_reverse_complement = should_enable_reverse_complement(&args.sep, args.rc);
    let umi_on_head = !enable_reverse_complement;
    let pattern = build_regex_pattern(&args.sep, args.fuzzy_sep != 0, umi_on_head);

    // 构建执行命令
    let command = std::env::args().collect::<Vec<String>>().join(" ");

    // Hardcode max_distance
    let max_distance = 3;

    // 创建统计信息
    let stats = ProcessingStats::new(
        command,
        format!("{}, {}", input_r1, input_r2),
        output_table.clone(),
        args.codebook
            .clone()
            .unwrap_or_else(|| "Built-in".to_string()),
        pattern,
        args.sep.clone(),
        enable_reverse_complement,
        max_distance,
        true, // 强制开启优化模式
        args.threads,
        args.batch_size,
    );

    // 强制使用优化模式 (OptimizedProcessor)
    println!("使用高性能双端处理模式...");
    let processor = OptimizedProcessor::new(
        codebook,
        &args.sep,
        args.rc,
        args.fuzzy_sep != 0,
        args.anchor,
        args.strict_anchor != 0,
        args.strategy,
    );
    let final_stats = processor.process_paired_fastq_streaming(
        input_r1,
        input_r2,
        &output_table,
        &output_r2,
        args.batch_size,
        args.threads,
        stats,
        args.head,
        args.shape_qc != 0,
        // New Args
        match args.output_mode {
            OutputMode::All => 0,
            OutputMode::Good => 1,
            OutputMode::Warning => 2,
            OutputMode::Merge => 3,
        },
        args.readable_debug_example != 0,
        args.dump_fail,
    )?;

    let elapsed = start_time.elapsed();

    // 输出统计信息
    println!("\n=== 处理完成 ===");
    println!("总处理时间: {:.2}秒", elapsed.as_secs_f64());
    println!("总序列数: {}", final_stats.total_sequences);
    println!("正则匹配成功: {}", final_stats.regex_matched);
    println!("正则匹配失败: {}", final_stats.regex_failed);
    println!("组1匹配成功: {}", final_stats.group1_matched);
    println!("组2匹配成功: {}", final_stats.group2_matched);
    println!("模棱两可匹配: {}", final_stats.ambiguous_matched);
    println!("最终输出序列: {}", final_stats.final_output);
    println!(
        "处理速度: {:.0} 序列/秒",
        final_stats.total_sequences as f64 / elapsed.as_secs_f64()
    );

    // 显示简单的距离分布
    println!("\n--- 编辑距离分布 (Code1) ---");
    let mut dists1: Vec<_> = final_stats.group1_dist_counts.iter().collect();
    dists1.sort_by_key(|k| k.0);
    for (d, c) in dists1 {
        println!("Distance {}: {}", d, c);
    }

    println!("\n--- 编辑距离分布 (Code2) ---");
    let mut dists2: Vec<_> = final_stats.group2_dist_counts.iter().collect();
    dists2.sort_by_key(|k| k.0);
    for (d, c) in dists2 {
        println!("Distance {}: {}", d, c);
    }

    // 保存详细统计信息
    let summary_json = serde_json::to_string_pretty(&final_stats)?;
    std::fs::write(&json_summary, summary_json)?;
    println!("详细统计信息已保存到: {}", json_summary);

    Ok(())
}

fn load_codebook(path: Option<&str>) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let mut codebook = Vec::new();

    if let Some(p) = path {
        let file = File::open(p)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() == 2 {
                let code = parts[0].trim().to_string();
                let sequence = parts[1].trim().to_string();
                codebook.push((code, sequence));
            }
        }
    } else {
        // Use embedded codebook
        for line in DEFAULT_CODEBOOK_CSV.lines() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }

            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() == 2 {
                let code = parts[0].trim().to_string();
                let sequence = parts[1].trim().to_string();
                codebook.push((code, sequence));
            }
        }
    }

    Ok(codebook)
}
