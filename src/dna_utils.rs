/// DNA序列处理工具函数

/// 构建正则表达式模式，使用捕获组直接提取三个组
pub fn build_regex_pattern(sep_config: &str, fuzzy_sep: bool, umi_on_head: bool) -> String {
    // 解析分隔符配置
    // 支持两种格式:
    // 1. 旧格式: 组1末尾,分隔符,组2开头 (例如: DD,CCC,DD)
    // 2. 新格式: 仅分隔符 (例如: GGG 或 CCC)
    
    let parts: Vec<&str> = sep_config.split(',').collect();
    
    let (group1_end, separator_raw, group2_start) = if parts.len() == 3 {
        (parts[0], parts[1], parts[2])
    } else {
        // 推断模式
        let sep = sep_config.trim();
        // 简单的推断逻辑：如果是 GGG 则用 HH，如果是 CCC 则用 DD
        if sep.contains('G') && !sep.contains('C') {
            ("HH", sep, "HH")
        } else if sep.contains('C') && !sep.contains('G') {
            ("DD", sep, "DD")
        } else {
            // 其他情况不加约束
            ("", sep, "")
        }
    };

    // If Right Mode (!umi_on_head), we are matching against RC sequence.
    // So separator must be RC'd.
    let separator_string = if !umi_on_head {
        reverse_complement(separator_raw)
    } else {
        separator_raw.to_string()
    };
    let separator = separator_string.as_str();

    let mut pattern = String::new();
    // pattern.push_str("^");

    // 定义 UMI 模式
    let umi_pattern = format!("({}{{9}})", convert_to_regex("N"));

    // 定义 Separator 模式 (包含 fuzzy 逻辑)
    let mut sep_pattern = String::new();
    if fuzzy_sep {
        // 允许 1 个错配
        if separator.len() <= 5 {
            sep_pattern.push_str("(?:"); // 非捕获组开始

            // 1. 原始精确匹配
            sep_pattern.push_str(&convert_to_regex(separator));

            // 2. 枚举每一位出错的情况
            for i in 0..separator.len() {
                sep_pattern.push('|');
                let prefix = &separator[0..i];
                let suffix = &separator[i + 1..];
                sep_pattern.push_str(&convert_to_regex(prefix));
                sep_pattern.push_str("[ATCGN]"); // 任意碱基
                sep_pattern.push_str(&convert_to_regex(suffix));
            }

            sep_pattern.push_str(")"); // 非捕获组结束
        } else {
            sep_pattern.push_str(&convert_to_regex(separator));
        }
    } else {
        sep_pattern.push_str(&convert_to_regex(separator));
    }

    // 定义 Group 1 模式
    let mut g1_pattern = String::new();
    g1_pattern.push_str("(");
    g1_pattern.push_str(&format!("{}{{13,14}}", convert_to_regex("N")));
    g1_pattern.push_str(&convert_to_regex(group1_end));
    g1_pattern.push_str(")");

    // 定义 Group 2 模式
    let mut g2_pattern = String::new();
    g2_pattern.push_str("(");
    g2_pattern.push_str(&convert_to_regex(group2_start));
    g2_pattern.push_str(&format!("{}{{13,14}}", convert_to_regex("N")));
    g2_pattern.push_str(")");

    // 根据 UMI 位置构建最终模式
    if umi_on_head {
        // 左模式: UMI + G1 + SEP + G2
        // G1 使用 group1_end, G2 使用 group2_start
        pattern.push_str(&umi_pattern);
        pattern.push_str(&g1_pattern);
        pattern.push_str(&sep_pattern);
        pattern.push_str(&g2_pattern);
    } else {
        // 右模式: G2' + SEP + G1' + UMI'
        // 注意：这里需要构造匹配 G2' 和 G1' 的模式
        // G2' 是 G2 的反向互补。G2 开头是 group2_start，所以 G2' 结尾是 RC(group2_start)
        // G1' 是 G1 的反向互补。G1 结尾是 group1_end，所以 G1' 开头是 RC(group1_end)

        let rc_g2_start = reverse_complement(group2_start);
        let rc_g1_end = reverse_complement(group1_end);

        // 构造 Part 1 (对应 G2')
        let mut p1_pattern = String::new();
        p1_pattern.push_str("(");
        p1_pattern.push_str(&format!("{}{{13,14}}", convert_to_regex("N")));
        p1_pattern.push_str(&convert_to_regex(&rc_g2_start));
        p1_pattern.push_str(")");

        // 构造 Part 2 (对应 G1')
        let mut p2_pattern = String::new();
        p2_pattern.push_str("(");
        p2_pattern.push_str(&convert_to_regex(&rc_g1_end));
        p2_pattern.push_str(&format!("{}{{13,14}}", convert_to_regex("N")));
        p2_pattern.push_str(")");

        pattern.push_str(&p1_pattern);
        pattern.push_str(&sep_pattern);
        pattern.push_str(&p2_pattern);
        pattern.push_str(&umi_pattern);
    }

    println!("pattern: {}", pattern);
    pattern
}

/// 将字符串转换为正则表达式
fn convert_to_regex(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            'A' => "A".to_string(),
            'T' => "T".to_string(),
            'C' => "C".to_string(),
            'G' => "G".to_string(),
            'N' => "[ATCGN]".to_string(),
            'D' => "[AGT]".to_string(),
            'H' => "[ATC]".to_string(),
            'W' => "[AT]".to_string(),
            'S' => "[CG]".to_string(),
            'M' => "[AC]".to_string(),
            'K' => "[GT]".to_string(),
            'R' => "[AG]".to_string(),
            'Y' => "[CT]".to_string(),
            'B' => "[CGT]".to_string(),
            'V' => "[ACG]".to_string(),
            _ => "[ATCGN]".to_string(),
        })
        .collect()
}

/// 计算反向互补序列
pub fn reverse_complement(seq: &str) -> String {
    let mut buffer = String::with_capacity(seq.len());
    reverse_complement_in_place(seq, &mut buffer);
    buffer
}

/// 将 DNA 序列编码为 u64
pub fn encode_dna(seq: &str) -> Option<u64> {
    if seq.len() > 30 {
        return None;
    }
    let mut val: u64 = 1; // Sentinel
    for b in seq.bytes() {
        val <<= 2;
        match b {
            b'A' => val |= 0,
            b'C' => val |= 1,
            b'G' => val |= 2,
            b'T' => val |= 3,
            _ => return None,
        }
    }
    Some(val)
}

/// 在原有buffer上计算反向互补序列，避免分配
pub fn reverse_complement_in_place(seq: &str, buffer: &mut String) {
    buffer.clear();
    buffer.reserve(seq.len());

    for c in seq.chars().rev() {
        let rc = match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            'R' => 'Y',
            'Y' => 'R',
            'M' => 'K',
            'K' => 'M',
            'W' => 'W',
            'S' => 'S',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            _ => c,
        };
        buffer.push(rc);
    }
}

/// 根据分隔符自动决定是否启用反向互补
pub fn should_enable_reverse_complement(sep_config: &str, rc_setting: i32) -> bool {
    match rc_setting {
        0 => false, // 明确禁用
        1 => true,  // 明确启用
        -1 => {
            // 自动决定
            let parts: Vec<&str> = sep_config.split(',').collect();
            let separator = if parts.len() >= 2 {
                parts[1]
            } else {
                sep_config.trim()
            };

            match separator {
                "CCC" => true,  // CCC分隔符启用反向互补
                "GGG" => false, // GGG分隔符不启用反向互补
                _ => {
                    // Try to guess based on content if not exact match
                    if separator.contains("CCC") {
                        true
                    } else if separator.contains("GGG") {
                        false
                    } else {
                        panic!("自动模式下，分隔符必须包含CCC或GGG，或者明确指定--rc参数")
                    }
                },
            }
        }
        _ => panic!("--rc参数必须是-1, 0, 或1"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use regex::Regex;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("CCCC"), "GGGG");
    }

    #[test]
    fn test_should_enable_reverse_complement() {
        // 测试自动决定功能
        assert!(should_enable_reverse_complement("DD,CCC,DD", -1));
        assert!(!should_enable_reverse_complement("DD,GGG,DD", -1));

        // 测试明确设置
        assert!(should_enable_reverse_complement("DD,CCC,DD", 1));
        assert!(!should_enable_reverse_complement("DD,CCC,DD", 0));
    }

    #[test]
    fn test_convert_to_regex() {
        assert_eq!(convert_to_regex("CCC"), "CCC");
        assert_eq!(convert_to_regex("DDD"), "[AGT][AGT][AGT]");
        assert_eq!(convert_to_regex("NNN"), "[ATCGN][ATCGN][ATCGN]");
    }
}
