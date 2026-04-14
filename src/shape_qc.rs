use image::{DynamicImage, ImageBuffer, Rgb};
use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

// 将常量公开，以便其他模块引用
pub const MATRIX_SIZE: usize = 5000;
const MATRIX_AREA: usize = MATRIX_SIZE * MATRIX_SIZE;

pub struct ShapeQC {
    pub good: Vec<AtomicU32>,
    pub ambiguous: Vec<AtomicU32>,
    pub low_quality: Vec<AtomicU32>,
    pub error_sep: Vec<AtomicU32>,
}

impl ShapeQC {
    pub fn new() -> Self {
        let mut good = Vec::with_capacity(MATRIX_AREA);
        let mut ambig = Vec::with_capacity(MATRIX_AREA);
        let mut low_q = Vec::with_capacity(MATRIX_AREA);
        let mut err_sep = Vec::with_capacity(MATRIX_AREA);
        for _ in 0..MATRIX_AREA {
            good.push(AtomicU32::new(0));
            ambig.push(AtomicU32::new(0));
            low_q.push(AtomicU32::new(0));
            err_sep.push(AtomicU32::new(0));
        }
        Self {
            good,
            ambiguous: ambig,
            low_quality: low_q,
            error_sep: err_sep,
        }
    }

    pub fn increment_good(&self, x: usize, y: usize) {
        if x < MATRIX_SIZE && y < MATRIX_SIZE {
            let idx = y * MATRIX_SIZE + x;
            self.good[idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn increment_ambiguous(&self, x: usize, y: usize) {
        if x < MATRIX_SIZE && y < MATRIX_SIZE {
            let idx = y * MATRIX_SIZE + x;
            self.ambiguous[idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn increment_low_quality(&self, x: usize, y: usize) {
        if x < MATRIX_SIZE && y < MATRIX_SIZE {
            let idx = y * MATRIX_SIZE + x;
            self.low_quality[idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn increment_error_sep(&self, x: usize, y: usize) {
        if x < MATRIX_SIZE && y < MATRIX_SIZE {
            let idx = y * MATRIX_SIZE + x;
            self.error_sep[idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn generate_images(
        &self,
        output_path: &str,
        mode: u8,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let dir = Path::new(output_path).parent().unwrap_or(Path::new("."));
        let prefix = "shape_qc";

        let good_counts: Vec<u32> = self
            .good
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect();
        let ambig_counts: Vec<u32> = self
            .ambiguous
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect();
        let low_q_counts: Vec<u32> = self
            .low_quality
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect();
        let err_sep_counts: Vec<u32> = self
            .error_sep
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect();

        let merged_counts: Vec<u32> = good_counts
            .iter()
            .zip(ambig_counts.iter())
            .map(|(a, b)| a + b)
            .collect();

        // 辅助函数：将 u32 转换为 f32 并应用高斯模糊，然后计算 Scale
        let process_and_save =
            |counts: &[u32], suffix: &str| -> Result<(), Box<dyn std::error::Error>> {
                // 1. 转换为 f32 矩阵
                let float_data: Vec<f32> = counts.iter().map(|&x| x as f32).collect();

                // 2. 应用高斯模糊 (Sigma = 2.0)
                // 我们直接在数据层面上模糊，而不是图像层面
                let blurred_data = gaussian_blur_f32(&float_data, MATRIX_SIZE, MATRIX_SIZE, 2.0);

                // 3. 计算 Scale (Quantile 0.99 of Log1p)
                // 注意：模糊后的数据可能是小数，也可能很小。
                // 应该在 Log1p 之前模糊？还是之后？
                // 用户说："在 counts 的时候计算 gaussian blur ... 在 blur 完了之后，计算 log1p"

                let mut log_vals = Vec::with_capacity(blurred_data.len());
                for &val in &blurred_data {
                    if val > 0.0 {
                        log_vals.push(val.ln_1p());
                    }
                }

                let max_val = if log_vals.is_empty() {
                    0.0
                } else {
                    log_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    let idx = (log_vals.len() as f64 * 0.99) as usize;
                    if idx >= log_vals.len() {
                        log_vals[log_vals.len() - 1]
                    } else {
                        log_vals[idx]
                    }
                };

                let scale = if max_val > 0.0 { 1.0 / max_val } else { 0.0 };

                // 4. 生成图像
                let img = ImageBuffer::from_fn(MATRIX_SIZE as u32, MATRIX_SIZE as u32, |x, y| {
                    let idx = (y as usize) * MATRIX_SIZE + (x as usize);
                    let val = blurred_data[idx]; // 已经模糊过的原始值
                    let log_val = val.ln_1p();

                    // Scale and Clip (0.0 - 1.0)
                    let mut scaled_f = log_val * scale;
                    if scaled_f > 1.0 {
                        scaled_f = 1.0;
                    }

                    get_reds_color(scaled_f as f64)
                });

                // 5. Flip Y (垂直翻转)
                let dynamic_img = DynamicImage::ImageRgb8(img);
                let flipped = dynamic_img.flipv();

                let file_name = if suffix == "MAIN" {
                    format!("{}.png", prefix)
                } else {
                    format!("{}_{}.png", prefix, suffix)
                };
                let out_path = dir.join(file_name);
                flipped.save(&out_path)?;
                println!("Saved Shape QC Image: {:?}", out_path);
                Ok(())
            };

        println!("Generating Shape QC images (Rust implementation)...");
        if let Err(e) = process_and_save(&good_counts, "good") {
            eprintln!("Warning: Failed to save shape_qc_good image: {}", e);
        }
        if let Err(e) = process_and_save(&ambig_counts, "ambiguous") {
            eprintln!("Warning: Failed to save shape_qc_ambiguous image: {}", e);
        }
        if let Err(e) = process_and_save(&merged_counts, "merged") {
            eprintln!("Warning: Failed to save shape_qc_merged image: {}", e);
        }
        if let Err(e) = process_and_save(&low_q_counts, "low_quality") {
            eprintln!("Warning: Failed to save shape_qc_low_quality image: {}", e);
        }
        if let Err(e) = process_and_save(&err_sep_counts, "error_sep") {
            eprintln!("Warning: Failed to save shape_qc_error_sep image: {}", e);
        }

        // Generate shape_qc.png based on output mode
        // 1=Good, 3=Merge. For All(0), we don't strictly define a "main" view, or maybe merged?
        // User requested explicitly for Merge and Good.
        let main_counts = match mode {
            1 => Some(&good_counts),   // Good
            3 => Some(&merged_counts), // Merge
            _ => None,
        };

        if let Some(counts) = main_counts {
            if let Err(e) = process_and_save(counts, "MAIN") {
                eprintln!("Warning: Failed to save shape_qc.png image: {}", e);
            }
        }

        Ok(())
    }
}

// 辅助函数：颜色映射
fn get_reds_color(value: f64) -> Rgb<u8> {
    let stops = [
        (0.0, Rgb([255, 245, 240])),
        (0.125, Rgb([254, 224, 210])),
        (0.25, Rgb([252, 187, 161])),
        (0.375, Rgb([252, 146, 114])),
        (0.5, Rgb([251, 106, 74])),
        (0.625, Rgb([239, 59, 44])),
        (0.75, Rgb([203, 24, 29])),
        (0.875, Rgb([165, 15, 21])),
        (1.0, Rgb([103, 0, 13])),
    ];

    let v = value.clamp(0.0, 1.0);
    for i in 0..stops.len() - 1 {
        let (p1, c1) = stops[i];
        let (p2, c2) = stops[i + 1];
        if v >= p1 && v <= p2 {
            let t = (v - p1) / (p2 - p1);
            let r = (c1[0] as f64 * (1.0 - t) + c2[0] as f64 * t) as u8;
            let g = (c1[1] as f64 * (1.0 - t) + c2[1] as f64 * t) as u8;
            let b = (c1[2] as f64 * (1.0 - t) + c2[2] as f64 * t) as u8;
            return Rgb([r, g, b]);
        }
    }
    stops[stops.len() - 1].1
}

use rayon::prelude::*;

/// 对扁平化的 2D 数组 (Vec<f32>) 进行高斯模糊
pub fn gaussian_blur_f32(data: &[f32], width: usize, height: usize, sigma: f32) -> Vec<f32> {
    if sigma <= 0.0 {
        return data.to_vec();
    }

    let radius = (sigma * 3.0).ceil() as usize;
    let kernel_size = 2 * radius + 1;
    let mut kernel = Vec::with_capacity(kernel_size);
    let mut sum = 0.0;

    let two_sigma_sq = 2.0 * sigma * sigma;

    for i in 0..kernel_size {
        let x = (i as isize - radius as isize) as f32;
        let val = (-x * x / two_sigma_sq).exp();
        kernel.push(val);
        sum += val;
    }

    for x in &mut kernel {
        *x /= sum;
    }

    // Horizontal Pass
    let temp: Vec<f32> = (0..height)
        .into_par_iter()
        .flat_map(|y| {
            let row_offset = y * width;
            let mut row_data = Vec::with_capacity(width);
            for x in 0..width {
                let mut val = 0.0;
                for k in 0..kernel_size {
                    let k_offset = k as isize - radius as isize;
                    let sample_x = (x as isize + k_offset).clamp(0, width as isize - 1) as usize;
                    val += data[row_offset + sample_x] * kernel[k];
                }
                row_data.push(val);
            }
            row_data
        })
        .collect();

    // Vertical Pass
    // We can't easily iterate cols in parallel and write to contiguous output without transposing or stride writes.
    // Easier: Iterate x (cols) in parallel, produce col data, then assemble?
    // Or iterate y (rows) of the *output* and sample vertically from temp.

    let output: Vec<f32> = (0..height)
        .into_par_iter()
        .flat_map(|y| {
            let mut row_data = Vec::with_capacity(width);
            for x in 0..width {
                let mut val = 0.0;
                for k in 0..kernel_size {
                    let k_offset = k as isize - radius as isize;
                    let sample_y = (y as isize + k_offset).clamp(0, height as isize - 1) as usize;
                    val += temp[sample_y * width + x] * kernel[k];
                }
                row_data.push(val);
            }
            row_data
        })
        .collect();

    output
}
