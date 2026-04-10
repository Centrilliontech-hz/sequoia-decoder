use std::io::Write;

/// Helper module to format ambiguous alignment visualization
pub struct DebugFormatter;

impl DebugFormatter {
    pub fn format<W: Write>(
        writer: &mut W,
        seq_id: &str,
        original_seq: &str,
        anchor_found: bool,
        g1_info: Option<(usize, usize)>, // (start, len)
        g2_info: Option<(usize, usize)>, // (start, len)
        g1_candidates: &[String],
        g1_dists: &[usize],
        g1_codes: &[String], // Add codes
        g2_candidates: &[String],
        g2_dists: &[usize],
        g2_codes: &[String], // Add codes
    ) -> std::io::Result<()> {
        writeln!(writer, "ID: {}", seq_id)?;
        writeln!(writer, "SEQ: {}", original_seq)?;

        if !anchor_found {
            writeln!(writer, "MAP: Separator not found")?;
            writeln!(writer, "---")?;
            return Ok(());
        }

        // Construct Marker Line
        let len = original_seq.len();
        // Determine max length needed (original seq or markers)
        // Markers might extend if G2 is at the very end.
        let mut marker_line = vec![' '; len + 5]; // +5 buffer

        // Draw G1 marker
        if let Some((start, len)) = g1_info {
            if start < marker_line.len() {
                let end = (start + len).min(marker_line.len());
                if end > start {
                    for i in start..end {
                        marker_line[i] = '-';
                    }
                    marker_line[start] = '<';
                    if end > start {
                        marker_line[end - 1] = '>';
                    }
                    // Label "G1" in middle
                    let mid = start + len / 2;
                    if mid + 1 < end {
                        marker_line[mid] = 'G';
                        marker_line[mid + 1] = '1';
                    }
                }
            }
        }

        // Draw G2 marker
        if let Some((start, len)) = g2_info {
            if start < marker_line.len() {
                let end = (start + len).min(marker_line.len());
                if end > start {
                    for i in start..end {
                        marker_line[i] = '-';
                    }
                    marker_line[start] = '<';
                    if end > start {
                        marker_line[end - 1] = '>';
                    }
                    // Label "G2" in middle
                    let mid = start + len / 2;
                    if mid + 1 < end {
                        marker_line[mid] = 'G';
                        marker_line[mid + 1] = '2';
                    }
                }
            }
        }

        writeln!(
            writer,
            "MAP: {}",
            marker_line.iter().collect::<String>().trim_end()
        )?;

        // Draw Candidates
        let max_cand = g1_candidates.len().max(g2_candidates.len());

        // If no candidates but anchor found, print empty line?
        if max_cand == 0 {
            writeln!(writer, "     (No matching candidates found)")?;
        }

        for i in 0..max_cand {
            let mut line = vec![' '; marker_line.len()];

            // G1 Candidate
            if i < g1_candidates.len() {
                if let Some((start, _)) = g1_info {
                    let seq = &g1_candidates[i];
                    for (j, char) in seq.chars().enumerate() {
                        if start + j < line.len() {
                            line[start + j] = char;
                        }
                    }
                }
            }

            // G2 Candidate
            if i < g2_candidates.len() {
                if let Some((start, _)) = g2_info {
                    let seq = &g2_candidates[i];
                    for (j, char) in seq.chars().enumerate() {
                        if start + j < line.len() {
                            line[start + j] = char;
                        }
                    }
                }
            }

            write!(
                writer,
                "     {}",
                line.iter().collect::<String>().trim_end()
            )?;

            // Add metadata
            let c1_info = if i < g1_candidates.len() {
                format!("(G1: #{} D{})", g1_codes[i], g1_dists[i])
            } else {
                String::new()
            };
            let c2_info = if i < g2_candidates.len() {
                format!("(G2: #{} D{})", g2_codes[i], g2_dists[i])
            } else {
                String::new()
            };

            if !c1_info.is_empty() && !c2_info.is_empty() {
                writeln!(writer, "  {} {}", c1_info, c2_info)?;
            } else {
                writeln!(writer, "  {}{}", c1_info, c2_info)?;
            }
        }

        writeln!(writer, "---")?;
        Ok(())
    }
}
