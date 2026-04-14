use crate::dna_utils::encode_dna;
use crate::optimized_processor::CodebookEntry;
use dashmap::DashMap;
use rayon::prelude::*;
use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Read, Write};
use std::path::Path;

#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct TableEntry {
    pub key: u64,      // Encoded DNA sequence
    pub code_idx: u16, // Index in codebook
    pub dist: u8,      // Edit distance
    pub _pad: u8,      // Alignment padding
}

pub struct FinalLookupTable {
    pub codebook_hash: u64,
    pub table: Vec<TableEntry>,
}

impl FinalLookupTable {
    #[inline(always)]
    pub fn get(&self, key: u64) -> Option<(u16, u8)> {
        if self.table.is_empty() {
            return None;
        }

        let mask = self.table.len() - 1;
        let mut pos = (self.hash_key(key) & mask as u64) as usize;

        loop {
            let entry = &self.table[pos];
            if entry.key == key {
                return Some((entry.code_idx, entry.dist));
            }
            if entry.key == 0 {
                return None;
            }
            pos = (pos + 1) & mask;
        }
    }

    #[inline(always)]
    fn hash_key(&self, mut x: u64) -> u64 {
        // Fast 64-bit hash (fxhash-like)
        x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
        x ^ (x >> 31)
    }
}

pub fn generate_neighbors(seq: &str, max_dist: usize) -> Vec<(Vec<u8>, u8)> {
    let mut neighbors = HashMap::new();
    let seq_bytes = seq.as_bytes().to_vec();
    neighbors.insert(seq_bytes.clone(), 0u8);

    let mut current_gen = Vec::new();
    current_gen.push(seq_bytes);

    let bases = [b'A', b'C', b'G', b'T'];

    for d in 1..=max_dist {
        let mut next_gen = HashSet::new();
        for s in current_gen {
            // Substitutions
            for i in 0..s.len() {
                let original_char = s[i];
                for &b in &bases {
                    if b != original_char {
                        let mut next_s = s.clone();
                        next_s[i] = b;
                        if !neighbors.contains_key(&next_s) {
                            neighbors.insert(next_s.clone(), d as u8);
                            if d < max_dist {
                                next_gen.insert(next_s);
                            }
                        }
                    }
                }
            }
            // Deletions
            if s.len() > 1 {
                for i in 0..s.len() {
                    let mut next_s = Vec::with_capacity(s.len() - 1);
                    next_s.extend_from_slice(&s[..i]);
                    next_s.extend_from_slice(&s[i + 1..]);
                    if !neighbors.contains_key(&next_s) {
                        neighbors.insert(next_s.clone(), d as u8);
                        if d < max_dist {
                            next_gen.insert(next_s);
                        }
                    }
                }
            }
            // Insertions
            if s.len() < 30 {
                for i in 0..=s.len() {
                    for &b in &bases {
                        let mut next_s = Vec::with_capacity(s.len() + 1);
                        next_s.extend_from_slice(&s[..i]);
                        next_s.push(b);
                        next_s.extend_from_slice(&s[i..]);
                        if !neighbors.contains_key(&next_s) {
                            neighbors.insert(next_s.clone(), d as u8);
                            if d < max_dist {
                                next_gen.insert(next_s);
                            }
                        }
                    }
                }
            }
        }
        current_gen = next_gen.into_iter().collect();
    }

    neighbors.into_iter().collect()
}

pub fn calculate_codebook_hash(entries: &[CodebookEntry]) -> u64 {
    let mut hasher = DefaultHasher::new();
    for entry in entries {
        entry.code.hash(&mut hasher);
        entry.seq.hash(&mut hasher);
    }
    hasher.finish()
}

pub fn build_final_table(entries: &[CodebookEntry], max_dist: usize) -> FinalLookupTable {
    println!("Building neighbor map with edit distance {}...", max_dist);
    // Use a more memory-efficient approach: avoid DashMap of DashMaps
    // Instead, use a single DashMap with a combined key or a more compact structure if possible.
    // For now, let's keep the logic but optimize the neighbor generation and insertion.

    let neighbor_map = DashMap::new(); // encoded_seq -> Vec<(code_idx, dist)>

    entries.par_iter().enumerate().for_each(|(idx, entry)| {
        let neighbors = generate_neighbors(&entry.seq, max_dist);
        for (n_seq_bytes, dist) in neighbors {
            // Safety: n_seq_bytes only contains ACGT
            let n_seq = unsafe { std::str::from_utf8_unchecked(&n_seq_bytes) };
            if let Some(encoded) = encode_dna(n_seq) {
                let mut code_vec = neighbor_map.entry(encoded).or_insert_with(|| Vec::new());
                code_vec.push((idx as u16, dist));
            }
        }
    });

    println!("Analyzing neighbor map and building final table...");
    let mut temp_results = Vec::new();

    for (encoded_seq, code_dists) in neighbor_map {
        // Find min distance
        let mut min_dist = u8::MAX;
        for &(_, dist) in &code_dists {
            if dist < min_dist {
                min_dist = dist;
            }
        }

        let min_code_indices: Vec<u16> = code_dists
            .iter()
            .filter(|&&(_, dist)| dist == min_dist)
            .map(|&(idx, _)| idx)
            .collect();

        if min_code_indices.len() == 1 {
            temp_results.push((encoded_seq, min_code_indices[0], min_dist));
        } else {
            let mut code_vals: Vec<(i64, u16)> = min_code_indices
                .iter()
                .filter_map(|&idx| {
                    entries[idx as usize]
                        .code
                        .parse::<i64>()
                        .ok()
                        .map(|v| (v, idx))
                })
                .collect();

            if code_vals.len() == min_code_indices.len() {
                let min_v = code_vals.iter().map(|(v, _)| *v).min().unwrap();
                let max_v = code_vals.iter().map(|(v, _)| *v).max().unwrap();
                let range = max_v - min_v;

                if range <= 4 {
                    code_vals.sort_by_key(|(v, _)| *v);
                    let median_idx = code_vals[code_vals.len() / 2].1;
                    temp_results.push((encoded_seq, median_idx, min_dist));
                }
            }
        }
    }

    // Build Flat Hash Table
    let num_entries = temp_results.len();
    let capacity = (num_entries as f64 / 0.7).round() as usize;
    let capacity = capacity.next_power_of_two();
    let mask = (capacity - 1) as u64;

    let mut table = vec![
        TableEntry {
            key: 0,
            code_idx: 0,
            dist: 0,
            _pad: 0
        };
        capacity
    ];
    let final_table_obj = FinalLookupTable {
        codebook_hash: calculate_codebook_hash(entries),
        table: Vec::new(), // Temporary
    };

    for (key, code_idx, dist) in temp_results {
        let mut pos = (final_table_obj.hash_key(key) & mask) as usize;
        while table[pos].key != 0 {
            pos = (pos + 1) & (mask as usize);
        }
        table[pos] = TableEntry {
            key,
            code_idx,
            dist,
            _pad: 0,
        };
    }

    FinalLookupTable {
        codebook_hash: final_table_obj.codebook_hash,
        table,
    }
}

pub fn save_cache(table: &FinalLookupTable, path: &Path) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, file); // 8MB buffer

    // Header: Magic Number (8 bytes) + Codebook Hash (8 bytes) + Table Capacity (8 bytes)
    writer.write_all(b"SEQDEC02")?; // Version 2
    writer.write_all(&table.codebook_hash.to_le_bytes())?;
    writer.write_all(&(table.table.len() as u64).to_le_bytes())?;

    // Body: raw bytes of Vec<TableEntry>
    let bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(
            table.table.as_ptr() as *const u8,
            table.table.len() * std::mem::size_of::<TableEntry>(),
        )
    };
    writer.write_all(bytes)?;

    writer.flush()?;
    Ok(())
}

pub fn load_cache(path: &Path, codebook_hash: u64) -> Option<FinalLookupTable> {
    if !path.exists() {
        return None;
    }
    let mut file = File::open(path).ok()?;

    let mut header = [0u8; 24];
    file.read_exact(&mut header).ok()?;

    if &header[0..8] != b"SEQDEC02" {
        println!("Invalid cache format magic or version mismatch");
        return None;
    }

    let cached_hash = u64::from_le_bytes(header[8..16].try_into().unwrap());
    if cached_hash != codebook_hash {
        println!("Cache hash mismatch, recomputing...");
        return None;
    }

    let capacity = u64::from_le_bytes(header[16..24].try_into().unwrap()) as usize;
    let mut table = Vec::with_capacity(capacity);
    unsafe {
        table.set_len(capacity);
    }

    let bytes: &mut [u8] = unsafe {
        std::slice::from_raw_parts_mut(
            table.as_mut_ptr() as *mut u8,
            capacity * std::mem::size_of::<TableEntry>(),
        )
    };

    file.read_exact(bytes).ok()?;

    Some(FinalLookupTable {
        codebook_hash,
        table,
    })
}
