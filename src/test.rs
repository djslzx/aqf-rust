use std::{
    collections::HashSet,
    fs,
    io,
};
use crate::{util,
            rsqf::{
                RankSelectQuotientFilter,
                rsqf::RSQF,
            },
            aqf::AQF,
            Filter};

const DEBUG: bool = true;

/// Counts the number of unique lines in the file at `path`.
fn num_uniq_lines(path: &str) -> Result<usize, io::Error> {
    eprint!("Counting the number of lines in {}...", path);
    let mut set = HashSet::new();
    for line in fs::read_to_string(path)?.lines() {
        set.insert(line.to_string());
    }
    eprintln!(" done!");
    Ok(set.len())
}

/// File-input false positive rate test.
///
/// Using each of the lines in `path` as an element, split the file into query set A
/// and member set S.
/// Query the AQF using A/S close to `a_s`, with slight adjustments
/// to make the AQF's load factor closer to `load`.
pub fn file_fpr_test(path: &str, n_uniq_lines: usize, a_s: f64, load: f64, rem_size: usize) {
    eprintln!("Running file_fpr_test for {} ({} unique lines) with A/S={}, load={}, rem_size={}...",
              path, n_uniq_lines, a_s, load, rem_size);
    // 1/(1 + A/S) = S/(A+S)
    // n_uniq_lines = A+S
    let n_uniq_lines = n_uniq_lines as f64; //num_uniq_lines(path).unwrap() as f64;
    let s = (1.0/(1.0 + a_s) * n_uniq_lines) as usize;
    let s = (util::nearest_pow_of_2(s) as f64 * load) as usize;

    // Initialize filter and pull lines from file
    eprintln!("Initializing filter and pulling lines from file...");
    let mut filter = AQF::new(s, rem_size);
    let mut set = HashSet::new();
    let fstr = fs::read_to_string(path).unwrap();
    let mut lines = fstr.lines();

    // Insert `s` items into the filter and set
    eprintln!("Inserting {} items into filter...", s);
    for _ in 0..s {
        let elt = lines.next().unwrap().to_string();
        set.insert(elt.clone());
        filter.insert(elt.clone());
        assert!(set.contains(&elt.clone()),
                "Set doesn't contain inserted elt={}", elt.clone());
        if !filter.query(elt.clone()) {
            let hash = filter.hash(&elt.clone());
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash);
            let block_neighborhood = [
                filter.block(quot/64 - 1),
                filter.block(quot/64),
                filter.block(quot/64 + 1),
            ];
            eprintln!(
                "Filter doesn't contain elt={}; quot={} (block_i={}, slot_i={}), rem=0x{:x}, blocks={:#?}",
                elt, quot, quot/64, quot%64, rem, block_neighborhood,
            );
        }
    }

    // Query the filter with the remaining items in the file
    eprintln!("Querying the filter with the remaining lines...", );
    let mut n_fps = 0;     // number of false positives
    let mut n_rfps = 0;    // number of repeated false positives
    let mut n_fns = 0;     // number of false negatives
    let mut n_queries = 0; // number of total queries
    let mut seen = HashSet::new();
    for line in lines {
        n_queries += 1;
        let elt = line.to_string();
        let in_filter = filter.query(elt.clone());
        let in_set = set.contains(&elt.clone());
        if in_filter && !in_set { // false positive
            n_fps += 1;
            if seen.contains(&elt.clone()) {
                n_rfps += 1;
            } else {
                seen.insert(elt.clone());
            }
        } else if !in_filter && in_set { // false negative
            n_fns += 1;
            if DEBUG {
                let hash = filter.hash(&elt);
                let quot = filter.calc_quot(hash);
                let rem = filter.calc_rem(hash);
                let block = filter.block(quot/64);
                panic!(
                    "False negative on {}; quot={} (block_i={}, slot_i={}), rem=0x{:x}, block={:#?}",
                    elt, quot, quot/64, quot%64, rem, block,
                );
            }
        }
    }
    println!("False positives: {} ({}%), repeated: {} ({}%)",
             n_fps, (n_fps as f64)/(n_queries as f64),
             n_rfps, (n_rfps as f64)/(n_queries as f64));
    println!("False negatives: {} ({}%)",
             n_fns, (n_fns as f64)/(n_queries as f64));
}