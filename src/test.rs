use std::{
    collections::HashSet,
    fs,
    io,
};
use rand::{
    Rng, SeedableRng, rngs::SmallRng
};
#[allow(unused_imports)]
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
    // let mut filter = RSQF::new(s, rem_size);
    let mut set = HashSet::new();
    let fstr = fs::read_to_string(path).unwrap();
    let mut lines = fstr.lines();

    // Insert `s` items into the filter and set
    eprintln!("Inserting {} items into filter...", s);
    for i in 0..s {
        let elt = lines.next().unwrap().to_string();
        set.insert(elt.clone());
        filter.insert(elt.clone());

        // Check that set and filter don't have false negatives right after insertions
        assert!(set.contains(&elt.clone()),
                "Set doesn't contain inserted elt={}", elt.clone());
        if !filter.query(elt.clone()) {
            eprintln!("Filter missing elt #{}", i);
            let hash = filter.hash(&elt.clone());
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash);
            // let ext = filter.calc_ext(hash, 1);
            let block_neighborhood = [
                filter.block(quot/64 - 1),
                filter.block(quot/64),
                filter.block(quot/64 + 1),
            ];
            panic!(
                "Filter doesn't contain elt={}, hash=0x{:x}; quot={} (block_i={}, slot_i={}), rem=0x{:x}, blocks={:#?};\
                 filter: q={}, r={}",
                elt, hash, quot, quot/64, quot%64, rem, block_neighborhood,
                filter.q(), filter.r(),
            );
        }
    }

    // Query the filter with the remaining items in the file
    eprintln!("Querying the filter with the remaining lines...", );
    let mut n_fps = 0;     // number of false positives
    let mut n_rfps = 0;    // number of repeated false positives
    let mut n_fns = 0;     // number of false negatives
    let mut n_queries = 0; // number of total queries
    let mut prev_fps = HashSet::new();
    for line in lines {
        n_queries += 1;
        let elt = line.to_string();
        let in_filter = filter.query(elt.clone());
        let in_set = set.contains(&elt.clone());
        if in_filter && !in_set { // false positive
            n_fps += 1;
            if prev_fps.contains(&elt.clone()) {
                n_rfps += 1;
            } else {
                prev_fps.insert(elt.clone());
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
             n_fps, (n_fps as f64)/(n_queries as f64) * 100_f64,
             n_rfps, (n_rfps as f64)/(n_queries as f64) * 100_f64);
    println!("False negatives: {} ({}%)",
             n_fns, (n_fns as f64)/(n_queries as f64) * 100_f64);
}

/// Test filter with randomly-generated queries
pub fn synthetic_test(a_s: f64, load: f64, queries_per_elt: usize, rem_size: usize, n_trials: usize) {
    println!("Running synthetic test with a_s={}, load={}, queries_per_elt={}, rem_size={}, n_trials={}",
             a_s, load, queries_per_elt, rem_size, n_trials);

    let n_slots = 1 << 14; // 16384
    let s = (util::nearest_pow_of_2(n_slots) as f64 * load) as usize;
    let a = (s as f64 * a_s) as usize;
    let n_queries = a * queries_per_elt;

    // Initialize member set
    eprintln!("Initializing filter and generating member set...");
    let mut rng = SmallRng::seed_from_u64(0);
    let mut set = HashSet::with_capacity(s);

    // Perform query trials
    eprintln!("Starting query trials...");
    let mut n_fps = 0;  // false positives
    let mut n_rfps = 0; // repeated false positives
    let mut n_fns = 0;  // false negatives
    let tot_queries = n_trials * n_queries;
    for i in 0..n_trials {
        // Initialize filter
        let mut filter = AQF::new(s, rem_size);
        // let mut filter = RSQF::new(s, rem_size);
        for _ in 0..s {
            let elt = rng.gen::<usize>().to_string();
            set.insert(elt.clone());
            filter.insert(elt.clone());
        }

        eprintln!("Performing query trial #{}...", i);
        let mut prev_fps = HashSet::new();
        let mut query_set: Vec<String> = Vec::with_capacity(a);
        for _ in 0..a {
            query_set.push(rng.gen::<usize>().to_string());
        }

        // Perform queries
        for _ in 0..n_queries {
            let elt = query_set.get(rng.gen::<usize>() % a).unwrap();
            let in_filter = filter.query(elt.clone());
            let in_set = set.contains(&elt.clone());
            if in_filter && !in_set {
                n_fps += 1;
                if prev_fps.contains(&elt.clone()) {
                    n_rfps += 1;
                } else {
                    prev_fps.insert(elt.clone());
                }
            } else if !in_filter && in_set {
                n_fns += 1;
            }
        }
        let queries_done = n_queries * (i+1);
        // eprintln!("Trial results:");
        // eprintln!("False positives: {} ({}%), repeated: {} ({}%)",
        //           n_fps, (n_fps as f64)/(queries_done as f64) * 100_f64,
        //           n_rfps, (n_rfps as f64)/(queries_done as f64) * 100_f64);
        // eprintln!("False negatives: {} ({}%)",
        //           n_fns, (n_fns as f64)/(queries_done as f64) * 100_f64);
    }
    println!("Test results:");
    println!("False positives: {} ({}%), repeated: {} ({}%)",
             n_fps, (n_fps as f64)/(tot_queries as f64) * 100_f64,
             n_rfps, (n_rfps as f64)/(tot_queries as f64) * 100_f64);
    println!("False negatives: {} ({}%)",
             n_fns, (n_fns as f64)/(tot_queries as f64) * 100_f64);
}