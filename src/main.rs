use std::cmp;
mod util;

#[derive(Debug)]
struct Block {
    remainders: [u8; 64],
    occupieds: u64,
    runends: u64,
    offset: u8,
}

#[derive(Debug)]
struct Filter {
    // Data
    blocks: Vec<Block>,

    // Sizes
    p: usize,                   // Fingerprint size
    q: usize,                   // Quotient size
    r: usize,                   // Remainder size
    nslots: usize,
    
    // Misc
    seed: i32,
}

impl Block {
    fn new() -> Block {
        Block {
            remainders: [0; 64],
            occupieds: 0,
            runends: 0,
            offset: 0,
        }
    }
}

impl Filter {
    fn new(n: usize, r: usize, seed: i32) -> Filter {
        let nslots = cmp::max(64, util::nearest_pow_of_2(n));
        let nblocks = nslots/64;
        let q = (nslots as f64).log2() as usize;
        let mut blocks = Vec::new();
        for _ in 0..nblocks {
            blocks.push(Block::new());
        }

        Filter {
            blocks: blocks,
            q: q,
            r: r,
            p: q + r,
            nslots: nslots,
            seed,
        }
    }
}

fn main() {
    let block = Block::new();
    println!("block: {:?}", block);
    let filter = Filter::new(64*2+1, 4, 0);
    println!("filter: {:?}", filter);
}
