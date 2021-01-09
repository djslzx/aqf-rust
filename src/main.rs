mod util;
use std::cmp;
use util::B64;

#[derive(Debug)]
struct Block {
    remainders: [u8; 64],
    occupieds: B64,
    runends: B64,
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
            occupieds: B64::zero(),
            runends: B64::zero(),
            offset: 0,
        }
    }
    fn is_occupied(&self, i: usize) -> bool {
        self.occupieds.get(i)
    }
    fn set_occupied(&mut self, to: bool, i: usize) {
        self.occupieds.set_to(to, i);
    }
    fn is_runend(&self, i: usize) -> bool {
        self.runends.get(i)
    }
    fn set_runend(&mut self, to: bool, i: usize) {
        self.runends.set_to(to, i);
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
            blocks,
            q,
            r,
            p: q + r,
            nslots,
            seed,
        }
    }
}

fn main() {
    let block = Block::new();
    println!("block: {:?}", block);

    let filter = Filter::new(64*2+1, 4, 0);
    println!("filter: {:?}", filter);

    let mut b = B64::zero();
    b.set_to(true, 0);
    println!("bitarray: {:x?}", b);
    b.set_to(true, 63);
    println!("bitarray: {:x?}", b);
}
