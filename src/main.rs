#![feature(asm)]                // Use assembly implementations of rank and select

mod util;
use util::{ B64, B128, popcnt, bitrank, bitselect };
use std::cmp;
use murmur3;

type Rem = u8;

#[derive(Debug)]
struct Block {
    remainders: [Rem; 64],
    occupieds: u64,
    runends: u64,
    offset: u8,
}

#[derive(Debug)]
struct Filter {
    // Data
    blocks: Vec<Block>,
    nblocks: usize,
    nslots: usize,

    // Sizes
    p: usize,                   // Fingerprint size
    q: usize,                   // Quotient size
    r: usize,                   // Remainder size
    
    // Misc
    seed: u32,
}

/// Error type for when multiblock_select fails to find the rank-th bit
#[derive(Debug, Clone)]
struct SelectOverflow;

impl Block {
    fn new() -> Block {
        Block {
            remainders: [0; 64],
            occupieds: 0,
            runends: 0,
            offset: 0,
        }
    }
    fn is_occupied(&self, i: usize) -> bool {
        B64::get(self.occupieds, i)
    }
    fn set_occupied(&mut self, to: bool, i: usize) {
        self.occupieds = B64::set_to(self.occupieds, to, i);
    }
    fn is_runend(&self, i: usize) -> bool {
        B64::get(self.runends, i)
    }
    fn set_runend(&mut self, to: bool, i: usize) {
        self.runends = B64::set_to(self.runends, to, i);
    }
}

impl Filter {
    fn new(n: usize, r: usize, seed: u32) -> Filter {
        let nslots = cmp::max(64, util::nearest_pow_of_2(n));
        let nblocks = nslots/64;
        let q = (nslots as f64).log2() as usize;
        let mut blocks = Vec::new();
        for _ in 0..nblocks {
            blocks.push(Block::new());
        }

        Filter {
            blocks,
            nblocks, 
            nslots,
            q,
            r,
            p: q + r,
            seed,
        }
    }
    fn hash(&self, word: &str) -> u128 {
        let ref mut b = word.as_bytes();
        match murmur3::murmur3_x64_128(b, self.seed) {
            Ok(v) => v,
            Err(_) => panic!("Failed to hash word"),
        }
    }
    fn calc_quot(&self, hash: u128) -> usize {
        (hash & ((1_u128 << self.q)-1)) as usize
    }
    /// Get the k-th remainder for the hash
    fn calc_rem(&self, hash: u128, k: usize) -> Rem {
        // Clamp k to [0, r)
        let k = if k >= self.r {
            k % self.r
        } else {
            k
        };
        // Get k-th chunk of r bits:
        // bits in [a, b) where a=q+kr, b=q+(k+1)r
        let a = self.q + k * self.r;
        let b = a + self.r;
        assert!(b > 128, "Remainder chunk overflowed 128 bits");
        let rem: Rem = ((hash & B128::half_open(a,b)) >> a) as Rem;
        // Don' let hash be 0
        if rem == 0 {
            rem + 1
        } else {
            rem 
        }
    }

    /// Finds the absolute index of the rank-th 1 bit in runends past 
    /// the start of the block_i-th block.
    /// Note: rank indexes from 1.
    fn multiblock_select(&self, block_i: usize, rank: usize) -> Result<usize, SelectOverflow> {
        assert!(block_i < self.nblocks, "Block index out of bounds");
        assert!(rank > 0, "Rank must index from 1");

        let mut rank: u64 = rank as u64;
        let mut b: &Block;                 // block pointer
        let mut step: usize;               // select result
        let mut loc: usize = block_i * 64; // absolute index of runend for input

        // Step through each block, decrementing rank and incrementing loc
        // to count seen runends and advance to the correct block
        loop {
            b = &self.blocks[loc/64]; 
            rank -= popcnt(b.runends); // account for seen runends
            step = bitselect(b.runends, rank-1) as usize; // use rank-1 b/c bitselect indexes from 0
            loc += step;
            
            if step != 64 || loc >= self.nslots {
                break;
            }
        }

        if loc >= self.nslots {
            Err(SelectOverflow {})
        } else {
            Ok(loc)
        }
    }

    /// Performs the blocked equivalent of the unblocked operation
    ///   select(Q.runends, rank(Q.occupieds, x)).
    fn find_runend(&self, x: usize) -> usize {
        0
    }
}

fn main() {
    let block = Block::new();
    println!("block: {:?}", block);

    let filter = Filter::new(64*2+1, 4, 0);
    println!("filter: {:?}", filter);

    let b = 0_u64;
    B64::set_to(b, true, 0);
    println!("bitarray: {:x?}", b);
    B64::set_to(b, true, 63);
    println!("bitarray: {:x?}", b);
}
