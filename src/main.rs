#![feature(asm)]                // Use assembly implementations of rank and select
#![allow(dead_code)]

use std::cmp;
use murmur3;
mod util;
use util::{
    bitarr::{ b64, b128 },
    popcnt, bitrank, bitselect,
};

type Rem = u8;

#[derive(Debug)]
struct Block {
    remainders: [Rem; 64],
    occupieds: u64,
    runends: u64,
    offset: usize,
}

#[derive(Debug)]
struct Filter {
    blocks: Vec<Block>,
    nblocks: usize,
    nslots: usize,

    // Sizes
    p: usize,                   // Fingerprint size
    q: usize,                   // Quotient size
    r: usize,                   // Remainder size
    
    seed: u32,
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
    fn is_occupied(&self, i: usize) -> bool {
        b64::get(self.occupieds, i)
    }
    fn set_occupied(&mut self, to: bool, i: usize) {
        self.occupieds = b64::set_to(self.occupieds, to, i);
    }
    fn is_runend(&self, i: usize) -> bool {
        b64::get(self.runends, i)
    }
    fn set_runend(&mut self, to: bool, i: usize) {
        self.runends = b64::set_to(self.runends, to, i);
    }
}

impl Filter {
    fn new(n: usize, r: usize, seed: u32) -> Filter {
        let nslots = cmp::max(64, util::nearest_pow_of_2(n));
        let nblocks = nslots/64;
        let q = (nslots as f64).log2() as usize;
        let mut blocks = Vec::with_capacity(nblocks);
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
        // If hash can only make C remainders and k > C, let h_k = h_{k mod C}
        let capacity = (128 - self.q)/self.r;
        let k = k % capacity;
        // Get k-th chunk of r bits: bits in [a, b) where a=q+kr, b=q+(k+1)r
        let a = self.q + k * self.r;
        let b = a + self.r;
        assert!(b <= 128, "Remainder chunk overflowed 128 bits (b={})", b);
        let rem: Rem = ((hash & b128::half_open(a,b)) >> a) as Rem;
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
    fn multiblock_select(&self, block_i: usize, rank: usize) -> Option<usize> {
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
            None
        } else {
            Some(loc)
        }
    }
    // TODO: write multiblock_select using absolute index instead of block index

    /// Performs the blocked equivalent of the unblocked operation
    ///   select(Q.runends, rank(Q.occupieds, x)).
    fn find_runend(&self, x: usize) -> Option<usize> {
        assert!(x < self.nslots, "Absolute index out of range");
        let mut block_i = x/64;
        let slot_i = x%64;
        let mut b = &self.blocks[block_i];
        let mut rank = bitrank(b.occupieds, slot_i);

        // Exit early when 
        // (1) slot is unoccupied and
        // (2) block offset is before the slot (or offset is 0) and
        // (3) there are no runs before the slot in the block
        if !b.is_occupied(slot_i) &&
            (b.offset == 0 || b.offset < slot_i) &&
            rank == 0 
        {
            None
        } else {
            // Skip ahead to the block that offset is pointing inside
            let offset = b.offset%64;
            block_i += b.offset/64;
            b = &self.blocks[block_i];

            // Account for the runends before the offset in the current block
            rank += if offset > 0 { 
                bitrank(b.runends, offset-1)
            } else {
                0
            };
            assert!(rank > 0, "Rank should be positive for multiblock_select to work");
            match self.multiblock_select(block_i, rank as usize) {
                None => panic!("find_runend went off the edge"),
                Some(loc) => Some(loc)
            }
        }
    }

    
}

fn main() {
    let block = Block::new();
    println!("block: {:?}", block);

    let filter = Filter::new(64*2+1, 4, 0);
    println!("filter: {:?}", filter);

    println!("x={:x}", (1_i128 << 127) >> 10);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_rem() {
        let filter = Filter::new(64, 4, 0);
        println!("r={}, q={}", filter.r, filter.q); // r=4, q=6
        let hash = 0x1234_ABCD_0000_0000__0000_0000_0000_0000_u128;
        let first_rem = ((hash & b128::half_open(6,10)) >> 6) as Rem;
        let first_rem = if first_rem == 0 { 1 } else { first_rem };
        assert_eq!(filter.calc_rem(hash,0), first_rem, "first_rem=0x{:x}", first_rem);
        
    }
}
