#![feature(asm)] // Use assembly implementations of rank and select
#![allow(dead_code)]

use murmur3;
use std::cmp;
mod util;
use util::{
    bitarr::{b128, b64},
    bitrank, bitselect, popcnt,
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
    p: usize, // Fingerprint size
    q: usize, // Quotient size
    r: usize, // Remainder size

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
    fn new(n: usize, r: usize) -> Filter {
        Self::new_seeded(n, r, 0)
    }
    fn new_seeded(n: usize, r: usize, seed: u32) -> Filter {
        let nslots = cmp::max(64, util::nearest_pow_of_2(n));
        let nblocks = nslots / 64;
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
        (hash & ((1_u128 << self.q) - 1)) as usize
    }
    /// Get the k-th remainder for the hash
    fn calc_rem(&self, hash: u128, k: usize) -> Rem {
        // If hash can only make C remainders and k > C, let h_k = h_{k mod C}
        let capacity = (128 - self.q) / self.r;
        let k = k % capacity;
        // Get k-th chunk of r bits: bits in [a, b) where a=q+kr, b=q+(k+1)r
        let a = self.q + k * self.r;
        let b = a + self.r;
        assert!(b <= 128, "Remainder chunk overflowed 128 bits (b={})", b);
        let rem: Rem = ((hash & b128::half_open(a, b)) >> a) as Rem;
        // Don' let hash be 0
        if rem == 0 {
            rem + 1
        } else {
            rem
        }
    }

    /// Finds the absolute index of the rank-th runend past
    /// the start of the block_i-th block.
    /// Note: rank indexes from 0.
    fn multiblock_select(&self, block_i: usize, rank: usize) -> Option<usize> {
        assert!(block_i < self.nblocks, "Block index out of bounds");

        let mut rank: u64 = rank as u64;
        let mut step: usize; // select result
        let mut loc: usize = block_i * 64; // absolute index of runend for input
        let mut b: &Block; // block pointer

        // Step through each block, decrementing rank and incrementing loc
        // to count seen runends and advance to the correct block
        loop {
            b = &self.blocks[loc / 64];
            step = bitselect(b.runends, rank) as usize;
            loc += step;
            if step != 64 || loc >= self.nslots {
                break;
            }
            rank -= popcnt(b.runends); // account for seen runends
        }

        if loc >= self.nslots {
            None
        } else {
            Some(loc)
        }
    }

    /// Performs the blocked equivalent of the unblocked operation
    ///   select(Q.runends, rank(Q.occupieds, x)).
    fn rank_select(&self, x: usize) -> Option<usize> {
        assert!(x < self.nslots, "Absolute index out of range");
        let mut block_i = x / 64;
        let slot_i = x % 64;
        let mut b = &self.blocks[block_i];
        let mut rank = bitrank(b.occupieds, slot_i);

        // Exit early when the result of rank_select would be in a prior block.
        // This happens when
        // (1) slot is unoccupied and
        // (2) block offset is 0
        // (3) there are no runs before the slot in the block
        if !b.is_occupied(slot_i) && b.offset == 0 && rank == 0 {
            None
        } else {
            // Skip ahead to the block that offset is pointing inside
            let offset = b.offset % 64;
            block_i += b.offset / 64;
            b = &self.blocks[block_i];

            // Account for the runends before the offset in the current block
            rank += if offset > 0 {
                bitrank(b.runends, offset - 1)
            } else {
                0
            };
            // If rank(Q.occupieds, x) == 0, then select gives an invalid result
            if rank == 0 {
                None
            } else {
                // (rank-1) accounts for select's indexing from 0
                self.multiblock_select(block_i, (rank-1) as usize)
            }
        }
    }
}

fn main() {
    let block = Block::new();
    println!("block: {:?}", block);

    let filter = Filter::new(64 * 2 + 1, 4);
    println!("filter: {:?}", filter);
}

#[cfg(test)]
mod tests {
    use super::*;
    
    fn if_0_inc(x: Rem) -> Rem {
        if x == 0 { 1 } else { x }
    }
    fn mask_rem(x: u128, start: usize, end: usize) -> Rem {
        ((x & b128::half_open(start, end)) >> start) as Rem
    }

    #[test]
    fn test_calc_quot_rem() {
        let filter = Filter::new(64, 4);
        println!("r={}, q={}", filter.r, filter.q); // r=4, q=6
        let hash = 0x1234_ABCD_0000_0000__0000_0000_0000_0000_u128;
        for i in 0..(128-6)/4 {
            let rem = if_0_inc(mask_rem(hash, 6 + 4*i, 10 + 4*i));
            assert_eq!(filter.calc_rem(hash, i), rem);
        }
        // First remainder is 0, bumped up to 1 to satisfy nonzero rem invariant
        assert_eq!(filter.calc_rem(0, 0), 1);
        assert_eq!(filter.calc_rem(1, 0), 1);
        assert_eq!(filter.calc_rem(0b11_1111, 0), 1);
        // 7 set bits -> rem = 0b1
        assert_eq!(filter.calc_rem(0b111_1111, 0), 1);
        assert_eq!(filter.calc_rem(0b1111_1111, 0), 0b11);
        assert_eq!(filter.calc_rem(0b1_1111_1111, 0), 0b111);
        assert_eq!(filter.calc_rem(0b11_1111_1111, 0), 0b1111);

        let hash = 0b11_1100_0011_1100_0011_1111;
        assert_eq!(filter.calc_quot(hash), 0b11_1111);
        assert_eq!(filter.calc_rem(hash, 0), 1); // 0x0
        assert_eq!(filter.calc_rem(hash, 1), 0xf);
        assert_eq!(filter.calc_rem(hash, 2), 1);
        assert_eq!(filter.calc_rem(hash, 3), 0xf);
    }
    #[test]
    fn test_multiblock_select_single_block() {
        let mut filter = Filter::new(64, 4);
        let mut b;

        // Empty filter
        {
            b = &mut filter.blocks[0];
            b.occupieds = 0;
            b.runends = 0;
            assert_eq!(filter.multiblock_select(0, 0), None);
        }
        // Filter with one run
        {
            b = &mut filter.blocks[0];
            b.occupieds = 1;
            b.runends = 1;
            assert_eq!(filter.multiblock_select(0, 0), Some(0));
        }
        // Filter with multiple runs
        {
            b = &mut filter.blocks[0];
            // First two runs each have two elts, third run has one elt
            b.occupieds = 0b1101;
            b.runends  = 0b11010;
            assert_eq!(filter.multiblock_select(0, 0), Some(1));
            assert_eq!(filter.multiblock_select(0, 1), Some(3));
            assert_eq!(filter.multiblock_select(0, 2), Some(4));
            assert_eq!(filter.multiblock_select(0, 3), None);
        }
    }
    #[test]
    fn test_multiblock_select_multiple_blocks() {
        let mut filter = Filter::new(64*3, 4);
        
        // Filter with run starting in block 0 and ending in block 1
        {
            let b0 = &mut filter.blocks[0];
            b0.occupieds = 1;
            b0.runends = 0;

            let b1 = &mut filter.blocks[1];
            b1.occupieds = 0;
            b1.runends = 1;
            
            assert_eq!(filter.multiblock_select(0, 0), Some(64));
            assert_eq!(filter.multiblock_select(1, 0), Some(64));
            assert_eq!(filter.multiblock_select(2, 0), None);
        }

        // Filter with two runs:
        // A run starts in b0 and ends in b1, making b1 have nonzero offset
        {
            let b0 = &mut filter.blocks[0];
            b0.occupieds = 0b11;
            b0.runends = 0b01;

            let b1 = &mut filter.blocks[1];
            b1.occupieds = 0b10;
            b1.runends = 0b11;
            b1.offset = 1;
            
            assert_eq!(filter.multiblock_select(0, 0), Some(0));
            assert_eq!(filter.multiblock_select(0, 1), Some(64));
            assert_eq!(filter.multiblock_select(1, 0), Some(64));
            assert_eq!(filter.multiblock_select(1, 1), Some(65));
        }
    }
    #[test]
    fn test_rank_select_single_block() {
        let mut filter = Filter::new(64, 4);

        // Empty filter
        {
            // No setup needed: filter blank by default
            for i in 0..64 {
                assert_eq!(filter.rank_select(i), None)
            }
        }
        // Filter with one run
        {
            let b = &mut filter.blocks[0];
            b.occupieds = 1;
            b.runends = 1;
            assert_eq!(bitrank(1, 63), 1);
            for i in 0..64 {
                assert_eq!(filter.rank_select(i), Some(0), "i={}", i);
            }
        }
        println!("filter:{:?}", filter);
        // Filter with multiple runs
        {
            let b = &mut filter.blocks[0];
            // First two runs each have two elts, third run has one elt
            b.occupieds = 0b101001;
            b.runends   = 0b110010;
            assert_eq!(filter.rank_select(0), Some(1));
            assert_eq!(filter.rank_select(1), Some(1));
            assert_eq!(filter.rank_select(2), Some(1));
            assert_eq!(filter.rank_select(3), Some(4));
            assert_eq!(filter.rank_select(4), Some(4));
            assert_eq!(filter.rank_select(5), Some(5));
            for i in 6..64 {
                assert_eq!(filter.rank_select(i), Some(5));
            }
        }
    }
    #[test]
    fn test_rank_select_multi_block() {
        // Filter with run starting in block 0 and ending in block 1
        {
            let mut filter = Filter::new(64*3, 4);
            let b0 = &mut filter.blocks[0];
            b0.occupieds = 1;
            b0.runends = 0;
            b0.offset = 64;     // position of b0[0]'s runend

            let b1 = &mut filter.blocks[1];
            b1.occupieds = 0;
            b1.runends = 1;
            b1.offset = 1;      // position where b1's first run's end should go,
                                // i.e., position after runend from b0
            
            for i in 0..128 {
                assert_eq!(filter.rank_select(i), Some(64), "i={}", i);
            }
            for i in 128..filter.nslots {
                assert_eq!(filter.rank_select(i), None, "i={}", i);
            }
        }
        // Filter with two runs:
        // A run starts in b0 and ends in b1, making b1 have nonzero offset
        {
            let mut filter = Filter::new(64*3, 4);
            let b0 = &mut filter.blocks[0];
            b0.occupieds = 0b11;
            b0.runends = 0b01;
            b0.offset = 0;

            let b1 = &mut filter.blocks[1];
            b1.occupieds = 0b10;
            b1.runends = 0b11;
            b1.offset = 1;
            
            assert_eq!(filter.rank_select(0), Some(0));
            assert_eq!(filter.rank_select(1), Some(64));
            for i in 2..65 {
                assert_eq!(filter.rank_select(i), Some(64), "i={}", i);
            }
            assert_eq!(filter.rank_select(65), Some(65));
            for i in 66..128 {
                assert_eq!(filter.rank_select(i), Some(65), "i={}", i);
            }
            for i in 128..filter.nslots {
                assert_eq!(filter.rank_select(i), None, "i={}", i);
            }
        }
    }
}
