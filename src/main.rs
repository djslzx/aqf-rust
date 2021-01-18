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

/// The possible results of the rank_select function
#[derive(Debug, PartialEq, Eq)]
enum RankSelectResult {
    Empty,                      // home slot unoccupied
    Full(usize),                // home slot occupied, last filled loc
    Overflow,                   // search went off the edge
}
use RankSelectResult::*;

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
    fn set_occupied(&mut self, i: usize, to: bool) {
        self.occupieds = b64::set_to(self.occupieds, to, i);
    }
    fn is_runend(&self, i: usize) -> bool {
        b64::get(self.runends, i)
    }
    fn set_runend(&mut self, i: usize, to: bool) {
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
    fn remainder(&self, x: usize) -> Rem {
        assert!(x < self.nslots, "x={}", x);
        self.blocks[x/64].remainders[x%64]
    }
    fn is_runend(&self, x: usize) -> bool {
        assert!(x < self.nslots, "x={}", x);
        self.blocks[x/64].is_runend(x%64)
    }
    fn is_occupied(&self, x: usize) -> bool {
        assert!(x < self.nslots, "x={}", x);
        self.blocks[x/64].is_occupied(x%64)
    }
    fn set_runend(&mut self, x: usize, to: bool) {
        assert!(x < self.nslots, "x={}", x);
        self.blocks[x/64].set_runend(x%64, to);
    }
    fn set_occupied(&mut self, x: usize, to: bool) {
        assert!(x < self.nslots, "x={}", x);
        self.blocks[x/64].set_occupied(x%64, to);
    }
    fn set_remainder(&mut self, x: usize, rem: Rem) {
        assert!(x < self.nslots, "x={}", x);
        self.blocks[x/64].remainders[x%64] = rem;
    }
    fn hash(&self, word: &str) -> u128 {
        let ref mut b = word.as_bytes();
        match murmur3::murmur3_x64_128(b, self.seed) {
            Ok(v) => v,
            Err(_) => panic!("Failed to hash word"),
        }
    }
    // Take first q bits of hash
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
    /// Returns None if no appropriate runend exists.
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
    ///   y = select(Q.runends, rank(Q.occupieds, x)).
    /// Note: x indexes from 0.
    ///
    /// Return behavior:
    /// - If y <= x, returns Empty
    /// - If y > x, returns Full(y)
    /// - If y runs off the edge, returns Overflow 
    fn rank_select(&self, x: usize) -> RankSelectResult {
        // Exit early if x is obviously out of range
        if x >= self.nslots {
            return Overflow
        }

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
            Empty
        } else {
            // Skip ahead to the block that offset is pointing inside
            let offset = b.offset % 64;
            block_i += b.offset / 64;
            // Handle overflowing offset (offset runs off the edge)
            if block_i >= self.nblocks {
                return Overflow;
            }
            b = &self.blocks[block_i];
            // Account for the runends before the offset in the current block
            rank += if offset > 0 {
                bitrank(b.runends, offset - 1)
            } else {
                0
            };
            // If rank(Q.occupieds, x) == 0, then there's nothing to see here
            if rank == 0 {
                Empty
            } else {
                // (rank-1) accounts for select's indexing from 0
                match self.multiblock_select(block_i, (rank-1) as usize) {
                    Some(loc) =>
                        if loc < x { 
                            Empty 
                        } else { 
                            Full(loc) 
                        },
                    None => Overflow,
                }
            }
        }
    }
    /// Finds the first unused slot at or after absloc x.
    /// Returns None if no slot found; 
    /// otherwise, returns Some(y) where y is the first open slot
    fn first_unused_slot(&self, x: usize) -> Option<usize> {
        let mut x = x;
        loop {
            match self.rank_select(x) {
                Empty => break Some(x),
                Full(loc) =>
                    if x <= loc {
                        x = loc + 1;
                    } else {
                        break Some(x);
                    }
                Overflow => break None,
            }
        }
    }
    /// Returns 1 if word is in the filter, 0 otherwise.
    fn contains(&self, word: &str) -> bool {
        let hash = self.hash(word);
        let quot = self.calc_quot(hash);
        let rem = self.calc_rem(hash, 0); // TODO: get 0-th rem for now
        
        if !self.is_occupied(quot) {
            false
        } else {
            if let Full(mut loc) = self.rank_select(quot) {
                loop {
                    // If matching remainder found, return true
                    if self.remainder(loc) == rem {
                        break true; 
                    }
                    // Stop when l < 0, l < quot, or l is a runend
                    if loc == 0 {
                        break false;
                    } else {
                        loc -= 1;
                        if loc < quot || self.is_runend(loc) {
                            break false;
                        }
                    }
                }
            } else {
                false
            }
        }
    }
    /// Shift the remainders and runends in [a,b] forward by 1 into [a+1, b+1]
    fn shift_remainders_and_runends(&mut self, a: usize, b: usize) {
        // TODO: use bit shifts instead of repeatedly masking
        for i in (a..=b).rev() {
            self.set_remainder(i+1, self.remainder(i));
            self.set_runend(i+1, self.is_runend(i));
            self.set_runend(i, false);
        }
    }
    /// Increment direct/indirect offsets with targets in [a,b]
    /// to reflect shifting remainders/runends in [a,b]
    fn inc_offsets(&mut self, a: usize, b: usize) {
        assert!(a <= b && b < self.nslots,
                "Parameters out of bounds: a={}, b={}, nslots={}",
                a, b, self.nslots);
        // Start block_i at the first block after b, clamping it so it doesn't go off the end,
        // and work backwards
        let mut block_i = cmp::min(b/64 + 1, self.nblocks - 1);
        loop {
            // Account for direct/indirect offsets:
            // If direct, b[0] is occupied and target = x + offset
            // If indirect, b[0] is unoccupied and target = x + offset - 1
            let block_start = block_i * 64;
            let block = &mut self.blocks[block_i];
            let target = block_start + block.offset - if block.is_occupied(0) { 0 } else { 1 };
            if target < a {
                // If we've stepped back far enough, exit
                break;
            } else if target <= b {
                // If a <= target <= b, increment offset
                block.offset += 1;
            }
            // If we're on the first block, we can't step back further: exit 
            if block_i == 0 {
                break;
            }
            // Step back by one block
            block_i -= 1;
        }
    }
    /// Increment indirect offsets targeting loc for quot's run
    fn inc_indirect_offsets(&mut self, quot: usize, loc: usize) {
        assert!(loc < self.nslots && quot < self.nslots,
                "Parameters out of bounds: quot={}, x={}",
                quot, loc);
        // Start block_i at the first block after b, clamping it so it doesn't go off the end
        let mut block_i = cmp::min(loc/64 + 1, self.nblocks - 1);
        // Walk through backwards from the first block after b
        loop {
            let block_start = block_i * 64;
            let block = &mut self.blocks[block_i];
            let target = block_start + block.offset - 1;
            // If we've stepped back far enough, exit
            if target < block_start {
                break;
            }
            // If target == loc, b[0] isn't occupied (offset is indirect), 
            // and quot < block_start (indirect offsets target runends 
            // that start in earlier blocks)
            if target == loc && !block.is_occupied(0) && quot < block_start {
                block.offset += 1;
            }
            // If we're on the first block, we can't step back further: exit 
            if block_i == 0 {
                break;
            }
            // Step back by one block
            block_i -= 1;
        }
    }
    /// Insert a (quot, rem) pair into filter
    fn raw_insert(&mut self, quot: usize, rem: Rem) {
        assert!(quot < self.nslots);
        assert!(rem > 0);

        // Find the appropriate runend
        match self.rank_select(quot) {
            Empty => {
                self.set_occupied(quot, true);
                self.set_runend(quot, true);
                self.set_remainder(quot, rem);
            }
            Full(r) => {
                // Find u, the first open slot after r, and
                // shift everything in [r+1, u-1] forward by 1 into [r+2, u], 
                // leaving r+1 writable
                let u = match self.first_unused_slot(r) {
                    Some(loc) => loc,
                    None => panic!("Couldn't find an unused slot"),
                };
                self.shift_remainders_and_runends(r+1, u-1);
                self.inc_offsets(r+1, u-1);
                // Start a new run or extend an existing one
                if !self.is_occupied(quot) {
                    // Set occupied, add runend, add rem, shift indirect offsets
                    self.set_occupied(quot, true);
                    self.set_runend(r+1, true);
                    self.set_remainder(r+1, rem);
                    self.inc_indirect_offsets(quot, r);
                } else {
                    // Don't need to set occupied
                    // Shift runend, add rem, shift offsets
                    self.set_runend(r, false);
                    self.set_runend(r+1, true);
                    self.set_remainder(r+1, rem);
                    self.inc_offsets(r, r);
                }
            }
            Overflow => 
                panic!(
                    "Filter ran out of space (nslots={}, quot=(block={}, slot={}))",
                    self.nslots, quot/64, quot%64,
                ),
        }
    }

    /// 
    fn insert(&mut self, word: &str) {}
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
                assert_eq!(filter.rank_select(i), Empty)
            }
        }
        // Filter with one run
        {
            let b = &mut filter.blocks[0];
            b.occupieds = 1;
            b.runends = 1;
            assert_eq!(filter.rank_select(0), Full(0));
            for i in 1..64 {
                assert_eq!(filter.rank_select(i), Empty, "i={}", i);
            }
        }
        // Filter with multiple runs
        {
            let b = &mut filter.blocks[0];
            // First two runs each have two elts, third run has one elt
            b.occupieds = 0b101001;
            b.runends   = 0b110010;
            assert_eq!(filter.rank_select(0), Full(1));
            assert_eq!(filter.rank_select(1), Full(1));
            assert_eq!(filter.rank_select(2), Empty);
            assert_eq!(filter.rank_select(3), Full(4));
            assert_eq!(filter.rank_select(4), Full(4));
            assert_eq!(filter.rank_select(5), Full(5));
            for i in 6..64 {
                assert_eq!(filter.rank_select(i), Empty);
            }
        }
    }
    #[test]
    fn test_rank_select_multi_block_1() {
        // Filter with run starting in block 0 and ending in block 1
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
        
        for i in 0..64 {
            assert_eq!(filter.rank_select(i), Full(64), "i={}", i);
        }
        for i in 65..filter.nslots {
            assert_eq!(filter.rank_select(i), Empty, "i={}", i);
        }
    }
    #[test]
    fn test_rank_select_multi_block_2() {
        // Filter with two runs:
        // A run starts in b0 and ends in b1, making b1 have nonzero offset
        let mut filter = Filter::new(64*3, 4);
        let b0 = &mut filter.blocks[0];
        b0.occupieds = 0b11;
        b0.runends = 0b01;
        b0.offset = 0;

        let b1 = &mut filter.blocks[1];
        b1.occupieds = 0b10;
        b1.runends = 0b11;
        b1.offset = 1;
        
        assert_eq!(filter.rank_select(0), Full(0));
        assert_eq!(filter.rank_select(1), Full(64));
        for i in 2..65 {
            assert_eq!(filter.rank_select(i), Full(64), "i={}", i);
        }
        assert_eq!(filter.rank_select(65), Full(65));
        for i in 66..filter.nslots {
            assert_eq!(filter.rank_select(i), Empty, "i={}", i);
        }
    }
    #[test]
    fn test_first_unused_empty() {
        // Empty filter
        // The first unused slot for a slot x should be x itself
        let filter = Filter::new(128, 4);
        for i in 0..filter.nslots {
            assert_eq!(filter.first_unused_slot(i), Some(i), "i={}", i);
        }
    }
    #[test]
    fn test_first_unused_single() {
        // Filter with one element at k
        // The first unused slot for a slot x should be 
        // k+1 if x == k, x otherwise (except when x=k=nslots-1)
        let nslots = 128;
        for k in 0..nslots {
            let mut filter = Filter::new(nslots, 4);
            filter.set_occupied(k, true);
            filter.set_runend(k, true);
            for i in 0..nslots {
                assert_eq!(
                    filter.first_unused_slot(i), 
                    if i == k {
                        if i == nslots-1 && k == nslots-1 { 
                            None 
                        } else { 
                            Some(k+1) 
                        }
                    } else {
                        Some(i)
                    },
                    "k={}, i={}", k, i,
                );
            }
        }
    }
    /// Insert a single run [a,b] into a filter of 128 slots 
    /// (doesn't handle overlaps with existing state)
    fn insert_run(filter: &mut Filter, a: usize, b: usize) {
        // Setup filter
        filter.set_occupied(a, true);
        filter.set_runend(b, true);
        // Set offset
        if a == 0 {
            // Set first block's offset if a = 0
            filter.blocks[0].offset = b;
        }
        if a < 64 && b >= 64 {
            // Set second block's offset if a in B[0] and b in B[1]
            filter.blocks[1].offset = b-63;
        }
    }
    #[test]
    fn test_first_unused_one_run() {
        // Filter with one multi-elt run from a to b
        // First unused for run [a,b] at slot x should be
        // b+1 if x in [a,b], otherwise x
        // Except when b=nslots-1, in which case result should be None
        let nslots = 128;
        for a in 0..nslots {
            for b in a..nslots {
                // Setup filter w/ one run in [a,b]
                let mut filter = Filter::new(nslots, 4);
                insert_run(&mut filter, a, b);
                // Test
                for i in 0..nslots {
                    assert_eq!(
                        filter.first_unused_slot(i),
                        if i < a || b < i {
                            Some(i)
                        } else {
                            if b == nslots-1 {
                                None
                            } else {
                                Some(b+1)
                            }
                        },
                        "run=[{}, {}], i={}", a, b, i,
                    );
                }
            }
        }
    }
    #[test]
    fn test_first_unused_two_runs() {
        // Two runs [a,b] and [c,d]
        let nslots = 128;

        // Test first_unused_slot for runs [a,b] and [c,d] nonoverlapping
        let test_two_runs = |a: usize, b: usize, c: usize, d: usize| {
            let mut filter = Filter::new(nslots, 4);
            insert_run(&mut filter, a, b);
            insert_run(&mut filter, c, d);

            for i in 0..nslots {
                assert_eq!(
                    filter.first_unused_slot(i),
                    if i < a || (i > b && i < c) || i > d {
                        // If i isn't in either interval
                        Some(i)
                    } else {
                        if i >= a && i <= b {
                            // If i is in the first interval
                            if c > b+1 {
                                // If there's a gap between the two intervals
                                Some(b+1)
                            } else if d < nslots - 1 {
                                // If the two intervals are connected
                                // and there's a gap after the second
                                Some(d+1)
                            } else {
                                None
                            }
                        } else if i >= c && i <= d && d < nslots-1 {
                            // If i is in the second interval and
                            // there's a gap after the second interval 
                            Some(d+1)
                        } else {
                            None
                        }
                    },
                    "r1=[{}, {}], r2=[{}, {}], i={}, filter={:?}",
                    a, b, c, d, i, filter,
                );
            }
        };
        // Run tests
        test_two_runs(0, 63, 65, 127); // center gap
        test_two_runs(0, 64, 66, 127); // center gap
        test_two_runs(1, 63, 64, 127); // left gap
        test_two_runs(1, 64, 65, 127); // left gap
        test_two_runs(0, 63, 64, 126); // right gap
        test_two_runs(0, 64, 65, 126); // right gap
        test_two_runs(1, 63, 65, 127); // left, center gap
        test_two_runs(1, 64, 66, 127); // left, center gap
        test_two_runs(1, 63, 64, 126); // left, right gap
        test_two_runs(1, 64, 65, 126); // left, right gap
        test_two_runs(0, 63, 65, 126); // center, right gap
        test_two_runs(0, 64, 66, 126); // center, right gap
        test_two_runs(1, 63, 65, 126); // left, center, right gap
        test_two_runs(1, 64, 66, 126); // left, center, right gap
    }
    #[test]
    /// Empty single-block filter 
    fn test_contains_empty_single() {
        let filter = Filter::new(64, 4);
        for i in 0..64 {
            let word = i.to_string();
            assert!(!filter.contains(&word), "word='{}'", word);
        }
    }
    #[test]
    /// 1-elt single-block filter
    fn test_contains_one_single() {
        let mut filter = Filter::new(64, 4);
        let word = "apples";
        
        // Manually insert 'apples'
        let hash = filter.hash(word);
        let quot = filter.calc_quot(hash);
        let rem = filter.calc_rem(hash, 0);
        let b = &mut filter.blocks[0];
        b.remainders[quot] = rem;
        b.set_occupied(quot, true);
        b.set_runend(quot, true);
        
        assert!(filter.contains("apples"));
    }
    #[test]
    /// Multi-elt single-block filter
    fn test_contains_multi_single() {
        let mut filter = Filter::new(64, 4);
        // I checked that these words' quot/rems don't conflict
        let words = ["apples", "bananas", "oranges"];

        // Manually insert words
        for word in words.iter() {
            let hash = filter.hash(word);
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash, 0);
            // println!("word={}, quot={:x}, rem={:x}", word, quot, rem);

            let b = &mut filter.blocks[0];
            b.remainders[quot] = rem;
            b.set_occupied(quot, true);
            b.set_runend(quot, true);
        }
        // Check that words are contained
        for word in words.iter() {
            assert!(filter.contains(word));
        }
    }
    #[test]
    // Empty multi-block filter
    fn test_contains_empty_multi() {
        let filter = Filter::new(64*3, 4);
        for i in 0..filter.nslots {
            let word = i.to_string();
            assert!(!filter.contains(&word), "word={}", word);
        }
    }
    #[test]
    // 1-elt multi-block filter
    fn test_contains_one_multi() {
        let mut filter = Filter::new(64*3, 4);
        let word = "apples";
        
        // Manually insert 'apples'
        let hash = filter.hash(word);
        let quot = filter.calc_quot(hash);
        let rem = filter.calc_rem(hash, 0);
        let b = &mut filter.blocks[quot/64];
        b.remainders[quot%64] = rem;
        b.set_occupied(quot%64, true);
        b.set_runend(quot%64, true);
        
        assert!(filter.contains("apples"));
    }
    #[test]
    // Multi-elt multi-block filter
    fn test_contains_multi_multi() {
        let mut filter = Filter::new(64*3, 4);
        let words = ["apples", "bananas", "oranges"];

        // Manually insert words
        for word in words.iter() {
            let hash = filter.hash(word);
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash, 0);
            // println!("word={}, quot={:x}, rem={:x}", word, quot, rem);

            let b = &mut filter.blocks[quot/64];
            b.remainders[quot%64] = rem;
            b.set_occupied(quot%64, true);
            b.set_runend(quot%64, true);
        }
        // Check that words are contained
        for word in words.iter() {
            assert!(filter.contains(word));
        }
    }
    #[test]
    fn test_shift_rem_and_runends() {
        let mut filter = Filter::new(128, 4);
        for i in 0..filter.nslots {
            filter.set_remainder(i, i as Rem);
            filter.set_runend(i, i % 3 == 0);
        }
        // Shift [0, nslots-2] to [1, nslots-1] and clear [0]
        filter.shift_remainders_and_runends(0, filter.nslots-2);
        assert_eq!(filter.remainder(0), 0);
        assert_eq!(filter.is_runend(0), false);
        for i in 1..filter.nslots {
            let j = (i as Rem) - 1;
            assert_eq!(filter.remainder(i), j);
            assert_eq!(filter.is_runend(i), j % 3 == 0);
        }
    }
    fn offset_state_init() -> Filter {
        let mut filter = Filter::new(64*4, 4);
        let b = &mut filter.blocks;
        // Run in b0: [0:(0,1)]
        b[0].set_occupied(0, true);
        b[0].set_runend(1, true);
        b[0].offset = 1;          // direct offset
        // Run in b0,b1: [63:(63,65)]
        b[0].set_occupied(63, true);
        b[1].set_runend(1, true);
        b[1].set_occupied(0, false);
        b[1].offset = 2;          // indirect offset
        // Run in b1: [67: (67,72)]
        b[1].set_occupied(3, true);
        b[1].set_runend(8, true);
        // Run in b1: [68: (73,73)]
        b[1].set_occupied(4, true);
        b[1].set_runend(9, true);
        // Run from b1 to b3: [94: (192)]
        b[1].set_occupied(30, true);
        b[2].offset = 65;     // dist from 128 to 192+1 (indirect)
        b[3].set_runend(0, true);
        b[3].offset = 1;      // dist from 192 to 192+1 (indirect)

        filter
    }
    #[test]
    fn test_inc_offsets() {
        // Inc all offsets [0, n-2] -> [1, n-1]
        {
            let mut filter = offset_state_init();
            filter.inc_offsets(0, filter.nslots-2);
            let b = filter.blocks;
            assert_eq!(b[0].offset, 2);
            assert_eq!(b[1].offset, 3);
            assert_eq!(b[2].offset, 66);
            assert_eq!(b[3].offset, 2);
        }
        // Inc ranges that nothing is pointing to
        {
            let mut filter = offset_state_init();
            let ranges = [(2,64), (66, 191), (193, 255)];
            for (start,end) in ranges.iter() {
                filter.inc_offsets(*start, *end);
                let b = &filter.blocks;
                // Check that for each inc_offsets call,
                // none of the offsets are changed
                assert_eq!(b[0].offset, 1);
                assert_eq!(b[1].offset, 2);
                assert_eq!(b[2].offset, 65);
                assert_eq!(b[3].offset, 1);
            }
        }
        // Inc the elts that things are pointing to
        {
            let mut filter = offset_state_init();
            filter.inc_offsets(1,1);
            assert_eq!(filter.blocks[0].offset, 2);
            assert_eq!(filter.blocks[1].offset, 2);
            assert_eq!(filter.blocks[2].offset, 65);
            assert_eq!(filter.blocks[3].offset, 1);
            filter.inc_offsets(65,65);
            assert_eq!(filter.blocks[0].offset, 2);
            assert_eq!(filter.blocks[1].offset, 3);
            assert_eq!(filter.blocks[2].offset, 65);
            assert_eq!(filter.blocks[3].offset, 1);
            filter.inc_offsets(192,192);
            assert_eq!(filter.blocks[0].offset, 2);
            assert_eq!(filter.blocks[1].offset, 3);
            assert_eq!(filter.blocks[2].offset, 66);
            assert_eq!(filter.blocks[3].offset, 2);
        }
    }
    #[test]
    fn test_inc_indirect_offsets() {
        // Inc direct targets does nothing
        {
            let mut filter = offset_state_init();
            filter.inc_indirect_offsets(0, 1);
            assert_eq!(filter.blocks[0].offset, 1);
            assert_eq!(filter.blocks[1].offset, 2);
            assert_eq!(filter.blocks[2].offset, 65);
            assert_eq!(filter.blocks[3].offset, 1);
        }
        // Inc indirect offsets but large quot does nothing
        {
            let mut filter = offset_state_init();
            filter.inc_indirect_offsets(66, 65);
            filter.inc_indirect_offsets(193, 192);
            assert_eq!(filter.blocks[0].offset, 1);
            assert_eq!(filter.blocks[1].offset, 2);
            assert_eq!(filter.blocks[2].offset, 65);
            assert_eq!(filter.blocks[3].offset, 1);
        }
        // Range with only indirect offsets does inc
        {
            let mut filter = offset_state_init();
            filter.inc_indirect_offsets(0, 65);
            filter.inc_indirect_offsets(0, 192);
            assert_eq!(filter.blocks[0].offset, 1);
            assert_eq!(filter.blocks[1].offset, 3);
            assert_eq!(filter.blocks[2].offset, 66);
            assert_eq!(filter.blocks[3].offset, 2);
        }
        // Range with both indirect/direct offsets
        // only increments indirect offsets
        {
            let mut filter = offset_state_init();
            for i in (0..filter.nslots).rev() {
                filter.inc_indirect_offsets(0, i);
            }
            assert_eq!(filter.blocks[0].offset, 1);
            assert_eq!(filter.blocks[1].offset, 3);
            assert_eq!(filter.blocks[2].offset, 66);
            assert_eq!(filter.blocks[3].offset, 2);
        }
    }
}
