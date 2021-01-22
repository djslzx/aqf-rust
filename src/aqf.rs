use super::*;

#[derive(Debug)]
struct Block {
    remainders: [Rem; 64],
    occupieds: u64,
    runends: u64,
    extensions: u64,
    offset: usize,
}

#[derive(Debug)]
struct AQF {
    blocks: Vec<Block>,
    nblocks: usize,
    nslots: usize,

    // Sizes
    p: usize, // Fingerprint size
    q: usize, // Quotient size
    r: usize, // Remainder size

    seed: u32,

    // Remote representation
    remote: HashMap<(usize, Rem), (String, u128)>,
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
            extensions: 0,
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
    // TODO: selector encoding/decoding
}

impl AQF {
    fn new(n: usize, r: usize) -> AQF {
        Self::new_seeded(n, r, 0)
    }
    fn new_seeded(n: usize, r: usize, seed: u32) -> AQF {
        // add extra blocks for overflow
        let nblocks = cmp::max(1, util::nearest_pow_of_2(n)/64);
        let nslots = nblocks * 64;
        let q = (nslots as f64).log2() as usize;
        let mut blocks = Vec::with_capacity(nblocks);
        for _ in 0..nblocks {
            blocks.push(Block::new());
        }

        AQF {
            blocks,
            nblocks,
            nslots,
            q,
            r,
            p: q + r,
            seed,
            remote: HashMap::new(),
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
    fn add_block(&mut self) {
        let b = Block::new();
        self.blocks.push(b);
        self.nslots += 64;
        self.nblocks += 1;
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
        // TODO: calc_rem by adding k extra bits instead of taking k-size chunks

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
            step = bitselect(b.runends, if rank >= 64 { 63 } else { rank }) as usize;
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
        assert!(a < self.nslots && b < self.nslots,
                "Parameters out of bounds: a={}, b={}, nslots={}",
                a, b, self.nslots);
        // Exit early if invalid range
        if a > b {
            return;
        }
        // Start block_i at the first block after b, clamping it so it doesn't go off the end,
        // and work backwards
        let mut block_i = cmp::min(b/64 + 1, self.nblocks - 1);
        loop {
            // Account for direct/indirect offsets:
            // If direct, b[0] is occupied and target = x + offset
            // If indirect, b[0] is unoccupied and target = x + offset - 1
            let block_start = block_i * 64;
            let block = &mut self.blocks[block_i];
            // If block_i == 0, then the offset must be direct, as there's no valid indirect target
            if block_i == 0 && !block.is_occupied(0) {
                break;
            }
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
            if block_i == 0 { break; }
            let block_start = block_i * 64;
            let block = &mut self.blocks[block_i];
            let target = block_start + block.offset - 1;
            // If we've stepped back far enough, exit
            if target < loc {
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
        assert!(rem > 0);       // FIXME: not necessary?

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
                    None => {
                        // Extend the filter by one block
                        // and return the first empty index
                        self.add_block();
                        self.nslots - 64
                    }
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
                    "AQF failed to find runend (nslots={}, quot=(block={}, slot={}))",
                    self.nslots, quot/64, quot%64,
                ),
        }
    }
    /// Computes filter load factor
    fn load(&self) -> f64 {
        // TODO: fix this so it doesn't rely on remainders being nonzero
        let mut count = 0;
        for b in self.blocks.iter() {
            for rem in b.remainders.iter() {
                count += (*rem != 0) as i32;
            }
        }
        (count as f64)/(self.nslots as f64)
    }
}

impl Filter<String> for AQF {
    fn query(&self, elt: String) -> bool {
        let hash = self.hash(&elt[..]);
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
    fn insert(&mut self, elt: String) {
        let hash = self.hash(&elt[..]);
        let quot = self.calc_quot(hash);
        let rem = self.calc_rem(hash, 0);
        self.raw_insert(quot, rem);
        self.remote.insert((quot,rem), (elt, hash));
    }
}

#[cfg(test)]
mod tests {

}
