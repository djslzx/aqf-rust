use murmur3;
use crate::{
    Filter, Rem,
};
use std::cmp::{max, min};
use crate::util::{
    bitarr::{b128, b64},
    bitrank, bitselect, popcnt,
    nearest_pow_of_2
};
use std::fmt;

/// Abstraction for blocks used in RSQF
pub trait RankSelectBlock: fmt::Debug {
    // Required methods:
    fn new() -> Self;
    fn occupieds(&self) -> u64;
    fn is_occupied(&self, i: usize) -> bool;
    fn set_occupied(&mut self, i: usize, to: bool);
    fn runends(&self) -> u64;
    fn is_runend(&self, i: usize) -> bool;
    fn set_runend(&mut self, i: usize, to: bool);
    fn remainder(&self, i: usize) -> Rem;
    fn set_remainder(&mut self, i: usize, to: Rem);
    fn offset(&self) -> usize;
    fn inc_offset(&mut self);
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rems = (0..64)
            .map(|x| format!("0x{:x}", self.remainder(x)))
            .collect::<Vec<String>>()
            .chunks(8)
            .map(|c| c.join(" "))
            .collect::<Vec<String>>();
        f.debug_struct("Block")
            .field("remainders", &rems)
            .field("occupieds", &format_args!("[{:064b}]", self.occupieds().reverse_bits()))
            .field("runends  ", &format_args!("[{:064b}]", self.runends().reverse_bits()))
            .field("offset   ", &self.offset()).finish()
    }
}

/// The possible results of the rank_select function
#[derive(Debug, PartialEq, Eq)]
pub enum RankSelectResult {
    Empty,                      // home slot unoccupied
    Full(usize),                // home slot occupied, last filled loc
    Overflow,                   // search went off the edge
}

/// RSQF abstraction:
/// Any filter that implements the first set of methods
/// gets the second set for free.
pub trait RankSelectQuotientFilter {
    type Block: RankSelectBlock;

    // Required methods:
    fn new(n: usize, r: usize) -> Self;
    fn new_seeded(n: usize, r: usize, seed: u32) -> Self;
    fn q(&self) -> usize;
    fn r(&self) -> usize;
    fn nslots(&self) -> usize;
    fn nblocks(&self) -> usize;
    fn seed(&self) -> u32;
    fn block(&self, i: usize) -> &Self::Block;
    fn mut_block(&mut self, i: usize) -> &mut Self::Block;
    fn add_block(&mut self);

    // Received methods:
    // Metadata getters/setters
    fn is_occupied(&self, i: usize) -> bool {
        debug_assert!(i < self.nslots(), "Index out of range: {}/{}", i, self.nslots());
        self.block(i/64).is_occupied(i%64)
    }
    fn set_occupied(&mut self, i: usize, to: bool) {
        debug_assert!(i < self.nslots(), "Index out of range: {}/{}", i, self.nslots());
        self.mut_block(i/64).set_occupied(i%64, to);
    }
    fn is_runend(&self, i: usize) -> bool {
        debug_assert!(i < self.nslots(), "Index out of range: {}/{}", i, self.nslots());
        self.block(i/64).is_runend(i%64)
    }
    fn set_runend(&mut self, i: usize, to: bool) {
        debug_assert!(i < self.nslots(), "Index out of range: {}/{}", i, self.nslots());
        self.mut_block(i/64).set_runend(i%64, to);
    }
    fn remainder(&self, i: usize) -> Rem {
        debug_assert!(i < self.nslots(), "Index out of range: {}/{}", i, self.nslots());
        self.block(i/64).remainder(i%64)
    }
    fn set_remainder(&mut self, i: usize, to: Rem) {
        debug_assert!(i < self.nslots(), "Index out of range: {}/{}", i, self.nslots());
        self.mut_block(i/64).set_remainder(i%64, to);
    }
    /// Checks that representation invariants are met
    fn check_rep(&self) {
        // Check offsets
        for i in 0..self.nblocks() {
            let b = self.block(i);
            let b_start = i * 64;
            let runend = self.rank_select(b_start);
            if b.is_occupied(0) {
                assert_eq!(
                    runend, RankSelectResult::Full(b_start + b.offset()),
                    "B[0] occupied => offset points to B[0]'s runend; blocks[{},{}]={:#?}",
                    i, i + b.offset()/64,
                    (i..=(i+b.offset()/64))
                        .map(|j| self.block(j))
                        .collect::<Vec<&Self::Block>>(),
                );
            } else { // b[0] is unoccupied
                if b.offset() == 0 {
                    if b.is_runend(0) {
                        assert_eq!(
                            runend, RankSelectResult::Full(b_start),
                            "B[0] unoccupied, B.offset = 0, B.runend[0] = 1 => \
                             last prior run ends at B[0]; blocks[{}]={:#?}",
                            i, b,
                        );
                    } else {
                        assert_eq!(
                            runend, RankSelectResult::Empty,
                            "B[0] unoccupied, B.offset = 0, B.runend[0] = 0 => \
                             last prior run ends before B[0]; blocks[{}]={:#?}",
                            i, b,
                        );
                    }
                } else { // b.offset > 0
                    assert_eq!(
                        runend, RankSelectResult::Full(b_start + b.offset()),
                        "B[0] unoccupied, B.offset > 0 => \
                         prior runend ends at B.start + B.offset; blocks[{}]={:#?}",
                        i, b
                    );
                }
            }
        }
    }
    fn hash(&self, word: &str) -> u128 {
        let ref mut b = word.as_bytes();
        match murmur3::murmur3_x64_128(b, self.seed()) {
            Ok(v) => v,
            Err(_) => panic!("Failed to hash word"),
        }
    }
    /// Use first q bits of quotient
    fn calc_quot(&self, hash: u128) -> usize {
        (hash & b128::ones(self.q())) as usize
    }
    fn calc_rem(&self, hash: u128) -> Rem {
        ((hash & b128::half_open(self.q(), self.q() + self.r())) >> self.q()) as Rem
    }
    /// Finds the absolute index of the rank-th 1 bit past the start of B[`block_i`]
    /// in the metadata bits accessed from a block `b` using `f(b)`.
    /// Note: `rank` indexes from 0.
    fn select(&self, block_i: usize, rank: usize, f: fn(&Self::Block) -> u64) -> Option<usize> {
        assert!(block_i < self.nblocks(), "Block index out of bounds");

        let mut rank: u64 = rank as u64;
        let mut step: usize; // select result
        let mut loc: usize = block_i * 64; // absolute index of runend for input
        let nslots = self.nslots();

        // Step through each block, decrementing rank and incrementing loc
        // to count seen runends and advance to the correct block
        loop {
            let b = self.block(loc/64);
            let meta = f(b);
            step = bitselect(meta, if rank >= 64 { 63 } else { rank }) as usize;
            loc += step;
            if step != 64 || loc >= nslots {
                break;
            }
            rank -= popcnt(meta); // account for seen runends
        }

        if loc >= nslots {
            None
        } else {
            Some(loc)
        }
    }
    /// Finds the absolute index of the rank-th runend past
    /// the start of the block_i-th block.
    /// Note: rank indexes from 0.
    /// Returns None if no appropriate runend exists.
    fn select_runend(&self, block_i: usize, rank: usize) -> Option<usize> {
        self.select(block_i, rank, |b| b.runends())
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
        if x >= self.nslots() {
            return RankSelectResult::Overflow
        }

        let mut block_i = x / 64;
        let slot_i = x % 64;
        let mut b = self.block(block_i);

        // Compute i + O_i where i = x - (x mod 64)
        if !b.is_occupied(0) && b.offset() == 0 && !b.is_runend(0) {
            // b[0] unoccupied, b.offset = 0, b[0] not a runend =>
            // Negative offset
            if slot_i == 0 {
                return RankSelectResult::Empty;
            }
        } else {
            // Non-negative offset
            if slot_i == 0 {
                return RankSelectResult::Full(block_i * 64 + b.offset());
            }
            block_i += b.offset()/64;
        }

        // Handle overflowing offset (offset runs off the edge)
        if block_i >= self.nblocks() {
            return RankSelectResult::Overflow;
        }
        // Count the number of occupieds between i+1 (b.start+1) and j (x),
        // excluding the (potential) quot at b[0] by subtracting off 1 if b[0] is occupied
        let mut d = bitrank(b.occupieds(), slot_i) - (b.is_occupied(0) as u64);
        // Advance offset to relevant value for the block that b.offset points to
        let offset = b.offset() % 64;
        b = self.block(block_i);

        // Account for the runends in [0, offset]
        d += bitrank(b.runends(), offset);

        // If rank(Q.occupieds, x) == 0, then there's nothing to see here
        if d == 0 {
            RankSelectResult::Empty
        } else {
            // (rank-1) accounts for select's indexing from 0
            match self.select_runend(block_i, (d -1) as usize) {
                Some(loc) =>
                    if loc < x {
                        RankSelectResult::Empty
                    } else {
                        RankSelectResult::Full(loc)
                    },
                None => RankSelectResult::Overflow,
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
                RankSelectResult::Empty => break Some(x),
                RankSelectResult::Full(loc) =>
                    if x <= loc {
                        x = loc + 1;
                    } else {
                        break Some(x);
                    }
                RankSelectResult::Overflow => break None,
            }
        }
    }
    fn apply_to_run_with_start<F>(&mut self, quot: usize, loc: usize, mut f: F)
        where F: FnMut(&mut Self, usize) {
        let mut i = loc;
        loop {
            f(self, i);
            if i == 0 || i == quot || self.is_runend(i-1) {
                break;
            } else {
                i -= 1;
            }
        }
        self.check_rep();
    }
    /// Applies `f(quot, pos)` to all slots in the run associated with `quot`.
    /// TODO: turn into a macro?
    fn apply_to_run<F>(&mut self, quot: usize, f: F)
        where F: FnMut(&mut Self, usize) {
        if self.is_occupied(quot) {
            if let RankSelectResult::Full(loc) = self.rank_select(quot) {
                self.apply_to_run_with_start(quot, loc, f);
            }
        }
    }
    /// Applies `f(quot, pos)` to all slots in the `block_i`-th block.
    fn apply_to_block<F>(&mut self, block_i: usize, mut f: F)
        where F: FnMut(&mut Self, usize, usize) {
        if let Some((mut q, mut i)) = self.last_intersecting_run(block_i) {
            let block_start = block_i * 64;

            // Skip forward until we hit the block boundaries
            while i >= block_start + 64 {
                i -= 1;
                if self.is_runend(i) {
                    // If we encounter a runend at i < end,
                    // then there must be another quotient q' < q
                    q = self.prev_quot(q).unwrap();
                }
            }
            // Step backwards through each slot and apply f
            'block: loop {
                'run: loop {
                    f(self, q, i); // apply f
                    if i == block_start {
                        break 'block;
                    } else if i == q || self.is_runend(i-1) {
                        break 'run;
                    } else {
                        i -= 1;
                    }
                }
                // Jump backwards to the previous run,
                // exiting if the runend of the run is in an earlier block
                match self.prev_pair(q, i, block_i) {
                    Some((prev_q, prev_end)) => {
                        q = prev_q;
                        i = prev_end;
                    }
                    None => break 'block
                }
            }
        }
        self.check_rep();
    }
    /// Determines the quotient and runend of the last run to intersect
    /// the `block_i`-th block.
    fn last_intersecting_run(&self, block_i: usize) -> Option<(usize, usize)> {
        let b = self.block(block_i);
        let block_start = block_i * 64;
        let n = popcnt(b.occupieds());
        let mut q: usize;
        if n == 0 {
            if b.offset() == 0 && !b.is_runend(0) {
                // Negative offset
                return None
            } else {
                q = match self.prev_quot(block_start) {
                    Some(loc) => loc,
                    None => panic!("Couldn't find previous quotient for a block with nonzero offset"),
                };
            }
        } else {
            // n > 0 => get last occupied quotient
            q = block_start + b64::highest_set_bit(b.occupieds());
        }
        // Get q's corresponding runend
        let mut end = match self.rank_select(q) {
            RankSelectResult::Full(loc) => loc,
            _ => panic!("q={} should have a runend", q),
        };
        if end <= block_start + 63 {
            // Exit if the first run we find (the rightmost valid candidate)
            // is guaranteed to overlap with the block
            Some((q, end))
        } else {
            // Otherwise, work backwards to find the first run whose runend
            // doesn't exceed the end of the block, then choose the run
            // after that one
            let mut last_q;
            let mut last_end;
            loop {
                // Cache current (q, end)
                last_q = q;
                last_end = end;
                // Get previous occupied pair
                match self.prev_pair(q, end, block_i) {
                    Some((new_q, new_end)) => {
                        q = new_q;
                        end = new_end;
                    }
                    // If there is no valid previous pair, then the current one
                    // is the last intersecting run
                    None => return Some((last_q, last_end))
                };
                // Exit w/ the old run when the new run ends before the end of the block
                if end < block_start + 63 {
                    break Some((last_q, last_end))
                }
            }
        }
    }
    /// Find the previous occupied quotient and runend pair.
    /// Exits early and returns `None` if the search goes past `bound`,
    /// a lower bound on block index.
    fn prev_pair(&self, q: usize, end: usize, bound_i: usize) -> Option<(usize, usize)> {
        debug_assert!(end/64 >= bound_i);
        match self.prev_end(end, bound_i) {
            Some(prev_end) => {
                match self.prev_quot(q) {
                    Some(prev_quot) => Some((prev_quot, prev_end)),
                    None => panic!("Found prev_end but failed to find prev_quot; \
                                    q={}, end={}, prev_end={}, bound_i={}",
                                   q, end, prev_end, bound_i),
                }
            }
            None => None
        }
    }
    /// Find the position of the last runend before position `end`,
    /// returning `None` if no runend can be found in or after the `bound_i`-th block.
    fn prev_end(&self, end: usize, bound_i: usize) -> Option<usize> {
        let block_i = end/64;
        debug_assert!(block_i >= bound_i, "block_i={}, bound_i={}", block_i, bound_i);
        let mut b = self.block(block_i);
        // Count the number of runends in the block excluding `end`
        let n = if end%64 == 0 { 0 } else { bitrank(b.runends(), end%64 - 1) };
        if n == 0 {
            for i in (bound_i..block_i).rev() {
                b = self.block(i);
                // Exit loop if b has any runends
                if b.runends() != 0 {
                    // Absolute index of last runend in b
                    return Some(i*64 + b64::highest_set_bit(b.runends()))
                }
            }
            // Reached bound without finding runends
            None
        } else {
            // When n > 0, we pick the second-to-last runend
            Some(block_i*64 + bitselect(b.runends(), n-1) as usize)
        }
    }
    /// Find the position of the last occupied quotient before position `q`.
    /// Returns `None` if `q` is the first occupied quotient.
    fn prev_quot(&self, q: usize) -> Option<usize> {
        let block_i = q/64;
        let mut b = self.block(block_i);
        // Count the number of occupied quotients excluding q
        let n = if q%64 == 0 { 0 } else { bitrank(b.occupieds(), q%64 - 1) };
        if n == 0 {
            // When there are no other occupied quotients in the block, look through
            // prior blocks.
            for i in (0..block_i).rev() {
                b = self.block(i);
                // Exit if b has occupied quotients
                if b.occupieds() != 0 {
                    // Absolute index of last occupied quotient in b
                    return Some(i*64 + b64::highest_set_bit(b.occupieds()));
                }
            }
            // Reached filter start without finding quotients
            None
        } else {
            // When n > 0, we pick the second-to-last quotient
            Some(block_i*64 + bitselect(b.occupieds(), n-1) as usize)
        }
    }
    /// Shift the remainders and runends in [a,b] forward by 1 into [a+1, b+1]
    fn shift_remainders_and_runends(&mut self, a: usize, b: usize) {
        // TODO: use bit shifts instead of repeatedly masking
        for i in (a..=b).rev() {
            self.set_remainder(i+1, self.remainder(i));
            self.set_runend(i+1, self.is_runend(i));
        }
        self.set_runend(a, false);
        self.check_rep();
    }
    /// Increment direct/indirect offsets with targets in [a,b]
    /// to reflect shifting remainders/runends in [a,b]
    fn inc_offsets(&mut self, a: usize, b: usize) {
        assert!(a < self.nslots() && b < self.nslots(),
                "Parameters out of bounds: a={}, b={}, nslots={}",
                a, b, self.nslots());
        // Exit early if invalid range
        if a > b {
            return;
        }
        // Start block_i at the first block after b, clamping it so it doesn't go off the end,
        // and work backwards
        let mut block_i = min(b/64 + 1, self.nblocks() - 1);
        loop {
            // Account for direct/indirect offsets:
            // If direct, b[0] is occupied and target = x + offset
            // If indirect, b[0] is unoccupied and target = x + offset - 1
            let block_start = block_i * 64;
            let block = self.mut_block(block_i);
            // If block_i == 0, then the offset must be direct, as there's no valid indirect target
            if block_i == 0 && !block.is_occupied(0) {
                break;
            }
            let target = block_start + block.offset() - if block.is_occupied(0) { 0 } else { 1 };
            if target < a {
                // If we've stepped back far enough, exit
                break;
            } else if target <= b {
                // If a <= target <= b, increment offset
                block.inc_offset();
            }
            // If we're on the first block, we can't step back further: exit
            if block_i == 0 {
                break;
            }
            // Step back by one block
            block_i -= 1;
        }
        self.check_rep();
    }
    /// Increment indirect offsets targeting loc for quot's run
    fn inc_indirect_offsets(&mut self, quot: usize, loc: usize) {
        assert!(loc < self.nslots() && quot < self.nslots(),
                "Parameters out of bounds: quot={}, x={}",
                quot, loc);
        // Start block_i at the first block after b, clamping it so it doesn't go off the end
        let mut block_i = min(loc/64 + 1, self.nblocks() - 1);
        // Walk through backwards from the first block after b
        loop {
            if block_i == 0 { break; }
            let block_start = block_i * 64;
            let block = self.mut_block(block_i);
            let target = block_start + block.offset() - 1;
            // If we've stepped back far enough, exit
            if target < loc {
                break;
            }
            // If target == loc, b[0] isn't occupied (offset is indirect),
            // and quot < block_start (indirect offsets target runends
            // that start in earlier blocks)
            if target == loc && !block.is_occupied(0) && quot < block_start {
                block.inc_offset();
            }
            if block_i == 0 {
                // If we're on the first block, we can't step back further: exit
                break;
            } else {
                // Step back by one block
                block_i -= 1;
            }

        }
        self.check_rep();
    }
}

pub mod rsqf {
    use super::*;
    use std::fmt;

    pub struct Block {
        remainders: [Rem; 64],
        occupieds: u64,
        runends: u64,
        offset: usize,
    }

    #[derive(Debug)]
    pub struct RSQF {
        blocks: Vec<Block>,
        nblocks: usize,
        nslots: usize,   // number of slots (nblocks * 64)
        nelts: usize,    // number of inserted elements

        // Sizes
        p: usize, // Fingerprint size
        q: usize, // Quotient size
        r: usize, // Remainder size

        seed: u32,
    }

    impl RankSelectBlock for Block {
        fn new() -> Block {
            Block {
                remainders: [0; 64],
                occupieds: 0,
                runends: 0,
                offset: 0,
            }
        }
        fn occupieds(&self) -> u64 {
            self.occupieds
        }
        fn is_occupied(&self, i: usize) -> bool {
            b64::get(self.occupieds, i)
        }
        fn set_occupied(&mut self, i: usize, to: bool) {
            self.occupieds = b64::set_to(self.occupieds, to, i);
        }
        fn runends(&self) -> u64 {
            self.runends
        }
        fn is_runend(&self, i: usize) -> bool {
            b64::get(self.runends, i)
        }
        fn set_runend(&mut self, i: usize, to: bool) {
            self.runends = b64::set_to(self.runends, to, i);
        }
        fn remainder(&self, i: usize) -> Rem {
            self.remainders[i]
        }
        fn set_remainder(&mut self, i: usize, to: Rem) {
            self.remainders[i] = to;
        }
        fn offset(&self) -> usize {
            self.offset
        }
        fn inc_offset(&mut self) {
            self.offset += 1;
        }
    }
    impl fmt::Debug for Block {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            RankSelectBlock::fmt(self, f)
        }
    }

    impl RankSelectQuotientFilter for RSQF {
        type Block = Block;

        fn new(n: usize, r: usize) -> RSQF {
            Self::new_seeded(n, r, 0)
        }
        fn new_seeded(n: usize, r: usize, seed: u32) -> RSQF {
            // add extra blocks for overflow
            let nblocks = max(1, nearest_pow_of_2(n) / 64);
            let nslots = nblocks * 64;
            let q = (nslots as f64).log2() as usize;
            let mut blocks = Vec::with_capacity(nblocks);
            for _ in 0..nblocks {
                blocks.push(Block::new());
            }

            RSQF {
                blocks,
                nblocks,
                nslots,
                nelts: 0,
                q,
                r,
                p: q + r,
                seed,
            }
        }
        fn q(&self) -> usize {
            self.q
        }
        fn r(&self) -> usize {
            self.r
        }
        fn nslots(&self) -> usize {
            self.nslots
        }
        fn nblocks(&self) -> usize {
            self.nblocks
        }
        fn seed(&self) -> u32 {
            self.seed
        }
        fn block(&self, i: usize) -> &Block {
            &self.blocks[i]
        }
        fn mut_block(&mut self, i: usize) -> &mut Block {
            &mut self.blocks[i]
        }
        fn add_block(&mut self) {
            let b = Block::new();
            self.blocks.push(b);
            self.nslots += 64;
            self.nblocks += 1;
            self.check_rep();
        }
    }

    impl RSQF {
        /// Get the k-th remainder for the hash
        fn calc_kth_rem(&self, hash: u128, k: usize) -> Rem {
            let q = self.q();
            let r = self.r();

            // If hash can only make C remainders and k > C, let h_k = h_{k mod C}
            let capacity = (128 - q) / r;
            let k = k % capacity;
            // Get k-th chunk of r bits: bits in [a, b) where a=q+kr, b=q+(k+1)r
            let a = q + k * r;
            let b = a + r;
            assert!(b <= 128, "Remainder chunk overflowed 128 bits (b={})", b);
            let rem: Rem = ((hash & b128::half_open(a, b)) >> a) as Rem;
            // Don' let hash be 0
            if rem == 0 {
                1
            } else {
                rem
            }
        }
        /// Insert a (quot, rem) pair into filter
        fn raw_insert(&mut self, quot: usize, rem: Rem) {
            assert!(quot < self.nslots);
            self.nelts += 1;

            // Find the appropriate runend
            match self.rank_select(quot) {
                RankSelectResult::Empty => {
                    self.set_occupied(quot, true);
                    self.set_runend(quot, true);
                    self.set_remainder(quot, rem);
                }
                RankSelectResult::Full(r) => {
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
                RankSelectResult::Overflow =>
                    panic!(
                        "RSQF failed to find runend (nslots={}, quot=(block={}, slot={}))",
                        self.nslots, quot/64, quot%64,
                    ),
            }
            self.check_rep();
        }
        /// Computes filter load factor
        pub fn load(&self) -> f64 {
            (self.nelts as f64)/(self.nslots as f64)
        }
    }

    impl Filter<String> for RSQF {
        fn query(&mut self, elt: String) -> bool {
            let hash = self.hash(&elt[..]);
            let quot = self.calc_quot(hash);
            let rem = self.calc_rem(hash); // TODO: get 0-th rem for now

            let result = if !self.is_occupied(quot) {
                false
            } else {
                if let RankSelectResult::Full(mut loc) = self.rank_select(quot) {
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
            };
            self.check_rep();
            result
        }
        fn insert(&mut self, elt: String) {
            let hash = self.hash(&elt[..]);
            let quot = self.calc_quot(hash);
            let rem = self.calc_rem(hash);
            self.raw_insert(quot, rem);
            self.check_rep();
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use rand::{
            Rng, SeedableRng, rngs::SmallRng,
        };
        use std::collections::{HashSet, HashMap};

        fn if_0_inc(x: Rem) -> Rem {
            if x == 0 { 1 } else { x }
        }
        fn mask_rem(x: u128, start: usize, end: usize) -> Rem {
            ((x & b128::half_open(start, end)) >> start) as Rem
        }

        #[test]
        fn test_calc_quot_kth_rem() {
            let filter = RSQF::new(64, 4);
            let (r, q) = (4, 6);
            let hash = 0x1234_ABCD_0000_0000__0000_0000_0000_0000_u128;
            for i in 0..(128-6)/4 {
                let rem = if_0_inc(mask_rem(hash, q + r*i, q+r*(i+1)));
                assert_eq!(filter.calc_kth_rem(hash, i), rem);
            }
            // First remainder is 0, bumped up to 1 to satisfy nonzero rem invariant
            assert_eq!(filter.calc_kth_rem(0, 0), 1);
            assert_eq!(filter.calc_kth_rem(1, 0), 1);
            assert_eq!(filter.calc_kth_rem(0b11_1111, 0), 1);
            // 7 set bits -> rem = 0b1
            assert_eq!(filter.calc_kth_rem(0b111_1111, 0), 1);
            assert_eq!(filter.calc_kth_rem(0b1111_1111, 0), 0b11);
            assert_eq!(filter.calc_kth_rem(0b1_1111_1111, 0), 0b111);
            assert_eq!(filter.calc_kth_rem(0b11_1111_1111, 0), 0b1111);

            let hash = 0b11_1100_00_11_1100_0011_1111;
            assert_eq!(filter.calc_quot(hash), 0b11_1111);
            assert_eq!(filter.calc_kth_rem(hash, 0), 1); // 0x0
            assert_eq!(filter.calc_kth_rem(hash, 1), 0xf);
            assert_eq!(filter.calc_kth_rem(hash, 2), 1);
            assert_eq!(filter.calc_kth_rem(hash, 3), 0xf);
        }
        #[test]
        fn test_calc_rem() {
            let filter = RSQF::new(64, 4);
            let hash = 0b1001_1100_111111;
            assert_eq!(filter.calc_quot(hash), 0b111111);
            assert_eq!(filter.calc_rem(hash), 0b1100, "left={:b}, right={:b}", filter.calc_rem(hash), 0b1100);
            assert_eq!(filter.calc_kth_rem(hash,0), 0b1100);
            assert_eq!(filter.calc_kth_rem(hash,1), 0b1001);
        }
        #[test]
        fn test_prev_q_single() {
            // Empty case
            {
                let filter = RSQF::new(64 * 3, 4);
                for i in 0..filter.nslots {
                    assert_eq!(filter.prev_quot(i), None);
                }
            }
            // Single quotient at i, query with q=j
            {
                let nslots = 64 * 3;
                for i in 0..nslots {
                    let mut filter = RSQF::new(nslots, 4);
                    filter.set_occupied(i, true);

                    for j in 0..nslots {
                        assert_eq!(
                            filter.prev_quot(j),
                            if j > i { Some(i) } else { None },
                            "i={}, j={}",
                            i, j,
                        );
                    }
                }
            }
        }
        #[test]
        fn test_prev_q_double() {
            // Two quotients at i <= j
            let nslots = 64 * 3;
            for i in 0..nslots {
                for j in i..nslots {
                    let mut filter = RSQF::new(nslots, 4);
                    filter.set_occupied(i, true);
                    filter.set_occupied(j, true);

                    assert_eq!(filter.prev_quot(i), None,
                               "q=i, i={}, j={}", i, j);
                    assert_eq!(filter.prev_quot(j),
                               if j > i { Some(i) } else { None },
                               "q=j, i={}, j={}", i, j);
                }
            }
        }
        #[test] //expensive test (~1.5s)
        #[ignore]
        fn test_prev_q_triple() {
            // Three quotients at i <= j <= k
            let nslots = 64*3;
            for i in 0..nslots {
                for j in i..nslots {
                    for k in j..nslots {
                        let mut filter = RSQF::new(nslots, 4);
                        filter.set_occupied(i, true);
                        filter.set_occupied(j, true);
                        filter.set_occupied(k, true);

                        assert_eq!(filter.prev_quot(i), None,
                                   "q=i, i={}, j={}, k={}",
                                   i, j, k);
                        assert_eq!(filter.prev_quot(j),
                                   if j > i { Some(i) } else { None },
                                   "q=j, i={}, j={}, k={}",
                                   i, j, k);
                        assert_eq!(filter.prev_quot(k),
                                   if k > j { Some(j) }
                                   else if k > i { Some(i) }
                                   else { None },
                                   "q=k, i={}, j={}, k={}",
                                   i, j, k);
                    }
                }
            }
        }
        #[test]
        fn test_prev_end_single() {
            // Empty case
            {
                let filter = RSQF::new(64 * 3, 4);
                for i in 0..filter.nslots {
                    assert_eq!(filter.prev_end(i, 0), None);
                }
            }
            // Single runend at i, query with end=j
            {
                let nslots = 64 * 3;
                for i in 0..nslots {
                    let mut filter = RSQF::new(nslots, 4);
                    filter.set_runend(i, true);

                    for j in 0..nslots {
                        assert_eq!(
                            filter.prev_end(j, 0),
                            if j > i { Some(i) } else { None },
                            "i={}, j={}",
                            i, j,
                        );
                        assert_eq!(
                            filter.prev_end(j, j / 64),
                            if j / 64 > i / 64 { None } // j is in a later block, i is before bound
                            else if j > i { Some(i) } // j is after i but in the same block
                            else { None }, // j <= i
                            "i={}, j={}",
                            i, j,
                        );
                    }
                }
            }
        }
        #[test]
        fn test_prev_end_double() {
            // Two runends at i <= j
            let nslots = 64*3;
            for i in 0..nslots {
                for j in i..nslots {
                    let mut filter = RSQF::new(nslots, 4);
                    filter.set_runend(i, true);
                    filter.set_runend(j, true);

                    // i <= j
                    assert_eq!(filter.prev_end(i, 0),
                               None,
                               "end=i, i={}, j={}", i, j);
                    assert_eq!(filter.prev_end(j, 0),
                               if j > i { Some(i) } else { None },
                               "end=j, i={}, j={}", i, j);
                    assert_eq!(filter.prev_end(i, i/64),
                               None,
                               "end=i, i={}, j={}", i, j);
                    assert_eq!(filter.prev_end(j, j/64),
                               if j > i && j/64 == i/64 { Some(i) } else { None },
                               "end=i, i={}, j={}", i, j);
                }
            }
        }
        #[ignore] //expensive
        #[test]
        fn test_prev_end_triple() {
            // Three quotients at i <= j <= k
            {
                let nslots = 64*3;
                for i in 0..nslots {
                    for j in i..nslots {
                        for k in j..nslots {
                            let mut filter = RSQF::new(nslots, 4);
                            filter.set_runend(i, true);
                            filter.set_runend(j, true);
                            filter.set_runend(k, true);

                            //bound_i = 0
                            assert_eq!(filter.prev_end(i, 0),
                                       None,
                                       "end=i, i={}, j={}, k={}",
                                       i, j, k);
                            assert_eq!(filter.prev_end(j, 0),
                                       if j > i { Some(i) } else { None },
                                       "end=j, i={}, j={}, k={}",
                                       i, j, k);
                            assert_eq!(filter.prev_end(k, 0),
                                       if k > j { Some(j) }
                                       else if k > i { Some(i) }
                                       else { None },
                                       "end=k, i={}, j={}, k={}",
                                       i, j, k);

                            //bound_i = end/64
                            assert_eq!(filter.prev_end(i, i/64),
                                       None,
                                       "end=i, i={}, j={}, k={}",
                                       i, j, k);
                            assert_eq!(filter.prev_end(j, j/64),
                                       if j > i && j/64 == i/64 { Some(i) } else { None },
                                       "end=j, i={}, j={}, k={}",
                                       i, j, k);
                            assert_eq!(filter.prev_end(k, k/64),
                                       if k > j && k/64 == j/64 { Some(j) }
                                       else if k > i && k/64 == i/64 { Some(i) }
                                       else { None },
                                       "end=k, i={}, j={}, k={}",
                                       i, j, k);
                        }
                    }
                }
            }
        }
        #[test]
        fn test_last_intersecting_run_empty() {
            // Empty block -> no last intersecting run
            let filter = RSQF::new(64*3, 4);
            for i in 0..filter.nblocks {
                assert_eq!(filter.last_intersecting_run(i), None);
            }
        }
        #[test]
        fn test_last_intersecting_run_prev() {
            // Block with an intersecting run from an earlier block
            for i in 0..64 {
                let mut filter = RSQF::new(64*3, 4);
                // Run: [i, [i, i+64]]
                filter.set_occupied(i, true);
                filter.set_runend(i+64, true);
                filter.blocks[0].offset = if i == 0 { 64 } else { 0 };
                filter.blocks[1].offset = i;
                for block_i in 0..=1 {
                    assert_eq!(filter.last_intersecting_run(block_i),
                               Some((i, i+64)),
                               "i={}, block_i={}", i, block_i);
                }
            }
        }
        #[test]
        fn test_last_intersecting_run_prev_spanning() {
            // Block completely covered by a run from an earlier block
            for i in 0..64 {
                let mut filter = RSQF::new(64*3, 4);
                // Run spanning three blocks: [i: [i, i+128]]
                filter.set_occupied(i, true);
                filter.set_runend(i+128, true);
                filter.blocks[0].offset = if i == 0 { 128 } else { 0 };
                filter.blocks[1].offset = i+64;
                filter.blocks[2].offset = i;
                for block_i in 0..=2 {
                    assert_eq!(filter.last_intersecting_run(block_i),
                               Some((i, i+128)));
                }
            }
        }
        #[test]
        fn test_last_intersecting_run_mult_quots() {
            // 3-block filter, multiple quotients in the middle block.
            let mut filter = RSQF::new(64*3, 4);
            // (1) Run in only block 0 [0,32]
            filter.set_occupied(0, true);
            filter.set_runend(32, true);
            filter.blocks[0].offset = 32;
            // (2) Run between block 0 and block 1 [16,64]
            filter.set_occupied(16, true);
            filter.set_runend(64, true);
            filter.blocks[1].offset = 0;
            // (3) Run in only block 1 [65, 69]
            filter.set_occupied(65, true);
            filter.set_runend(69, true);
            // (4) Run going from block 1 to 2 [127, 132]
            filter.set_occupied(127, true);
            filter.set_runend(132, true);
            filter.blocks[2].offset = 4;
            // (5) Run in only block 2 [129, 135]
            filter.set_occupied(129, true);
            filter.set_runend(135, true);

            assert_eq!(
                filter.last_intersecting_run(0),
                Some((16, 64)),
                "block[0]={:#?}",
                filter.blocks[0],
            );
            assert_eq!(
                filter.last_intersecting_run(1),
                Some((127, 132)),
                "block[0]={:#?}\nblock[1]={:#?}",
                filter.blocks[0], filter.blocks[1],
            );
            assert_eq!(
                filter.last_intersecting_run(2),
                Some((129, 135)),
                "block[1]={:#?}\nblock[2]={:#?}",
                filter.blocks[1], filter.blocks[2],
            );
        }
        fn apply_collect_runs(filter: &mut RSQF, block_i: usize) -> HashMap<usize, usize> {
            // Store mappings of the form (slot -> quot)
            let mut runs: HashMap<usize, usize> = HashMap::new();
            let collect_runs = |_: &mut RSQF, quot: usize, i: usize| {
                runs.insert(i, quot);
            };
            filter.apply_to_block(block_i, collect_runs);
            runs
        }
        #[test]
        fn test_apply_to_block() {
            let mut filter = RSQF::new(64*3, 4);
            // Run from end of b0 to start of b1 [60: [60, 68]]
            filter.set_occupied(60, true);
            filter.set_runend(68, true);
            filter.blocks[1].offset = 4;
            // Run in b1 [64: [69,76]]
            filter.set_occupied(64, true);
            filter.set_runend(76, true);
            filter.blocks[1].offset = 12;
            // Run in b1 [70: [77,110]]
            filter.set_occupied(70, true);
            filter.set_runend(110, true);
            // Run in b1 [90: [111,126]]
            filter.set_occupied(90, true);
            filter.set_runend(126, true);
            // Run in b1 [120: [127,129]]
            filter.set_occupied(120, true);
            filter.set_runend(129, true);
            filter.blocks[2].offset = 1;

            let runs = apply_collect_runs(&mut filter, 1);
            //println!("runs={:#?}", runs);

            // Check that the expected slot->quot mappings are in place
            for i in 64..128 {
                assert_eq!(runs[&i],
                           if i <= 68 { 60 }
                           else if i <= 76 { 64 }
                           else if i <= 110 { 70 }
                           else if i <= 126 { 90 }
                           else { 120 });
            }
            // Check that f is only applied to slots in the block
            for i in 0..64 {
                assert!(!runs.contains_key(&i))
            }
            for i in 128..filter.nslots {
                assert!(!runs.contains_key(&i))
            }
        }
        #[test]
        fn test_last_intersecting_run_overshooting() {
            // Multiple runs in block 0 and 1,
            // but all runs in block i have runends in block i+1
            let mut filter = RSQF::new(64*3, 4);
            // (1) Run from block 0 to 1 [10: [10,70]]
            filter.set_occupied(10, true);
            filter.set_runend(70, true);
            filter.blocks[1].offset = 6;
            // (2) Run from block 0 to 1 [20: [71,80]
            filter.set_occupied(20, true);
            filter.set_runend(80, true);
            filter.blocks[1].offset = 16;
            // (3) Run from block 0 to block 2 [30:[81,130]]
            filter.set_occupied(30, true);
            filter.set_runend(130, true);
            filter.blocks[1].offset = 66;
            filter.blocks[2].offset = 2;
            // (4) Run from block 1 to block 2 [65: [131]]
            filter.set_occupied(65, true);
            filter.set_runend(131, true);
            filter.blocks[2].offset = 3;
            // (5) Run from block 1 to block 2 [120: [132]]
            filter.set_occupied(120, true);
            filter.set_runend(132, true);
            filter.blocks[2].offset = 4;

            assert_eq!(
                filter.last_intersecting_run(0),
                Some((10, 70)),
                "block[0]={:#?}",
                filter.blocks[0],
            );
            assert_eq!(
                filter.last_intersecting_run(1),
                Some((30, 130)),
                "block[0]={:#?}\nblock[1]={:#?}",
                filter.blocks[0], filter.blocks[1],
            );
            assert_eq!(
                filter.last_intersecting_run(2),
                Some((120, 132)),
                "block[1]={:#?}\nblock[2]={:#?}",
                filter.blocks[1], filter.blocks[2],
            );
        }
        #[test]
        fn test_multiblock_select_single_block() {
            let mut filter = RSQF::new(64, 4);
            let mut b;

            // Empty filter
            {
                b = &mut filter.blocks[0];
                assert_eq!(filter.select_runend(0, 0), None);
            }
            // Filter with one run
            {
                b = &mut filter.blocks[0];
                b.occupieds = 1;
                b.runends = 1;
                b.offset = 0;
                assert_eq!(filter.select_runend(0, 0), Some(0));
            }
            // Filter with multiple runs
            {
                b = &mut filter.blocks[0];
                // First two runs each have two elts, third run has one elt
                b.occupieds = 0b01101;
                b.runends   = 0b11010;
                b.offset    = 1;
                assert_eq!(filter.select_runend(0, 0), Some(1));
                assert_eq!(filter.select_runend(0, 1), Some(3));
                assert_eq!(filter.select_runend(0, 2), Some(4));
                assert_eq!(filter.select_runend(0, 3), None);
            }
        }
        #[test]
        fn test_multiblock_select_multiple_blocks() {
            let mut filter = RSQF::new(64*3, 4);

            // Filter with run starting in block 0 and ending in block 1
            {
                let b0 = &mut filter.blocks[0];
                b0.set_occupied(0, true);

                let b1 = &mut filter.blocks[1];
                b1.set_runend(0, true);
                b1.offset = 0;

                assert_eq!(filter.select_runend(0, 0), Some(64));
                assert_eq!(filter.select_runend(1, 0), Some(64));
                assert_eq!(filter.select_runend(2, 0), None);
            }

            // Filter with two runs:
            // A run starts in b0 and ends in b1, making b1 have nonzero offset
            {
                let b0 = &mut filter.blocks[0];
                b0.occupieds = 0b11;
                b0.runends   = 0b01;
                b0.offset    = 0;

                let b1 = &mut filter.blocks[1];
                b1.occupieds = 0b10;
                b1.runends   = 0b11;
                b1.offset    = 0;

                assert_eq!(filter.select_runend(0, 0), Some(0));
                assert_eq!(filter.select_runend(0, 1), Some(64));
                assert_eq!(filter.select_runend(1, 0), Some(64));
                assert_eq!(filter.select_runend(1, 1), Some(65));
            }
        }
        #[test]
        fn test_rank_select_single_block() {
            let mut filter = RSQF::new(64, 4);

            // Empty filter
            {
                for i in 0..64 {
                    assert_eq!(filter.rank_select(i), RankSelectResult::Empty)
                }
            }
            // Filter with one singleton run
            {
                for i in 0..64 {
                    let b = &mut filter.blocks[0];
                    b.set_occupied(i, true);
                    b.set_runend(i, true);
                    assert_eq!(filter.rank_select(i), RankSelectResult::Full(i));
                    for j in (i+1)..64 {
                        assert_eq!(filter.rank_select(j), RankSelectResult::Empty,
                                   "i={}, j={}", i, j);
                    }
                }

            }
            // Filter with two runs
            {
                let b = &mut filter.blocks[0];
                // First two runs each have two elts, third run has one elt
                b.occupieds = 0b101001;
                b.runends   = 0b110010;
                b.offset    = 1;
                assert_eq!(filter.rank_select(0), RankSelectResult::Full(1));
                assert_eq!(filter.rank_select(1), RankSelectResult::Full(1));
                assert_eq!(filter.rank_select(2), RankSelectResult::Empty);
                assert_eq!(filter.rank_select(3), RankSelectResult::Full(4));
                assert_eq!(filter.rank_select(4), RankSelectResult::Full(4));
                assert_eq!(filter.rank_select(5), RankSelectResult::Full(5));
                for i in 6..64 {
                    assert_eq!(filter.rank_select(i), RankSelectResult::Empty);
                }
            }
        }
        #[test]
        fn test_rank_select_multi_block_1() {
            // Filter with run starting in block 0 and ending in block 1
            let mut filter = RSQF::new(64*3, 4);
            let b0 = &mut filter.blocks[0];
            b0.set_occupied(0, true);
            b0.offset = 64;

            let b1 = &mut filter.blocks[1];
            b1.set_runend(0, true);
            b1.offset = 0;

            for i in 0..=64 {
                assert_eq!(filter.rank_select(i), RankSelectResult::Full(64), "i={}", i);
            }
            for i in 65..filter.nslots {
                assert_eq!(filter.rank_select(i), RankSelectResult::Empty, "i={}", i);
            }
        }
        #[test]
        fn test_rank_select_multi_block_2() {
            // Run 1: [0, [0]]
            // Run 2: [1, [1,64]]
            // Run 3: [65, [65, 68]]
            // Run 4: [66, [69, 130]]
            let mut filter = RSQF::new(64*3, 4);
            let b0 = &mut filter.blocks[0];
            b0.set_occupied(0, true); // start run 1
            b0.set_runend(0, true);   // end run 1
            b0.offset = 0;
            b0.set_occupied(1, true); // start run 2

            let b1 = &mut filter.blocks[1];
            b1.set_runend(0, true);   // end run 2
            b1.offset = 0;
            b1.set_occupied(1, true); // start run 3
            b1.set_runend(4, true);   // end run 3
            b1.set_occupied(2, true); // start run 4

            let b2 = &mut filter.blocks[2];
            b2.set_runend(2, true);   // end of run 4
            b2.offset = 2;

            assert_eq!(filter.rank_select(0), RankSelectResult::Full(0));
            for i in 1..=64 {
                assert_eq!(filter.rank_select(i),
                           RankSelectResult::Full(64),
                           "i={}", i);
            }
            assert_eq!(filter.rank_select(65),
                       RankSelectResult::Full(68));
            for i in 66..=130 {
                assert_eq!(filter.rank_select(i),
                           RankSelectResult::Full(130),
                           "i={}", i);
            }
            for i in 131..filter.nslots {
                assert_eq!(filter.rank_select(i),
                           RankSelectResult::Empty,
                           "i={}", i);
            }
        }
        #[test]
        fn test_first_unused_empty() {
            // Empty filter
            // The first unused slot for a slot x should be x itself
            let filter = RSQF::new(128, 4);
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
                let mut filter = RSQF::new(nslots, 4);
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
        fn insert_run(filter: &mut RSQF, a: usize, b: usize) {
            debug_assert!(a < 128 && b < 128);

            // Setup filter
            filter.set_occupied(a, true);
            filter.set_runend(b, true);
            // Set offsets
            if a == 0 {
                filter.blocks[0].offset = b;
            }
            if a == 64 || (a < 64 && b >= 64) {
                filter.blocks[1].offset = b-64;
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
                    let mut filter = RSQF::new(nslots, 4);
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
                debug_assert!(a < b && b < c && c < d);

                let mut filter = RSQF::new(nslots, 4);
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
                        "r1=[{}, {}], r2=[{}, {}], i={}, filter={:#?}",
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
        fn test_query_empty_single() {
            let mut filter = RSQF::new(64, 4);
            for i in 0..64 {
                let word = i.to_string();
                assert!(!filter.query(word.clone()), "word='{}'", word);
            }
        }
        #[test]
        /// 1-elt single-block filter
        fn test_query_one_single() {
            let mut filter = RSQF::new(64, 4);
            let word = "apples";

            // Manually insert 'apples'
            let hash = filter.hash(word);
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash);
            let b = &mut filter.blocks[0];
            b.remainders[quot] = rem;
            b.set_occupied(quot, true);
            b.set_runend(quot, true);

            assert!(filter.query("apples".to_string()));
        }
        #[test]
        /// Multi-elt single-block filter
        fn test_query_multi_single() {
            let mut filter = RSQF::new(64, 4);
            // I checked that these words' quot/rems don't conflict
            let words = ["apples", "bananas", "oranges"];

            // Manually insert words
            for word in words.iter() {
                let hash = filter.hash(word);
                let quot = filter.calc_quot(hash);
                let rem = filter.calc_rem(hash);
                // println!("word={}, quot={:x}, rem={:x}", word, quot, rem);

                let b = &mut filter.blocks[0];
                b.remainders[quot] = rem;
                b.set_occupied(quot, true);
                b.set_runend(quot, true);
            }
            // Check that words are contained
            for word in words.iter() {
                assert!(filter.query(word.to_string()));
            }
        }
        #[test]
        // Empty multi-block filter
        fn test_query_empty_multi() {
            let mut filter = RSQF::new(64*3, 4);
            for i in 0..filter.nslots {
                let word = i.to_string();
                assert!(!filter.query(word.to_string()), "word={}", word);
            }
        }
        #[test]
        // 1-elt multi-block filter
        fn test_query_one_multi() {
            let mut filter = RSQF::new(64*3, 4);
            let word = "apples";

            // Manually insert 'apples'
            let hash = filter.hash(word);
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash);
            let b = &mut filter.blocks[quot/64];
            b.remainders[quot%64] = rem;
            b.set_occupied(quot%64, true);
            b.set_runend(quot%64, true);

            assert!(filter.query("apples".to_string()));
        }
        #[test]
        // Multi-elt multi-block filter
        fn test_query_multi_multi() {
            let mut filter = RSQF::new(64*3, 4);
            let words = ["apples", "bananas", "oranges"];

            // Manually insert words
            for word in words.iter() {
                let hash = filter.hash(word);
                let quot = filter.calc_quot(hash);
                let rem = filter.calc_rem(hash);
                // println!("word={}, quot={:x}, rem={:x}", word, quot, rem);

                let b = &mut filter.blocks[quot/64];
                b.remainders[quot%64] = rem;
                b.set_occupied(quot%64, true);
                b.set_runend(quot%64, true);
            }
            // Check that words are contained
            for word in words.iter() {
                assert!(filter.query(word.to_string()));
            }
        }
        #[test]
        fn test_shift_rem_and_runends() {
            let mut filter = RSQF::new(128, 4);
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
        fn offset_state_init() -> RSQF {
            let mut filter = RSQF::new(64*4, 4);
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
        fn test_inc_offset_negative_target() {
            // Check that target doesn't go negative (block_i == 0 and b[0][0] unoccupied)
            let mut filter = RSQF::new(64, 4); // target is indirect
            filter.inc_offsets(0,0);
            assert_eq!(filter.blocks[0].offset, 0);
            filter.blocks[0].set_occupied(0, true); // target is now direct
            filter.inc_offsets(0,0);                // should bump up offset to 1
            assert_eq!(filter.blocks[0].offset, 1);
            filter.inc_offsets(1,1);                // should bump up offset to 2
            assert_eq!(filter.blocks[0].offset, 2);

            // Check for multi-block case
            let mut filter = RSQF::new(64*5, 4);
            filter.inc_offsets(0, filter.nslots-1);
            assert_eq!(filter.blocks[0].offset, 0);
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
        #[test]
        fn test_inc_indirect_offset_negative_target() {
            // Shouldn't ever affect first block b/c first block
            // can't have an indirect target
            let mut filter = RSQF::new(64, 4);
            filter.inc_offsets(0, 0);
            assert_eq!(filter.blocks[0].offset, 0);

            // Check for case where indirect offset
            // points to run in previous block at slot 63
            // This also tests our check for block_i == 0 in the loop
            let mut filter = RSQF::new(128, 4);
            filter.blocks[0].set_occupied(63, true);
            filter.blocks[0].set_runend(63, true);
            filter.inc_indirect_offsets(0, 63);
            assert_eq!(filter.blocks[1].offset, 1);
        }
        // Checks if x is a power of 2
        fn is_pow_2(x: usize) -> bool {
            let mask = if x == 0 { 0 } else { x-1 };
            return (x & mask) == 0;
        }
        #[test]
        fn test_raw_insert_new_run() {
            // Insert new runs that don't overlap with anything
            let mut filter = RSQF::new(128, 4);
            // Insert new runs at powers of 2
            for i in 1..filter.nslots {
                if is_pow_2(i) {
                    filter.raw_insert(i, i as Rem);
                }
            }
            // Check offsets
            assert_eq!(filter.blocks[0].offset, 0); // indirect offset b/c b[0][0] empty
            assert_eq!(filter.blocks[1].offset, 0); // direct offset b/c b[1][0] has 1-elt run
            // Check occupieds/runends
            for i in 1..filter.nslots {
                assert_eq!(filter.is_occupied(i), is_pow_2(i));
                assert_eq!(filter.is_runend(i), is_pow_2(i));
                assert_eq!(filter.remainder(i),
                           if is_pow_2(i) { i } else { 0 } as Rem);
            }
        }
        #[test]
        fn test_raw_insert_overlapping_run() {
            // Insert a new run that overlaps with an existing run
            // and shifts multiple offsets

            fn one_long_run() -> RSQF {
                // Setup: one existing run: [0:(0,130)]
                let mut filter = RSQF::new(64*3, 4);
                let b = &mut filter.blocks;
                b[0].set_occupied(0, true);
                b[0].offset = 130;      // direct offset
                b[1].offset = 130-64+1; // indirect offset
                b[2].offset = 3;        // indirect offset
                b[2].set_runend(2, true);
                for i in 0..=130 {
                    filter.set_remainder(i, i as Rem);
                }
                filter
            }
            // Insert after the run ends
            let mut filter = one_long_run();
            filter.raw_insert(131, 0xff);
            for i in 0..filter.nslots {
                assert_eq!(filter.is_occupied(i), i==0 || i==131, "i={}", i);
                assert_eq!(filter.is_runend(i), i==130 || i==131, "i={}", i);
                assert_eq!(filter.remainder(i),
                           if i < 131 { i as Rem }
                           else if i == 131 { 0xff }
                           else { 0 },
                           "i={}",
                           i);
            }
            // Insert into a new run with quot where 0 < quot < 130
            let mut filter = one_long_run();
            filter.raw_insert(10, 0xff);
            for i in 0..filter.nslots {
                assert_eq!(filter.is_occupied(i), i==0 || i==10, "i={}", i);
                assert_eq!(filter.is_runend(i), i==130 || i==131, "i={}", i);
                assert_eq!(filter.remainder(i),
                           if i < 131 { i as Rem }
                           else if i == 131 { 0xff }
                           else { 0 },
                           "i={}",
                           i);
            }
            // Extend the run (insert with quot=0)
            let mut filter = one_long_run();
            filter.raw_insert(0, 131);
            for i in 0..filter.nslots {
                assert_eq!(filter.is_occupied(i), i==0, "i={}", i);
                assert_eq!(filter.is_runend(i), i==131, "i={}", i);
                assert_eq!(filter.remainder(i),
                           if i <= 131 { i as Rem }
                           else { 0 },
                           "i={}",
                           i);
            }
        }
        #[test]
        fn test_raw_insert_extend() {
            let mut filter = RSQF::new(128, 4);
            for i in 0..filter.nslots {
                filter.set_occupied(i, true);
                filter.set_runend(i, true);
                filter.set_remainder(i, i as Rem);
            }
            filter.raw_insert(0, 0xff);
            for i in 0..filter.nslots {
                assert_eq!(filter.is_occupied(i), i < 128, "i={}", i);
                assert_eq!(filter.is_runend(i), i > 0 && i <= 128, "i={}", i);
                assert_eq!(filter.remainder(i),
                           if i==0 { 0 }
                           else if i==1 { 0xff }
                           else if i <= 128 { (i-1) as Rem }
                           else { 0 },
                           "i={}",
                           i);
            }
        }
        #[test]
        fn test_raw_insert_repeated() {
            let n = 1 << 10;
            let mut filter = RSQF::new(n, 4);
            for i in 0..n {
                filter.insert("apple".to_string());
                assert!(filter.query("apple".to_string()));
                // if i % 1000 == 0 {
                //     println!("i={}", i);
                // }
            }
            println!("filter={:#?}", filter);
        }
        #[test]
        fn test_insert_and_query() {
            // Insert and query elts, ensure that there are no false negatives
            let a: usize = 1 << 20;
            let ratio = 100.0;      // a/s
            let s = nearest_pow_of_2((a as f64/ratio) as usize);
            let s = ((s as f64) * 0.95) as usize;
            let mut filter = RSQF::new(s, 4);
            // Generate query set
            let mut rng = SmallRng::seed_from_u64(0);
            let mut set = HashSet::with_capacity(s);
            for _i in 0..s {
                let elt = (rng.gen::<usize>() % a).to_string();
                set.insert(elt.clone());
                filter.insert(elt.clone());
            }
            println!("|set|: {}, s={}, a={}, load={}",
                     set.len(), s, a, filter.load());
            // Query [0, a] and ensure that all items in the set return true
            let mut fps = 0;
            for i in 0..a {
                let elt = &i.to_string();
                if set.contains(elt) {
                    assert!(filter.query(elt.clone()), "elt={}", elt);
                } else {
                    fps += filter.query(elt.clone()) as usize;
                }
            }
            println!("FP rate: {}", (fps as f64)/(a as f64));
        }
        #[test]
        fn test_load() {
            // Insert some elts in filter and check load factor
            let mut filter = RSQF::new(1000, 4);
            for i in 0..500 {
                filter.insert(i.to_string());
            }
            assert_eq!(filter.load(), 500_f64/(filter.nslots as f64));
        }
    }
}

