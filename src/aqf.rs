use std::{cmp, collections::{HashSet, HashMap}};
use rand::{
    Rng, SeedableRng, rngs::SmallRng
};
use crate::{Rem, Filter};
use crate::util::{bitarr::{b64, b128}, nearest_pow_of_2, };
use crate::rsqf::{
    RankSelectBlock,
    RankSelectResult,
    RankSelectQuotientFilter,
};
use crate::arcd::{
    Arcd,
    Ext,
    ext_arcd::ExtensionArcd,
    // selector_arcd::SelectorArcd,
};

#[derive(Debug)]
struct Block {
    remainders: [Rem; 64],
    occupieds: u64,
    runends: u64,
    extensions: u64,
    offset: usize,
}

impl RankSelectBlock for Block {
    fn new() -> Block {
        Block {
            remainders: [0; 64],
            occupieds: 0,
            runends: 0,
            extensions: 0,
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

impl Block {
    // TODO: extension encoding/decoding
}

#[derive(Debug)]
struct AQF {
    blocks: Vec<Block>,
    nblocks: usize,
    nslots: usize,
    nelts: usize,

    // Sizes
    p: usize, // Fingerprint size
    q: usize, // Quotient size
    r: usize, // Remainder size

    seed: u32,

    // Remote representation
    // (quot, rem, ext) -> (elt, hash)
    remote: HashMap<(usize, Rem, Ext), (String, u128)>,
}

impl RankSelectQuotientFilter for AQF {
    type Block = Block;

    fn new(n: usize, r: usize) -> AQF {
        Self::new_seeded(n, r, 0)
    }
    fn new_seeded(n: usize, r: usize, seed: u32) -> AQF {
        let nblocks = cmp::max(1, nearest_pow_of_2(n)/64);
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
            nelts: 0,
            q,
            r,
            p: q + r,
            seed,
            remote: HashMap::new(),
        }
    }
    fn q(&self) -> usize { self.q }
    fn r(&self) -> usize { self.r }
    fn nslots(&self) -> usize { self.nslots }
    fn nblocks(&self) -> usize { self.nblocks }
    fn seed(&self) -> u32 { self.seed }
    fn block(&self, i: usize) -> &Block { &self.blocks[i] }
    fn mut_block(&mut self, i: usize) -> &mut Block { &mut self.blocks[i] }
    fn add_block(&mut self) {
        let b = Block::new();
        self.blocks.push(b);
        self.nslots += 64;
        self.nblocks += 1;
    }
}

impl AQF {
    fn calc_ext(&self, hash: u128, k: usize) -> u64 {
        debug_assert_ne!(k, 0);
        let a = self.q + self.r;
        let b = a + k;
        assert!(b <= 64, "Extension overflowed 64 bits, b={}", b);
        ((hash & b128::half_open(a, b)) >> a) as u64
    }
    /// Generate the shortest extension from the member's hash that doesn't conflict with
    /// the non-member's hash. 
    /// `prev_ext` is the extension previously associated with the member's hash.
    fn shortest_diff_ext(&self, member_hash: u128, non_member_hash: u128) -> Ext {
        let a = member_hash >> (self.q + self.r);
        let b = non_member_hash >> (self.q + self.r);
        if a == b {
            Ext::None
        } else {
            // Find fewest LSBs needed to distinguish member from non-member hash:
            // Determine number of common LSBs and add 1
            let len = (a ^ b).trailing_zeros() + 1;
            let bits = a & ((1 << len) - 1);
            Ext::Some {
                bits: bits as u64,
                len: len as usize,
            }
        }
    }
    /// Adapt on a false match for a fingerprint at loc
    fn adapt(&mut self, loc: usize, member_hash: u128, non_member_hash: u128, mut exts: [Ext; 64]) {
        let new_ext = self.shortest_diff_ext(member_hash, non_member_hash);
        assert_ne!(new_ext, Ext::None,
                   "Hashes were identical, member_hash={}, non_member_hash={}",
                   member_hash, non_member_hash);
        // Write encoding to the appropriate block
        exts[loc%64] = new_ext;
        match ExtensionArcd::encode(exts) {
            Ok(code) =>
                self.blocks[loc/64].extensions = code,
            Err(_) => {
                // Encoding failed: rebuild (TODO)
                unimplemented!("Need to rebuild for ext={:?}, exts={:?}", new_ext, exts);
            }
        }
    }
    /// Rebuild a block's extensions, determining the block's bounds heuristically.
    /// This function updates extensions & the remote rep.
    fn rebuild_block(&mut self, block_i: usize) {
        // Edit remote representation:
        // Go through all (quot, rem) pairs in this block and shorten their remainders
        let block_start = block_i * 64;
        // Find the quot of the last runend before this block
        let mut q = self.last_prior_run(block_i);
        let exts = ExtensionArcd::decode(self.blocks[block_i].extensions);

        // Go through all the runs in this block and clear their extensions in the remote rep
        'block: loop {
            if self.is_occupied(q) {
                // Get the end of the run associated with q
                let mut i = match self.rank_select(q) {
                    RankSelectResult::Full(loc) => loc,
                    RankSelectResult::Empty => panic!("Occupied slot should not be empty"),
                    RankSelectResult::Overflow => panic!("Occupied slot should not go off the edge"),
                };
                // Set the extensions in the run to the empty extension
                'run: loop {
                    // Clear ext if in current block
                    // (The run that i is going through might end in the next block
                    // but still go through the block we're interested in)
                    if i/64 == block_i {
                        let rem = self.remainder(i);
                        self.clear_ext_from_remote(q, rem, exts[i%64]);
                    }
                    // If i reaches the block's start, then exit the outer loop.
                    // If i reaches the start of the run (i=quot or i-1 is a runend),
                    // then exit the inner loop. Otherwise, walk backwards through the run.
                    if i <= block_start {
                        break 'block;
                    } else if i <= q || self.is_runend(i-1) {
                        break 'run;
                    } else {
                        i -= 1;
                    }
                }
            }
            q += 1;
            // TODO: advance immediately to the next occupied quotient:
            // q = select(rank(q)+1)
        }
        // Set all block extensions to 0
        self.blocks[block_i].extensions = 0;
    }
    /// Removes the fingerprint extension from the remote rep
    fn clear_ext_from_remote(&mut self, quot: usize, rem: Rem, ext: Ext) {
        // Remove previous (quot, rem, ext) triple
        let val = self.remote.remove(&(quot, rem, ext)).unwrap();
        // Insert new (quot, rem, empty_ext) triple
        self.remote.insert((quot, rem, Ext::None), val);
    }
    /// Finds the quotient and runend of the last run whose runend is before this block
    fn last_prior_run(&self, block_i: usize) -> usize {
        let block_start = block_i * 64;
        let mut q = block_start;
        loop {
            match self.rank_select(q) {
                RankSelectResult::Full(loc) =>
                // Exit when the end of q's run is at or before the start of the block
                // or if q is the very first quotient; otherwise, step backwards
                    if loc <= block_start || q == 0 {
                        break q;
                    } else {
                        q -= 1;
                    }
                RankSelectResult::Empty =>
                // Exit when the home slot for q is free and q is at or before
                // the start of the block; otherwise, step backwards
                    if q <= block_start {
                        break q;
                    } else {
                        q -= 1;
                    }
                RankSelectResult::Overflow => panic!("Rebuilding went off the edge"),
                // We should never hit this case: if we do, that means that we don't
                // have a 1-1 between occupieds and runends
            }
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
                // TODO: also shift fingerprint extensions
                self.inc_offsets(r+1, u-1);
                // Start a new run or extend an existing one
                if !self.is_occupied(quot) {
                    // Set occupied, add runend, add rem, shift indirect offsets
                    self.set_occupied(quot, true);
                    self.set_runend(r+1, true);
                    self.set_remainder(r+1, rem);
                    // TODO: add fingerprint extension
                    self.inc_indirect_offsets(quot, r);
                } else {
                    // Don't need to set occupied
                    // Shift runend, add rem, shift offsets
                    self.set_runend(r, false);
                    self.set_runend(r+1, true);
                    self.set_remainder(r+1, rem);
                    // TODO: add fingerprint extension
                    self.inc_offsets(r, r);
                }
            }
            RankSelectResult::Overflow =>
                panic!(
                    "AQF failed to find runend (nslots={}, quot=(block={}, slot={}))",
                    self.nslots, quot/64, quot%64,
                ),
        }
    }
    /// Computes filter load factor
    fn load(&self) -> f64 {
        (self.nelts as f64)/(self.nslots as f64)
    }
}

impl Filter<String> for AQF {
    fn query(&mut self, elt: String) -> bool {
        let query_hash = self.hash(&elt[..]);
        let quot = self.calc_quot(query_hash);
        let rem = self.calc_rem(query_hash);

        if !self.is_occupied(quot) {
            false
        } else {
            if let RankSelectResult::Full(mut loc) = self.rank_select(quot) {
                // Cache decoded letters as (block_i, letters), using block_i
                // to check whether we need to update letters
                let mut decode: Option<(usize, [Ext; 64])> = None;
                loop {
                    // Matching remainder found => compare extensions
                    if self.remainder(loc) == rem {
                        // Check cached code:
                        match decode {
                            // If cached block index is the same as the current index,
                            // leave cache as is
                            Some((block_i, _)) if block_i == loc/64 => {}
                            // Otherwise, decode and store result in cache
                            _ => {
                                let exts = self.blocks[loc/64].extensions;
                                decode = Some((loc/64, ExtensionArcd::decode(exts)));
                            }
                        }
                        // Check if extensions match:
                        // We should be able to unwrap cached decode w/o error b/c of the previous match
                        let exts = decode.unwrap().1;
                        let ext = exts[loc/64];
                        match ext {
                            // If extensions exist and don't match, move on
                            Ext::Some {bits, len} if self.calc_ext(query_hash, len) != bits => {}
                            // Otherwise, extension doesn't exist (automatic match) or exists and matches
                            _ => {
                                // Check remote to see if true match
                                // NOTE: think about only storing the hash in remote and ditching
                                // the word (comparing word to elt is potentially slow)
                                if let Some((word, remote_hash)) = self.remote.get(&(quot, rem, ext)) {
                                    if *word != elt {
                                        // false match => adapt
                                        self.adapt(loc, *remote_hash, query_hash, exts);
                                    }
                                }
                                return true;
                            }
                        }
                    }
                    // Stop when l < 0, l-1 < quot, or l-1 is a runend
                    if loc == 0 || loc-1 < quot || self.is_runend(loc-1) {
                        return false;
                    } else {
                        loc -= 1;
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
        let rem = self.calc_rem(hash);
        self.raw_insert(quot, rem);
        self.remote.insert((quot,rem,Ext::None), (elt, hash));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_extension() {
        let filter = AQF::new(128, 4);
        let hash = filter.hash("apple");
        let start = filter.q + filter.r;
        for i in 1..(64 - start) {
            assert_eq!(
                filter.calc_ext(hash, i),
                ((hash & b128::half_open(start,start+i)) >> start) as u64,
            );
        }
    }
    #[test]
    fn test_shortest_diff_ext() {
        // Identical hashes:
        let filter = AQF::new(64, 4);
        assert_eq!(filter.shortest_diff_ext(0, 0), Ext::None);
        for i in 1..64 {
            // h1, h2 only differ at position j
            let h1 = 0;
            let h2 = 1 << (i + filter.r + filter.q);
            assert_eq!(
                filter.shortest_diff_ext(h1, h2),
                Ext::Some {
                    bits: 0,
                    len: i+1,
                }
            );
            assert_eq!(
                filter.shortest_diff_ext(h2, h1),
                Ext::Some {
                    bits: 1 << i,
                    len: i+1,
                }
            );
        }
    }
    #[test]
    fn test_insert_and_query() {
        // Insert and query elts, ensure that there are no false negatives
        let a: usize = 1 << 14; // use set size 2^14
        let ratio = 100.0;      // a/s
        let s = nearest_pow_of_2((a as f64/ratio) as usize);
        let s = ((s as f64) * 0.95) as usize;
        let mut filter = AQF::new(s, 4);
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
}
