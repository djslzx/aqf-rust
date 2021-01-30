use std::{
    cmp::max,
    collections::HashMap,
    fmt,
};
use crate::{Rem, Filter};
use crate::util::{
    bitarr::{b64, b128},
    join_displayable,
    nearest_pow_of_2,
};
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

impl fmt::Debug for Block {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        RankSelectBlock::fmt(self, f)
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
        let nblocks = max(1, nearest_pow_of_2(n)/64);
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
    fn block(&self, i: usize) -> &Block {
        debug_assert!(i < self.nblocks, "Block index out of bounds: ({}/{})", i, self.nblocks);
        &self.blocks[i]
    }
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
    fn adapt(&mut self, loc: usize, quot: usize, rem: Rem,
             member_hash: u128, non_member_hash: u128, exts: [Ext; 64]) {
        let new_ext = self.shortest_diff_ext(member_hash, non_member_hash);
        assert_ne!(new_ext, Ext::None,
                   "Hashes were identical, member_hash={}, non_member_hash={}",
                   member_hash, non_member_hash);
        // Write encoding to the appropriate block
        let mut exts = exts;
        //println!("exts={:?}", join_displayable(&exts, ","));
        let old_ext = exts[loc%64];
        exts[loc%64] = new_ext;
        match ExtensionArcd::encode(exts) {
            Ok(code) => {
                // Update code and add new extension to remote
                self.blocks[loc/64].extensions = code;
                self.update_remote_ext(quot, rem, old_ext, new_ext);
            }
            Err(_) => {
                // Encoding failed: rebuild
                // Clear all extensions in the offending block + in remote rep
                self.clear_block_remote_exts(loc/64);
                // Write new extension encoding
                exts = [Ext::None; 64];
                exts[loc%64] = new_ext;
                match ExtensionArcd::encode(exts) {
                    Ok(code) => {
                        self.blocks[loc/64].extensions = code;
                        // Add new extension back into the remote
                        self.update_remote_ext(quot, rem, Ext::None, new_ext);
                    }
                    Err(_) =>
                        panic!("Failed to encode after rebuilding block: \
                                block={:?}, new_ext={:?}",
                               self.blocks[loc/64], new_ext),
                }
            }
        }
    }
    /// Rebuild a block's extensions in the remote representation.
    #[allow(unused_variables)]
    fn clear_block_remote_exts(&mut self, block_i: usize) {
        let exts = ExtensionArcd::decode(self.blocks[block_i].extensions);
        // Clear extensions from remote representation
        self.apply_to_block(
            block_i,
            |filter: &mut Self, quot: usize, i: usize| {
                let rem = filter.remainder(i);
                let ext = exts[i%64];
                filter.clear_remote_ext(quot, rem, ext);
            });
    }
    /// Removes the fingerprint extension from the remote rep
    fn clear_remote_ext(&mut self, quot: usize, rem: Rem, ext: Ext) {
        self.update_remote_ext(quot, rem, ext, Ext::None);
    }
    /// Updates a fingerprint extension in the remote rep
    fn update_remote_ext(&mut self, quot: usize, rem: Rem, old_ext: Ext, new_ext: Ext) {
        // Remove previous (quot, rem, ext) triple
        let val = match self.remote.remove(&(quot, rem, old_ext)) {
            Some(v) => v,
            None => panic!("Tried to update a triple that doesn't exist: \
                            quot={}, rem={}, ext={}, new_ext={}, filter: {:#?}",
                           quot, rem, old_ext, new_ext, self)
        };
        // Insert new (quot, rem, empty_ext) triple
        self.remote.insert((quot, rem, new_ext), val);
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
                        let ext = exts[loc%64];
                        match ext {
                            // If extensions exist and don't match, move on
                            Ext::Some {bits, len} if self.calc_ext(query_hash, len) != bits => {}
                            // Otherwise, extension doesn't exist (automatic match) or exists and matches
                            _ => {
                                // Check remote to see if true match
                                // NOTE: think about only storing the hash in remote and ditching
                                // the word (comparing word to elt is potentially slow)
                                if let Some((word, remote_hash)) = self.remote.get(&(quot, rem, ext)) {
                                    let hash = *remote_hash; // finish immut borrow before mut borrow for adapt
                                    if *word != elt {
                                        // false match => adapt
                                        self.adapt(loc, quot, rem, hash, query_hash, exts);
                                    }
                                }
                                return true;
                            }
                        }
                    }
                    // Stop when l < 0, l = quot, or l-1 is a runend
                    if loc == 0 || loc == quot || self.is_runend(loc-1) {
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
        // Don't insert duplicates FIXME
        if !self.remote.contains_key(&(quot, rem, Ext::None)) {
            self.raw_insert(quot, rem);
            self.remote.insert((quot, rem, Ext::None), (elt, hash));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use rand::{
        Rng, SeedableRng, rngs::SmallRng
    };

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
    fn test_clear_block_remote_exts_empty() {
        // Empty extensions -> no change
        let mut filter = AQF::new(64*4, 4);
        for i in 0..filter.nblocks {
            filter.clear_block_remote_exts(i);
        }
        // Check remote rep
        assert!(filter.remote.is_empty(),
                "remote is populated by rebuild");
        // Check extensions
        for i in 0..filter.nblocks {
            assert_eq!(filter.blocks[i].extensions, 0,
                       "block extensions are nonzero after rebuild");
        }
    }
    #[test]
    fn test_clear_block_remote_exts() {
        // Insert elts with conflicting hashes raw
        let nslots = 64*2;
        let nelts = 10;
        let mut filter = AQF::new(nslots, 4);
        for i in 0..nelts {
            let elt = i.to_string();
            let hash = filter.hash(&elt);
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash);
            filter.raw_insert(quot, rem);
            filter.remote.insert((quot, rem, Ext::Some{bits: 1, len: 1}), (elt, hash));
        }
        // Edit codes
        for b in filter.blocks.iter_mut() {
            let mut raw = [Ext::None; 64];
            for j in 0..64 {
                raw[j] = if b.is_occupied(j) {
                    Ext::Some{bits: 1, len: 1}
                } else {
                    Ext::None
                };
            }
            let code = match ExtensionArcd::encode(raw) {
                Ok(code) => code,
                Err(_) => panic!("Failed to encode extensions={:?}", raw),
            };
            b.extensions = code;
            println!("exts={:?}, code={}", join_displayable(&raw, ""), code);
        }
        // Clear and check remote
        for i in 0..filter.nblocks {
            filter.clear_block_remote_exts(i);
        }
        for (_, _, ext) in filter.remote.keys() {
            assert_eq!(*ext, Ext::None);
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
                assert!(filter.query(elt.clone()), "False negative on elt={}", elt);
            } else {
                fps += filter.query(elt.clone()) as usize;
            }
        }
        println!("FP rate: {}", (fps as f64)/(a as f64));
    }
}
