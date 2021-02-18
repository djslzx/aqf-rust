use std::{
    cmp::max,
    collections::HashMap,
    fmt,
};
#[macro_use]
use crate::{
    hashmap,
    Rem,
    Filter,
    util::{
        bitarr::{b64, b128},
        nearest_pow_of_2,
    },
    rsqf::{
        RankSelectBlock,
        RankSelectResult,
        RankSelectQuotientFilter,
    },
    arcd::{
        Arcd,
        extensions::Ext,
        ext_arcd::ExtensionArcd,
    },
};

pub struct Block {
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
        let decoded = ExtensionArcd::decode(self.extensions);
        let exts = decoded
            .iter()
            .enumerate()
            .filter_map(|(i, ext)|
                match ext {
                    Ext::Some{bits, len} => Some(format!("{a}: {b:0c$b}", a=i, b=bits, c=len)),
                    Ext::None => None
                })
            .collect::<Vec<String>>()
            .join(", ");
        let rems = (0..64)
            .map(|x| format!("0x{:x}", self.remainder(x)))
            .collect::<Vec<String>>()
            .chunks(8)
            .map(|c| c.join(" "))
            .collect::<Vec<String>>();
        f.debug_struct("Block")
            .field("remainders", &rems)
            .field("occupieds ", &format_args!("[{:064b}]", self.occupieds().reverse_bits()))
            .field("runends   ", &format_args!("[{:064b}]", self.runends().reverse_bits()))
            .field("offset    ", &self.offset())
            .field("extensions", &format_args!("{}", &exts))
            .field("code      ", &self.extensions)
            .finish()
    }
}

impl Block {
    // TODO: extension encoding/decoding
}

mod remote {
    use super::*;

    type Quot = usize;

    pub struct Remote<T> {
        pub data: HashMap<(Quot, Rem), Vec<(Ext, T, u128)>>,
    }
    impl fmt::Debug for Remote<String> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.debug_map()
                .entries(self.data.iter().map(|(&(q, r), vec)| {(
                    format!("{}, {}", q, r),
                    vec.iter()
                        .map(|(e, s, h)| format!("{}, {}, {}", e, s, h))
                        .collect::<Vec<String>>()
                )}))
                .finish()
        }
    }
    impl Remote<String> {
        /// Constructs a new remote representation
        pub fn new() -> Self {
            let data: HashMap<(Quot, Rem), Vec<(Ext, String, u128)>> = HashMap::new();
            Remote {
                data
            }
        }
        pub fn is_empty(&self) -> bool {
            self.data.is_empty()
        }
        pub fn items(&self) -> Vec<(&(Quot, Rem), &Vec<(Ext, String, u128)>)> {
            self.data.iter().collect()
        }
        /// Inserts the mapping `(quot, rem, ext) -> (elt, hash)` into the remote representation
        pub fn add(&mut self, quot: Quot, rem: Rem, ext: Ext, elt: String, hash: u128) {
            match self.data.get_mut(&(quot, rem)) {
                Some(v) => {
                    v.push((ext, elt, hash));
                }
                None => {
                    self.data.insert((quot, rem), vec![(ext, elt, hash)]);
                }
            }
        }
        /// Returns true if `(quot, rem, ext) -> (elt, hash)` is in the remote, false otherwise
        pub fn contains(&self, quot: Quot, rem: Rem, ext: Ext, elt: &str, hash: u128) -> bool {
            match self.data.get(&(quot, rem)) {
                None => false,
                Some(v) => v.iter().any(|(a, b, c)|
                    ext == *a && elt == b && hash == *c
                )
            }
        }
        /// Returns the `(elt, hash)` pairs associated with a `(quot, rem, ext)` triple
        pub fn get(&self, quot: Quot, rem: Rem, ext: Ext) -> Vec<(String, u128)> {
            match self.data.get(&(quot, rem)) {
                None => Vec::new(),
                Some(vec) => {
                    vec.iter()
                        .filter_map(|(e, s, h)|
                            if *e == ext {
                                Some((s.clone(), *h))
                            } else {
                                None
                            })
                        .collect()
                }
            }
        }
        /// Updates the extension for a `(quot, rem, old_ext)` triple to `new_ext`.
        ///
        /// Finds the first `(old_ext, word, hash)` triple associated with `quot` and `rem`
        /// and updates its extension.
        pub fn update_ext(&mut self, quot: Quot, rem: Rem, old_ext: Ext, new_ext: Ext) {
            match self.data.get_mut(&(quot, rem)) {
                None => panic!("Couldn't find any entries associated with quot={}, rem={}",
                               quot, rem),
                Some(vec) => {
                    match vec.iter_mut().find(|(ext, _, _)| *ext == old_ext) {
                        None => panic!("Couldn't find any entries associated with quot={}, rem={}, ext={}",
                                       quot, rem, old_ext),
                        Some((ext, _, _)) => *ext = new_ext,
                    }
                }
            }
        }
        pub fn clear_ext(&mut self, quot: Quot, rem: Rem, ext: Ext) {
            self.update_ext(quot, rem, ext, Ext::None);
        }
    }
    #[cfg(test)]
    mod tests {
        use super::*;

        fn make_ext(len: usize) -> Ext {
            if len == 0 {
                Ext::None
            } else {
                Ext::Some {
                    bits: 0,
                    len,
                }
            }
        }
        fn make_elt(len: usize) -> String {
            "a".repeat(len)
        }

        #[test]
        fn test_add() {
            let mut r = Remote::new();
            r.add(0, 0, Ext::None, "a".to_string(), 0);
            for i in 1..=3 {
                r.add(1, 1, Ext::Some{bits: 0, len: i}, "b".repeat(i), 1);
            }

            assert_eq!(r.data, hashmap!(
              (0,0) => vec![
                (Ext::None, "a".to_string(), 0),
              ],
              (1,1) => vec![
                (Ext::Some{bits: 0, len: 1}, "b".to_string(), 1),
                (Ext::Some{bits: 0, len: 2}, "bb".to_string(), 1),
                (Ext::Some{bits: 0, len: 3}, "bbb".to_string(), 1)
              ]
            ));
        }
        #[test]
        fn test_contains() {
            let mut r = Remote::new();
            for i in 0..100 {
                r.add(i, i as Rem, make_ext(i),
                      make_elt(i), i as u128);
            }
            for i in 0..100 {
                assert!(r.contains(i, i as Rem, make_ext(i),
                                   &make_elt(i), i as u128));
            }
        }
        #[test]
        fn test_get() {
            let mut r = Remote::new();
            for i in 0..10 {
                for j in 0..10 {
                    r.add(i, i as Rem, make_ext(i),
                          make_elt(j), j as u128);
                }
            }
            for i in 0..10 {
                assert_eq!(
                  r.get(i, i as Rem, make_ext(i)),
                  (0..10)
                      .map(|j| (make_elt(j), j as u128))
                      .collect::<Vec<(String, u128)>>()
                );
            }
        }
        #[test]
        fn test_update_ext() {
            let mut r = Remote::new();
            r.add(0,0,Ext::None, "a".to_string(), 1);
            r.update_ext(0, 0, Ext::None, Ext::Some{bits: 0, len: 1});
            assert_eq!(r.get(0, 0, Ext::None), vec![]);
            assert_eq!(r.get(0, 0, Ext::Some{bits: 0, len: 1}), vec![("a".to_string(), 1)])
        }
    }
}
use remote::Remote;

#[derive(Debug)]
pub struct AQF {
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
    remote: Remote<String>,
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
            remote: Remote::new(),
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
    fn inc_nelts(&mut self) {
        self.nelts += 1;
    }
}

impl AQF {
    fn calc_ext(&self, hash: u128, k: usize) -> u64 {
        debug_assert_ne!(k, 0);
        debug_assert!(self.q + self.r + k <= 64,
                      "Extension overflowed 64 bits at {}",
                      self.q + self.r + k <= 64);
        let hash = hash >> (self.q + self.r);
        (hash & b128::half_open(0, k)) as u64
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
    /// Adapt on a false match for all remote rep entries in pairs
    fn adapt(&mut self, loc: usize, query_quot: usize, query_rem: Rem, query_hash: u128,
             pairs: Vec<(String, u128)>, block_exts: [Ext; 64]) {
        let mut i = loc;                 // position in filter
        let mut pairs = pairs.iter();    // position in pairs
        let mut hash;
        match pairs.next() {
            Some((_, h)) => {
                hash = *h;
            }
            None => panic!("Empty pairs!")
        }
        let mut block_exts = block_exts; // in case the run spans multiple blocks
        loop {
            // Check if rem, ext match (quot already matches b/c we're in its run)
            if self.remainder(i) == query_rem && self.ext_matches(block_exts[i%64], query_hash) {
                // If we get a match, adapt on the current (elt, hash) pair
                // and advance to the next elt, hash pair
                self.adapt_loc(i, query_quot, query_rem, hash, query_hash, &mut block_exts);
                if let Some((_, remote_hash)) = pairs.next() {
                    hash = *remote_hash;
                } else {
                    // Exit when we've run out of (elt, hash) pairs
                    break;
                }
            }
            // Make sure we don't go off the end of the run:
            // we should terminate before this point by running out of pairs
            debug_assert!(i > query_quot && !self.is_runend(i-1),
                          "Went off the end of the run before running out of remote pairs!");
            i -= 1;
            // Update decoded block_exts if we enter a new block
            if i % 64 == 63 {
                block_exts = ExtensionArcd::decode(self.blocks[i/64].extensions);
            }
        }
    }
    /// Adapt on a false match for a fingerprint at loc
    fn adapt_loc(&mut self, loc: usize, quot: usize, rem: Rem, member_hash: u128, non_member_hash: u128, exts: &mut [Ext; 64]) {
        let new_ext = self.shortest_diff_ext(member_hash, non_member_hash);
        assert_ne!(new_ext, Ext::None,
                   "Hashes were identical, member_hash={}, non_member_hash={}",
                   member_hash, non_member_hash);
        // Write encoding to the appropriate block
        let old_ext = exts[loc%64];
        exts[loc%64] = new_ext;
        match ExtensionArcd::encode(*exts) {
            Ok(code) => {
                // Update code and add new extension to remote
                self.blocks[loc/64].extensions = code;
                self.remote.update_ext(quot, rem, old_ext, new_ext);
            }
            Err(_) => {
                // Encoding failed: rebuild
                // Clear all extensions in the offending block + in remote rep
                self.clear_block_remote_exts(loc/64);
                // Write new extension encoding
                *exts = [Ext::None; 64];
                exts[loc%64] = new_ext;
                match ExtensionArcd::encode(*exts) {
                    Ok(code) => {
                        self.blocks[loc/64].extensions = code;
                        // Add new extension back into the remote
                        self.remote.update_ext(quot, rem, Ext::None, new_ext);
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
    fn clear_block_remote_exts(&mut self, block_i: usize) {
        let exts = ExtensionArcd::decode(self.blocks[block_i].extensions);
        // Clear extensions from remote representation
        self.apply_to_block(
            block_i,
            |filter: &mut Self, quot: usize, i: usize| {
                let rem = filter.remainder(i);
                let ext = exts[i%64];
                // Clear extension if it is not already Ext::None
                if let Ext::Some{bits: _, len: _} = ext {
                    filter.remote.clear_ext(quot, rem, ext);
                }
            });
    }
    fn raw_query(&mut self, hash: u128, elt: &str) -> bool {
        let query_hash = hash;
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
                        // // Check cached code:
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
                        if self.ext_matches(ext, query_hash) {
                            // Extensions match => check remote to see if true match
                            let pairs = self.remote.get(quot, rem, ext);
                            debug_assert!(
                                !pairs.is_empty(),
                                "Each stored fingerprint must have a corresponding entry in the \
                                 remote rep: [quot={}, rem={}, ext={}] has no entry",
                                quot, rem, ext
                            );
                            // Check if any of the (elt, hash) pairs stored in the remote rep
                            // match the query elt
                            match pairs.iter().find(|(e, _)| *e == elt) {
                                None => {
                                    // If false match, adapt on all elt/hash pairs
                                    // FIXME: I had to change pairs to use strings instead of slices b/c of borrowing checks,
                                    //  and that might slow things down; might make sense to put everything in this block
                                    //  in the same function
                                    self.adapt(loc, quot, rem, query_hash, pairs, exts);
                                }
                                Some(_) => {}
                            }
                            return true;
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
                panic!("Failed to find runend for occupied quotient ({}), filter={:#?}",
                       quot, self);
            }
        }
    }
    /// Returns true if the extension is consistent with the query hash, false otherwise.
    ///
    /// For an extension of length `len` with bits `bits`, checks if the bits in the query hash
    /// from `q + r` to `q + r + len - 1`
    fn ext_matches(&self, ext: Ext, query_hash: u128) -> bool {
        match ext {
            Ext::Some {bits, len} => {
                // Compute extension from query hash and check if it matches the stored ext
                self.calc_ext(query_hash, len) == bits
            }
            Ext::None => {
                // No extension: extension matches by default
                true
            }
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
        self.raw_query(query_hash, &elt[..])
    }
    // FIXME: insert should shift extensions
    fn insert(&mut self, elt: String) {
        let hash = self.hash(&elt[..]);
        let quot = self.calc_quot(hash);
        let rem = self.calc_rem(hash);
        // Don't insert duplicates FIXME
        if !self.remote.contains(quot, rem, Ext::None, &elt[..], hash) {
            self.raw_insert(quot, rem);
            self.remote.add(quot, rem, Ext::None, elt, hash);
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
            filter.remote.add(quot, rem, Ext::Some{bits: 1, len: 1}, elt, hash);
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
        }
        // Clear and check remote
        for i in 0..filter.nblocks {
            filter.clear_block_remote_exts(i);
        }
        // Check that all extensions are cleared
        for (_, vec) in filter.remote.items() {
            for (ext, _, _) in vec {
                assert_eq!(*ext, Ext::None);
            }
        }
        // eprintln!("after clearing remote: filter: {:#?}", filter);
    }
    /// Use filter's hash fn but set quotient bits to 0, remainder bits to 1s
    fn hash_trunc(filter: &AQF, elt: &str) -> u128 {
        let (r, q) = (filter.r, filter.q);
        let hash = filter.hash(elt);
        (hash & !b128::half_open(0, q)) | b128::half_open(q, q+r)
    }
    fn build_and_test_repeating<F>(n: usize, ext_policy: F)
        where F: FnOnce(&AQF, u128, usize) -> Ext + Copy {
        // Build filter with n extensions, following `ext_policy`
        let mut filter = AQF::new(64, 4);
        let mut exts = [Ext::None; 64];
        for i in 0..n {
            let elt = i.to_string();
            let hash = hash_trunc(&filter, &elt[..]);
            let quot = filter.calc_quot(hash);
            let rem = filter.calc_rem(hash);
            let ext = ext_policy(&filter, hash, i);
            // eprintln!("i={}, ext={:?}", i, ext);
            exts[i%64] = ext;
            let code = ExtensionArcd::encode(exts).unwrap();
            filter.raw_insert(quot, rem);
            filter.remote.add(quot, rem, ext, elt, hash);
            filter.blocks[i/64].extensions = code;
        }
        // eprintln!("filter after inserts: {:#?}", filter);

        // Check that the filter contains elts without extensions
        for i in 0..n {
            let elt = &i.to_string()[..];
            let mut hash = hash_trunc(&filter, elt);
            if let Ext::Some{bits, len} = ext_policy(&filter, hash, i) {
                hash |= (bits as u128) << (filter.q + filter.r + len);
            }
            // eprintln!("i={}, hash=0b{:b}, filter={:#?}", i, hash, filter);
            assert!(
                filter.raw_query(hash, elt),
                "Filter doesn't contain elt {}",
                elt);
        }
    }
    /// Get the `k`-th to `(k+n-1)`-th bits in `bits`.
    fn get_bits(bits: u128, k: usize, n: usize) -> u64 {
        ((b128::half_open(k, k+n) & bits) >> k) as u64
    }
    #[test]
    fn test_repeating_sawtooth() {
        // Try inserting multiple elts that have the same quot, rem
        // where selectors are on/off by parity
        build_and_test_repeating(10, |filter: &AQF, hash: u128, i: usize| {
            let k = filter.q + filter.r;
            if i % 2 == 0 {
                Ext::Some{
                    bits: get_bits(hash, k, 1),
                    len: 1
                }
            } else {
                Ext::None
            }
        })
    }
    #[test]
    fn test_repeating_ascending() {
        // Try inserting multiple elts with same quot, rem
        // where selectors grow with index
        build_and_test_repeating(6, |filter: &AQF, hash: u128, i: usize| {
            let k = filter.q + filter.r;
            if i < 2 {
                Ext::None
            } else if i < 4 {
                Ext::Some {
                    bits: get_bits(hash, k, 1),
                    len: 1
                }
            } else {
                Ext::Some {
                    bits: get_bits(hash, k, 2),
                    len: 2
                }
            }
        });
    }
    #[test]
    fn test_insert_and_query() {
        // Insert and query elts, ensure that there are no false negatives
        let a: usize = 1 << 20;
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
