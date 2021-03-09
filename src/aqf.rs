use std::{
    cmp::{max, min},
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
    use std::collections::HashSet;

    type Quot = usize;

    /// Remote representation mapping from position to (elt, hash)
    pub struct IndexRemote<T> {
        data: HashMap<usize, (T, u128)>, // Store (elt, hash) pairs
        inserted: HashSet<T>, // Track which elements have been inserted
    }
    impl IndexRemote<String> {
        pub fn new(size: usize) -> Self {
            let mut data: HashMap<usize, (String, u128)> = HashMap::with_capacity(size);
            for i in 0..size {
                data.insert(i, ("".to_string(), 0));
            }
            IndexRemote {
                data,
                inserted: HashSet::with_capacity(size),
            }
        }
        pub fn insert(&mut self, i: usize, s: String, hash: u128) {
            self.data.insert(i, (s.clone(), hash));
            self.inserted.insert(s);
        }
        pub fn get(&mut self, i: usize) -> Option<&(String, u128)> {
            self.data.get(&i).clone()
        }
        pub fn remove(&mut self, i: usize) -> (String, u128) {
            debug_assert!(self.data.contains_key(&i));
            let (elt, hash) = self.data.remove(&i).unwrap();
            self.inserted.remove(&elt);
            (elt, hash)
        }
        pub fn contains(&self, s: &str) -> bool {
            self.inserted.contains(s)
        }
        pub fn is_empty(&self) -> bool {
            self.data.is_empty()
        }
    }
    impl fmt::Debug for IndexRemote<String> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.debug_list()
                .entries(self.data.iter().filter_map(|(i, (elt, hash))|
                    if elt != "" || *hash != 0 {
                        Some(format!("i={}, elt={}, hash=0x{:x}", i, elt, hash))
                    } else {
                        None
                    }
                ))
                .finish()
        }
    }

    pub struct Remote<T> {
        pub data: HashMap<(Quot, Rem), Vec<(Ext, T, u128)>>,
    }
    impl fmt::Debug for Remote<String> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.debug_map()
                .entries(self.data.iter().map(|(&(q, r), vec)| {
                    (
                        format!("{}, {}", q, r),
                        vec.iter()
                            .map(|(e, s, h)| format!("{}, {}, 0x{:x}", e, s, h))
                            .collect::<Vec<String>>()
                    )
                }))
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
        /// Returns true if the remote rep contains no elements
        pub fn is_empty(&self) -> bool {
            self.data.is_empty()
        }
        /// Get a vector of all the items in the remote rep
        pub fn items(&self) -> Vec<(&(Quot, Rem), &Vec<(Ext, String, u128)>)> {
            self.data.iter().collect()
        }
        /// Inserts the mapping `(quot, rem, ext) -> (elt, hash)` into the remote representation
        pub fn add(&mut self, quot: Quot, rem: Rem, ext: Ext, elt: String, hash: u128) {
            // TODO: check that quot, rem, ext are consistent with hash (requires hooking up with global state)
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
        /// Returns true if `e1` and `e2` are consistent.
        fn consistent(e1: Ext, e2: Ext) -> bool {
            if let Ext::Some{bits: bits1, len: len1} = e1 {
                if let Ext::Some{bits: bits2, len: len2} = e2 {
                    let len = min(len1, len2);
                    return (bits1 ^ bits2) & b64::mask(len) == 0;
                }
            }
            true
        }
        /// Returns the `(elt, hash)` pairs associated with a `(quot, rem, ext)` triple
        /// by choosing elements associated with `(quot, rem)` that are consistent with `ext`
        pub fn get_consistent(&self, quot: Quot, rem: Rem, ext: Ext) -> Vec<(String, u128)> {
            match self.data.get(&(quot, rem)) {
                None => Vec::new(),
                Some(vec) => {
                    vec.iter()
                        .filter_map(|(xt, elt, hash)|
                            if Remote::consistent(*xt, ext) {
                                Some((elt.clone(), *hash))
                            } else {
                                None
                            })
                        .collect()
                }
            }
        }
        /// Returns all `(elt, hash)` pairs that match the fingerprint `quot:rem:ext`.`
        pub fn get_matching(&self, quot: Quot, rem: Rem, ext: Ext) -> Vec<(String, u128)> {
            match self.data.get(&(quot, rem)) {
                None => Vec::new(),
                Some(vec) => {
                    vec.iter()
                        .filter_map(|(xt, elt, hash)|
                            if *xt == ext {
                                Some((elt.clone(), *hash))
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

        fn zero_ext(len: usize) -> Ext {
            if len == 0 {
                Ext::None
            } else {
                Ext::Some {
                    bits: 0,
                    len,
                }
            }
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
                r.add(i, i as Rem, zero_ext(i),
                      "a".repeat(i), i as u128);
            }
            for i in 0..100 {
                assert!(r.contains(i, i as Rem, zero_ext(i),
                                   &"a".repeat(i), i as u128));
            }
        }
        #[test]
        fn test_get_matching() {
            let mut r = Remote::new();
            for i in 0..10 {
                for j in 0..10 {
                    r.add(i, i as Rem, zero_ext(i), "a".repeat(j), j as u128);
                }
            }
            for i in 0..10 {
                assert_eq!(
                    r.get_consistent(i, i as Rem, zero_ext(i)),
                    (0..10)
                      .map(|j| ("a".repeat(j), j as u128))
                      .collect::<Vec<(String, u128)>>()
                );
            }
        }
        #[test]
        fn test_consistent() {
            assert!(Remote::consistent(Ext::None,
                                       Ext::None));
            assert!(Remote::consistent(Ext::None,
                                       Ext::Some{bits: 0, len: 1}));
            assert!(Remote::consistent(Ext::Some{bits: 0, len: 1},
                                       Ext::None));
            assert!(Remote::consistent(Ext::Some{bits: 0, len: 1},
                                       Ext::Some{bits: 0, len: 10}));
            assert!(Remote::consistent(Ext::Some{bits: 0b000010101, len: 5},
                                       Ext::Some{bits: 0b111110101, len: 10}));
            assert!(!Remote::consistent(Ext::Some{bits: 0b000010101, len: 6},
                                        Ext::Some{bits: 0b111110101, len: 10}));
            assert!(!Remote::consistent(Ext::Some{bits: 0, len: 1},
                                        Ext::Some{bits: 1, len: 1}));
        }
        #[test]
        fn test_get_consistent() {
            // Insert elements with consistent but not matching fingerprints
            // and make sure we retrieve the consistent ones
            let mut r = Remote::new();
            for i in 0..3 {
                r.add(0, 0, zero_ext(i), "a".repeat(i), i as u128);
            }
            r.add(0, 0, Ext::Some{bits: 1, len: 1}, "b".to_string(), 0);

            assert_eq!(
                r.get_consistent(0, 0, Ext::None),
                vec![
                    ("a".repeat(0), 0),
                    ("a".repeat(1), 1),
                    ("a".repeat(2), 2),
                    ("b".to_string(), 0)
                ]
            );
            for i in 1..10 {
                assert_eq!(
                    r.get_consistent(0, 0, zero_ext(i)),
                    (0..3)
                        .map(|i| ("a".repeat(i), i as u128))
                        .collect::<Vec<(String, u128)>>()
                );
            }
        }
        #[test]
        fn test_update_ext() {
            let mut r = Remote::new();
            r.add(0, 0, Ext::None, "a".to_string(), 1);
            r.add(0, 0, Ext::None, "b".to_string(), 2);

            // Update ext for one of a or b
            r.update_ext(0, 0, zero_ext(0), zero_ext(1));
            assert_eq!(*r.data.get(&(0, 0)).unwrap(),
                       vec![
                           (zero_ext(1), "a".to_string(), 1),
                           (zero_ext(0), "b".to_string(), 2)
                       ]);

            // Remove ext for whichever of a or b was extended
            r.update_ext(0, 0, Ext::Some{bits: 0, len: 1}, Ext::None);
            assert_eq!(*r.data.get(&(0, 0)).unwrap(),
                       vec![
                         (zero_ext(0), "a".to_string(), 1),
                         (zero_ext(0), "b".to_string(), 2)
                       ]);
        }
    }
}
use remote::IndexRemote;

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
    pub remote: IndexRemote<String>,
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
            remote: IndexRemote::new(nslots),
            //remote: Remote::new(),
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
    pub fn calc_ext(&self, hash: u128, length: usize) -> u64 {
        debug_assert_ne!(length, 0, "Can only calculate ext for non-zero lengths");
        debug_assert!(self.q + self.r + length <= 64,
                      "Extension overflowed 64 bits at {}",
                      self.q + self.r + length <= 64);
        let hash = hash >> (self.q + self.r);
        (hash & b128::mask(length)) as u64
    }
    /// Generate the shortest extension from the member's hash that doesn't conflict with
    /// the non-member's hash. 
    /// `prev_ext` is the extension previously associated with the member's hash.
    fn shortest_diff_ext(&self, member_hash: u128, non_member_hash: u128) -> Ext {
        // Shift to get rid of first q+r bits in both hashes
        let a = member_hash >> (self.q + self.r);
        let b = non_member_hash >> (self.q + self.r);
        if a == b {
            Ext::None
        } else {
            // Find fewest LSBs needed to distinguish member from non-member hash:
            // Determine number of common LSBs and add 1
            let len = ((a ^ b).trailing_zeros() + 1) as usize;
            let bits = a & b128::mask(len); // mask len bits from member hash
            Ext::Some {
                bits: bits as u64,
                len: len as usize,
            }
        }
    }
    /// Returns true if `elt` is in the filter, false otherwise.
    ///
    /// Determines whether `elt` is in the filter by checking whether there is
    /// a remote element associated with the triple `(quot, rem, ext)`.
    // fn is_true_match(&self, elt: &str, quot: usize, rem: Rem, ext: Ext) -> bool {
    //     unimplemented!();
    //     // let matched_elts = self.remote.get_matching(quot, rem, ext);
    //     // debug_assert!(!matched_elts.is_empty());
    //     // let exact_match = matched_elts.iter().find(|(e, _)| *e == elt);
    //     // exact_match != None
    // }
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
    /// Adapt on a false match at `loc`.
    fn adapt(&mut self, query_elt: &str, loc: usize, quot: usize, rem: Rem, ext: Ext, hash: u128, exts: [Ext; 64]) {
        // eprintln!("Calling adapt with loc={}, quot={}, rem={:x}, ext={}, hash={:x}, exts={:?}",
        //           loc, quot, rem, ext, hash, exts);
        // eprintln!("before adapt, block#{}={:#?}", loc/64, self.block(loc/64));
        let mut exts = exts;
        // Collect locations that match rem, ext
        // TODO: put into helper fn
        let locs: Vec<usize> = (quot..=loc)
            .rev()// iterate backwards over [loc, quot]
            .take_while(|&i| i == loc || !self.is_runend(i))// don't step back further if we see a runend
            .filter(|&i| {
                // Re-decode at i if we're in a new block
                if i != loc && i % 64 == 63 {
                    exts = ExtensionArcd::decode(self.blocks[i/64].extensions);
                }
                let i_rem = self.remainder(i);
                let i_ext = exts[i%64];
                i_rem == rem && self.ext_matches(i_ext, hash)
            })
            .collect();
        // eprintln!("is={:?}", (quot..=loc).rev().take_while(|&i| i == loc || !self.is_runend(i)).collect::<Vec<usize>>());
        // eprintln!("locs={:?}", locs);

        // Make sure that the query elt isn't mapped to an earlier index in the sequence
        for &loc in &locs {
            if let Some((elt, _)) = self.remote.get(loc) {
                if elt == query_elt {
                    return;
                }
            }
        }
        for loc in locs {
            let stored_hash = match self.remote.get(loc) {
                None => {
                    panic!("Failed to get hash associated with loc")
                }
                Some((_, h)) => *h,
            };
            // eprintln!("loc={}, hash=0x{:x}", loc, stored_hash);
            // eprintln!("adapt_loc with loc={}, stored_hash=0x{:x}, query_hash=0x{:x}", loc, stored_hash, hash);
            self.adapt_loc(loc, stored_hash, hash);
        }
        // eprintln!("after adapt, block#{}={:#?}", loc/64, self.block(loc/64));
    }
    /// Adapt on a false match for a fingerprint at loc
    fn adapt_loc(&mut self, loc: usize, in_hash: u128, out_hash: u128) {
        let new_ext = self.shortest_diff_ext(in_hash, out_hash);
        if new_ext == Ext::None {
            println!("Hashes were identical, member_hash={:x}, non_member_hash={:x}, block={:#?}",
                     in_hash, out_hash, self.block(loc/64),
            );
            println!("remote={:#?}",
                     (0..64)
                         .map(|i| {
                             let (elt, hash) = self.remote.get(loc/64 + i).unwrap();
                             (elt.to_string(), *hash)
                         })
                         .filter(|(e, _h)| *e != "")
                         .map(|(e, h)| format!("elt={}, h=0x{:x}", e, h))
                         .collect::<Vec<String>>()
            );
            panic!("hashes were identical!!!")
        }
        // Write encoding to the appropriate block
        let mut exts = ExtensionArcd::decode(self.blocks[loc/64].extensions);
        exts[loc%64] = new_ext;
        match ExtensionArcd::encode(exts) {
            Ok(code) => {
                // Update arithmetic code
                self.blocks[loc/64].extensions = code;
            }
            Err(_) => {
                // Encoding failed: rebuild
                // Clear all extensions in the offending block in remote rep
                //self.clear_block_remote_exts(loc/64);
                // Write new extension encoding where loc has the new extension
                // and all other extensions are cleared
                exts = [Ext::None; 64];
                exts[loc%64] = new_ext;
                match ExtensionArcd::encode(exts) {
                    Ok(code) => {
                        // Write new arithmetic code to block
                        self.blocks[loc/64].extensions = code;
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
                // Clear extension if it is not already Ext::None
                // if let Ext::Some{bits: _, len: _} = ext {
                //     filter.remote.clear_ext(quot, rem, ext);
                // }
            });
    }
    /// Shift the remote elements in [a,b] forward by 1 into [a+1, b+1]
    fn shift_remote_elts(&mut self, a: usize, b: usize) {
        for i in (a..=b).rev() {
            let (elt, hash) = self.remote.remove(i);
            self.remote.insert(i+1, elt, hash);
        }
    }
    fn raw_query(&mut self, elt: &str, hash: u128) -> bool {
        let quot = self.calc_quot(hash);
        let rem = self.calc_rem(hash);

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
                        let exts = decode.unwrap().1;
                        let ext = exts[loc%64];
                        if self.ext_matches(ext, hash) {
                            if let Some((e, _)) = self.remote.get(loc) {
                                if *e != elt {
                                    // eprintln!("Calling adapt for elt={} b/c it doesn't match e={}, hash=0x{:x}", elt, *e, hash);
                                    self.adapt(elt, loc, quot, rem, ext, hash, exts);
                                }
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
    /// Insert a (quot, rem) pair into filter
    // TODO: take in ext?
    // FIXME: insert should shift extensions
    fn raw_insert(&mut self, quot: usize, rem: Rem, elt: String, hash: u128) {
        // Check that quot, rem are consistent with hash
        debug_assert_eq!(self.calc_quot(hash), quot);
        debug_assert_eq!(self.calc_rem(hash), rem);

        // Check that elt hasn't already been inserted
        if self.remote.contains(&elt) {
            return;
        }

        // Insert
        assert!(quot < self.nslots());
        self.inc_nelts();

        // Find the appropriate runend
        match self.rank_select(quot) {
            RankSelectResult::Empty => {
                // Insert a new singleton run at its home slot
                // (Doesn't need to modify offsets)
                self.set_occupied(quot, true);
                self.set_runend(quot, true);
                self.set_remainder(quot, rem);
                self.remote.insert(quot, elt.clone(), hash);
                debug_assert!(self.remote.contains(&elt), "Remote doesn't contain elt!");
            }
            RankSelectResult::Full(r) => {
                // Find u, the first open slot after r, and
                // shift everything in [r+1, u-1] forward by 1 into [r+2, u],
                // leaving r+1 writable
                let u = match self.first_unused_slot(r + 1) {
                    Some(loc) => loc,
                    None => {
                        // Extend the filter by one block
                        // and return the first empty index
                        self.add_block();
                        self.nslots() - 64
                    }
                };
                self.inc_offsets(r + 1, u - 1);
                self.shift_remainders_and_runends(r + 1, u - 1);
                self.shift_remote_elts(r + 1, u - 1);

                // Start a new run or extend an existing one
                if self.is_occupied(quot) {
                    // quot occupied: extend an existing run
                    self.inc_offsets(r, r);
                    self.set_runend(r, false);
                    self.set_runend(r + 1, true);
                    self.set_remainder(r + 1, rem);
                } else {
                    // quot unoccupied: start a new run
                    self.inc_offsets_for_new_run(quot, r);
                    self.set_occupied(quot, true);
                    self.set_runend(r + 1, true);
                    self.set_remainder(r + 1, rem);
                }
                // Insert element into remote rep at r+1
                self.remote.insert(r+1, elt.clone(), hash);
                debug_assert!(self.remote.contains(&elt), "Remote doesn't contain elt!");
            }
            RankSelectResult::Overflow =>
                panic!(
                    "RSQF failed to find runend (nslots={}, quot=(block={}, slot={}))",
                    self.nslots(), quot / 64, quot % 64,
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
        let query_hash = self.hash(&elt);
        self.raw_query(&elt, query_hash)
    }
    fn insert(&mut self, elt: String) {
        let hash = self.hash(&elt[..]);
        let quot = self.calc_quot(hash);
        let rem = self.calc_rem(hash);
        AQF::raw_insert(self, quot, rem, elt, hash);
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
    #[ignore] // index remote
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
    #[ignore] // index remote
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
            RankSelectQuotientFilter::raw_insert(&mut filter, quot, rem);
            //filter.remote.add(quot, rem, Ext::Some{bits: 1, len: 1}, elt, hash);
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
        // for (_, vec) in filter.remote.items() {
        //     for (ext, _, _) in vec {
        //         assert_eq!(*ext, Ext::None);
        //     }
        // }
        // eprintln!("after clearing remote: filter: {:#?}", filter);
    }
    /// Use filter's hash fn but set quotient bits to 0, remainder bits to 1s
    fn hash_trunc(filter: &AQF, elt: &str) -> u128 {
        let (r, q) = (filter.r, filter.q);
        let hash = filter.hash(elt);
        (hash & !b128::half_open(0, q)) | b128::half_open(q, q+r)
    }
    #[ignore]
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
            RankSelectQuotientFilter::raw_insert(&mut filter, quot, rem);
            //filter.remote.add(quot, rem, ext, elt, hash);
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
                filter.raw_query(elt, hash),
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
    /// Make a hash for elt where quot is 0 and rem is all 1s
    fn fake_hash(filter: &AQF, elt: &str) -> u128 {
        let hash = filter.hash(elt);
        (hash | b128::mask(filter.r)) << filter.q
    }
    #[test]
    fn test_adapt_1() {
        // Add multiple elements with the same quot, rem
        let mut filter = AQF::new(64, 4);
        for i in 0..10 {
            let elt = format!("elt[{}]", i);
            let hash = fake_hash(&filter, &elt);
            let quot = filter.calc_quot(hash);
            assert_eq!(quot, 0);
            let rem = filter.calc_rem(hash);
            assert_eq!(rem, 0b1111);
            filter.raw_insert(quot, rem, elt, hash);
        }
        eprintln!("filter={:#?}", filter);
        // Query on elts w/ same quot, rem
        for i in 11..20 {
            let elt = format!("elt[{}]", i);
            let hash = fake_hash(&filter, &elt);
            eprintln!("querying with elt={}, hash={:x}", elt, hash);
            filter.raw_query(&elt, hash);
            eprintln!("after querying elt={}, filter={:#?}", elt, filter);
            for j in 0..10 {
                let elt = format!("elt[{}]", j);
                let hash = fake_hash(&filter, &elt);
                assert!(filter.raw_query(&elt, hash),
                        "filter should still contain elt={}, hash={:x}",
                        elt, hash);
            }
        }
    }
    #[test]
    fn test_adapt_2() {
        let q = 6;
        let r = 4;
        let mut filter = AQF::new(1 << 6, r);
        assert_eq!(q, filter.q);
        let elt = "elt[0]".to_string();
        filter.raw_insert(0, 0, elt, 0);
        eprintln!("filter={:#?}", filter);

        // Adapt by adding a 0 to ext
        filter.raw_query("elt[1]", !0 << (q + r));
        eprintln!("filter={:#?}", filter);
        let exts = ExtensionArcd::decode(filter.blocks[0].extensions);
        assert_eq!(exts[0], Ext::Some{bits: 0, len: 1});

        // Adapt by adding two 0s to ext
        filter.raw_query("elt[2]", !0 << (q + r + 1));
        eprintln!("filter={:#?}", filter);
        let exts = ExtensionArcd::decode(filter.blocks[0].extensions);
        assert_eq!(exts[0], Ext::Some{bits: 0, len: 2});
    }
    #[test]
    fn test_adapt_3() {
        let q = 7;
        let r = 4;
        let mut filter = AQF::new(128, r);
        assert_eq!(q, filter.q);
        for i in 0..5 {
            let elt = format!("elt[{}]", i);
            let hash = (filter.hash(&elt) << 11) | 0b111_1000_0000;
            let quot = filter.calc_quot(hash);
            assert_eq!(quot, 0);
            let rem = filter.calc_rem(hash);
            assert_eq!(rem, 0b1111);
            filter.raw_insert(quot, rem, elt, hash);
        }
        eprintln!("filter={:#?}", filter);
        {
            let elt = "elt[6]";
            let hash = (filter.hash(elt) << 11) | 0b111_1000_0000;
            filter.raw_query(elt, hash);
            eprintln!("filter={:#?}", filter);
        }


        // eprintln!("filter={:#?}", filter);
        // // Query on elts w/ same quot, rem
        // for i in 11..20 {
        //     let elt = format!("elt[{}]", i);
        //     let hash = fake_hash(&filter, &elt);
        //     eprintln!("querying with elt={}, hash={:x}", elt, hash);
        //     filter.raw_query(&elt, hash);
        //     eprintln!("after querying elt={}, filter={:#?}", elt, filter);
        //     for j in 0..10 {
        //         let elt = format!("elt[{}]", j);
        //         let hash = fake_hash(&filter, &elt);
        //         assert!(filter.raw_query(&elt, hash),
        //                 "filter should still contain elt={}, hash={:x}",
        //                 elt, hash);
        //     }
        // }
    }

}
