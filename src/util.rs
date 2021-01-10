// Utility functions

use std::num::Wrapping;

// Bit arrays
pub mod bitarr {
    // 64-bit bit arrays
    pub mod b64 {
        pub fn get(x: u64, at: usize) -> bool {
            (x & (1 << at)) != 0
        }
        fn set(x: u64, at: usize) -> u64 {
            x | (1 << at)
        }
        fn unset(x: u64, at: usize) -> u64 {
            x & !(1 << at)
        }
        pub fn set_to(x: u64, on: bool, at: usize) -> u64 {
            if on {
                self::set(x, at)
            } else {
                self::unset(x, at)
            }
        }
        #[cfg(test)]
        mod tests {
            use super::*;

            #[test]
            fn test_b64_get() {
                let zeros = 0_u64;
                let ones = !0_u64;
                for i in 0..64 {
                    assert_eq!(get(zeros, i), false, "i={}", i);
                    assert_eq!(get(ones, i), true, "i={}", i);
                }
                assert_eq!(get(0b1101_1001_1000, 0), false);
                assert_eq!(get(0b1101_1001_1000, 1), false);
                assert_eq!(get(0b1101_1001_1000, 2), false);
                assert_eq!(get(0b1101_1001_1000, 3), true);
                assert_eq!(get(0b1101_1001_1000, 4), true);
                assert_eq!(get(0b1101_1001_1000, 5), false);
                assert_eq!(get(0b1101_1001_1000, 6), false);
                assert_eq!(get(0b1101_1001_1000, 7), true);
                assert_eq!(get(0b1101_1001_1000, 8), true);
                assert_eq!(get(0b1101_1001_1000, 9), false);
                assert_eq!(get(0b1101_1001_1000, 10), true);
                assert_eq!(get(0b1101_1001_1000, 11), true);
            }
            #[test]
            fn test_b64_set() {
                // 0 -> 1
                assert_eq!(set(0, 0), 1);
                assert_eq!(set(0, 1), 2);
                assert_eq!(set(0, 2), 4);
                assert_eq!(set(0, 4), 16);
                assert_eq!(set(0, 63), 1 << 63);
                // 0 -> 0
                assert_eq!(set(!0, 0), !0);
                assert_eq!(set(!0, 4), !0);
                assert_eq!(set(!0, 63), !0);
            }
            #[test]
            fn test_unset() {
                let ones = !0_u64;
                let zeros = 0_u64;
                for i in 0..64 {
                    assert_eq!(unset(ones, i), !(1 << i));
                    assert_eq!(unset(zeros, i), zeros);
                }
            }
            #[test]
            fn test_set_to() {
                let zeros = 0_u64;
                let ones = !0_u64;
                for i in 0..64 {
                    assert_eq!(set_to(zeros, true, i), 1 << i, "i={}", i);
                    assert_eq!(set_to(zeros, false, i), 0, "i={}", i);
                    assert_eq!(set_to(ones, true, i), !0, "i={}", i);
                    assert_eq!(set_to(ones, false, i), !(1 << i), "i={}", i);
                }
            }
        }
    }
    // 128-bit bit arrays
    pub mod b128 {
        /// Bits in [a,b) set to 1, others set to 0;
        /// Indexing from right (LSB has index 0)
        pub fn half_open(a: usize, b: usize) -> u128 {
            assert!(a != b);
            if b - a == 128 {
                !0_u128
            } else {
                ((1_u128 << (b - a)) - 1) << a
            }
        }
        /// Bits in [a,b] set to 1, others set to 0
        /// Indexing from right (LSB has index 0)
        pub fn closed(a: usize, b: usize) -> u128 {
            self::half_open(a, b + 1)
        }
        #[cfg(test)]
        mod tests {
            use super::*;

            #[test]
            fn test_half_open() {
                assert_eq!(half_open(0, 1), 1, "{:x}", half_open(0, 1));
                assert_eq!(half_open(127, 128), 1 << 127, "{:x}", half_open(127, 128));
                assert_eq!(half_open(0, 128), !0, "{:x}", half_open(0, 128));
                assert_eq!(half_open(1, 128), !1, "{:x}", half_open(1, 128));
                assert_eq!(half_open(0, 127), !(1 << 127), "{:x}", half_open(0, 127));
                assert_eq!(half_open(1, 127), !0 << 2 >> 1, "{:x}", half_open(1, 128));
            }
            #[test]
            fn test_closed() {
                assert_eq!(closed(0, 0), 1, "{:x}", closed(0, 0));
                assert_eq!(closed(127, 127), 1 << 127, "{:x}", closed(127, 127));
                assert_eq!(closed(0, 127), !0, "{:x}", closed(0, 127));
                assert_eq!(closed(1, 127), !1, "{:x}", closed(1, 127));
                assert_eq!(closed(0, 126), !(1 << 127), "{:x}", closed(0, 126));
                assert_eq!(closed(1, 126), !0 << 2 >> 1, "{:x}", closed(1, 126));
            }
        }
    }
}

// Round v to nearest power of 2
// Pre: v >= 0
// https://graphics.stanford.edu/~seander/bithacks.html
pub fn nearest_pow_of_2(mut v: usize) -> usize {
    v = v - 1;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v + 1
}

// Fast assembly rank and select
pub fn popcnt(val: u64) -> u64 {
    unsafe {
        let o: u64;
        asm!("popcnt {0}, {1}",
             out(reg) o,
             in(reg) val,
             options(preserves_flags),
        );
        o
    }
}

/// Counts the number of 1 bits in val.
/// Note: pos indexes from 0.
pub fn bitrank(val: u64, pos: usize) -> u64 {
    if pos >= 64 {
        popcnt(val)
    } else {
        let val = val & ((2 << pos) - 1);
        unsafe {
            let o: u64;
            asm!("popcnt {0}, {1}",
                 out(reg) o,
                 in(reg) val,
                 options(preserves_flags),
            );
            o
        }
    }
}

// Returns the position of the k-th 1 in the 64-bit word x.
// k is 0-based, so k=0 returns the position of the first 1.

// Uses the broadword selection algorithm by Vigna [1], improved by Gog
// and Petri [2] and Vigna [3].

// [1] Sebastiano Vigna. Broadword Implementation of Rank/Select
//    Queries. WEA, 2008

// [2] Simon Gog, Matthias Petri. Optimized succinct data
// structures for massive data. Softw. Pract. Exper., 2014

// [3] Sebastiano Vigna. MG4J 5.2.1. http://mg4j.di.unimi.it/
// The following code is adapted from
// https://github.com/facebook/folly/blob/b28186247104f8b90cfbe094d289c91f9e413317/folly/experimental/Select64.h

const K_SELECT_IN_BYTE: [u8; 2048] = [
    8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    8, 8, 8, 1, 8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    8, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    8, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 3, 8, 3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2,
    8, 8, 8, 5, 8, 5, 5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2,
    8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    8, 8, 8, 7, 8, 7, 7, 2, 8, 7, 7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2,
    8, 7, 7, 5, 7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2,
    7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3,
    8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8, 5, 5, 4, 5, 4, 4, 3,
    8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8, 8, 6, 8, 6, 6, 4, 8, 6, 6, 4, 6, 4, 4, 3,
    8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 3, 8, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7, 4, 4, 3,
    8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3,
    8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4, 6, 4, 4, 3,
    8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 4,
    8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 4,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 4,
    8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4,
    8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4,
    8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6,
    8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
];

fn _select64(x: u64, k: u64) -> u64 {
    if k >= popcnt(x) {
        64
    } else {
        let k = Wrapping(k);
        let k_ones_step4 = Wrapping(0x1111111111111111_u64);
        let k_ones_step8 = Wrapping(0x0101010101010101_u64);
        let k_msbs_step8 = Wrapping(0x80_u64) * k_ones_step8;

        let mut s = Wrapping(x);
        s = s - ((s & Wrapping(0xA) * k_ones_step4) >> 1);
        s = (s & Wrapping(0x3) * k_ones_step4) + ((s >> 2) & Wrapping(0x3) * k_ones_step4);
        s = (s + (s >> 4)) & Wrapping(0xF) * k_ones_step8;

        let byte_sums = s * k_ones_step8;
        let k_step8 = k * k_ones_step8;
        let geq_k_step8 = ((k_step8 | k_msbs_step8) - byte_sums) & k_msbs_step8;

        let place = popcnt(geq_k_step8.0) * 8;
        let byte_rank = k.0 - (((byte_sums.0 << 8) >> place) & 0xFF);
        place + (K_SELECT_IN_BYTE[(((x >> place) & 0xFF) | (byte_rank << 8)) as usize] as u64)
    }
}

/// Returns the position of the rank-th 1 bit in val.
/// Indexes from 0?
/// Equivalent to tzcnt(pdep(2^rank, val)) when SSE4.2 supported.
// Adapted from C:
// uint64_t i = 1ULL << rank;
// asm("pdep %[val], %[mask], %[val]"
//     : [val] "+r" (val)          // +r: register in a general register, read/write
//     : [mask] "r" (i));          // r: register in a general register
// asm("tzcnt %[bit], %[index]"
//     : [index] "=r" (i)          // =r: register in a general register, written to
//     : [bit] "g" (val)           // g: any registers except those that aren't general registers
//     : "cc");
pub fn bitselect(val: u64, rank: u64) -> u64 {
    if is_x86_feature_detected!("sse4.2") {
        let mut i = 1_u64 << rank; // 2^rank
        unsafe {
            asm!("pdep {0}, {1}, {0}",
                 "tzcnt {1}, {0}",
                 in(reg) val,
                 inout(reg) i,
                 options(preserves_flags),
            );
            i
        }
    } else {
        _select64(val, rank)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nearest_pow_2() {
        assert_eq!(nearest_pow_of_2(1), 1);
        assert_eq!(nearest_pow_of_2(2), 2);
        assert_eq!(nearest_pow_of_2(3), 4);
        assert_eq!(nearest_pow_of_2(4), 4);
        assert_eq!(nearest_pow_of_2(5), 8);
        assert_eq!(nearest_pow_of_2(8), 8);
        assert_eq!(nearest_pow_of_2(9), 16);
        assert_eq!(nearest_pow_of_2(16), 16);
        assert_eq!(nearest_pow_of_2(17), 32);
        assert_eq!(nearest_pow_of_2(32), 32);
        assert_eq!(nearest_pow_of_2(33), 64);
        assert_eq!(nearest_pow_of_2(65), 128);
        assert_eq!(nearest_pow_of_2(129), 256);
    }
    #[test]
    fn test_popcnt() {
        assert_eq!(popcnt(0), 0);
        assert_eq!(popcnt(1), 1);
        assert_eq!(popcnt(2), 1);
        assert_eq!(popcnt(3), 2);
        assert_eq!(popcnt(0x10101010101010), 7);
        assert_eq!(popcnt(0x10001000100010), 4);
    }
    #[test]
    fn test_bitrank() {
        assert_eq!(bitrank(0, 0), 0);
        assert_eq!(bitrank(0, 64), 0);
        assert_eq!(bitrank(1, 0), 1);
        assert_eq!(bitrank(1, 64), 1);
        assert_eq!(bitrank(0x100, 0), 0);
        assert_eq!(bitrank(0x100, 1), 0);
        assert_eq!(bitrank(0x100, 16), 1);
        assert_eq!(bitrank(0x100, 17), 1);
        assert_eq!(bitrank(0x100, 64), 1);
    }
    #[test]
    fn test_select64() {
        for i in 0..64 {
            // No 1s
            assert_eq!(_select64(0, i), 64);
            // One 1
            for j in 0..64 {
                assert_eq!(_select64(1 << j, i), if i < 1 { j } else { 64 });
            }
            // Two 1s
            for a in 0..64 {
                for b in (a + 1)..64 {
                    let x = (1 << a) | (1 << b);
                    assert_eq!(
                        _select64(x, i),
                        if i < 1 {
                            a
                        } else if i < 2 {
                            b
                        } else {
                            64
                        },
                        "x=0x{:x}, a={}, b={}",
                        x,
                        a,
                        b
                    );
                }
            }
            // Hard-coded
            assert_eq!(_select64(0b10_100_010_101_001, 0), 0);
            assert_eq!(_select64(0b10_100_010_101_001, 1), 3);
            assert_eq!(_select64(0b10_100_010_101_001, 2), 5);
            assert_eq!(_select64(0b10_100_010_101_001, 3), 7);
            assert_eq!(_select64(0b10_100_010_101_001, 4), 11);
            assert_eq!(_select64(0b10_100_010_101_001, 5), 13);
        }
    }
    #[test]
    fn test_bitselect() {
        assert_eq!(bitselect(0, 1), 64);
        for i in 0..64 {
            // No 1s
            assert_eq!(bitselect(0, i), 64, "i={}", i);
            // One 1
            for j in 0..64 {
                assert_eq!(bitselect(1 << j, i), if i < 1 { j } else { 64 });
            }
            // Two 1s
            for a in 0..64 {
                for b in (a + 1)..64 {
                    let x = (1 << a) | (1 << b);
                    assert_eq!(
                        bitselect(x, i),
                        if i < 1 {
                            a
                        } else if i < 2 {
                            b
                        } else {
                            64
                        },
                        "x=0x{:x}, a={}, b={}",
                        x,
                        a,
                        b
                    );
                }
            }
            // Hard-coded
            assert_eq!(bitselect(0b10_100_010_101_001, 0), 0);
            assert_eq!(bitselect(0b10_100_010_101_001, 1), 3);
            assert_eq!(bitselect(0b10_100_010_101_001, 2), 5);
            assert_eq!(bitselect(0b10_100_010_101_001, 3), 7);
            assert_eq!(bitselect(0b10_100_010_101_001, 4), 11);
            assert_eq!(bitselect(0b10_100_010_101_001, 5), 13);
        }
    }
}
