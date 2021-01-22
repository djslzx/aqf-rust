#![feature(asm)] // Use assembly implementations of rank and select
#![allow(dead_code)]

/// Arithmetic code
mod arcd;

/// Bit manipulation utilities
mod util;

// TODO: change based on rem size
type Rem = u8;

trait Filter<T> {
    fn query(&self, elt: T) -> bool;
    fn insert(&mut self, elt: T);
}

/// Rank-and-Select Quotient Filter
mod rsqf;

/// Adaptive Quotient Filter
mod aqf;

fn main() {

}
