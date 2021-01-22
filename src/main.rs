#![feature(asm)] // Use assembly implementations of rank and select
#![allow(dead_code)]

use murmur3;
use rand::{
    Rng, SeedableRng,
};
use rand::rngs::SmallRng;
use std::{cmp, collections::{HashSet, HashMap}};
mod arcd;
mod util;
use util::{
    bitarr::{b128, b64},
    bitrank, bitselect, popcnt,
};

// TODO: change based on rem size
type Rem = u8;

trait Filter<T> {
    fn query(&self, elt: T) -> bool;
    fn insert(&mut self, elt: T);
}

/// Adaptive Quotient Filter
mod aqf;

/// Rank-and-Select Quotient Filter
mod rsqf;

fn main() {

}
