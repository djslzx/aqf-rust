#![feature(asm)] // Use assembly implementations of rank and select
#![allow(dead_code)]

/// Arithmetic code
mod arcd;

/// Bit manipulation utilities
mod util;

// TODO: change based on rem size
type Rem = u8;

trait Filter<T> {
    fn query(&mut self, elt: T) -> bool;
    fn insert(&mut self, elt: T);
}

/// Rank-and-Select Quotient Filter
mod rsqf;

/// Adaptive Quotient Filter
mod aqf;

/// Experiments
mod test;

fn main() {
    // test::file_fpr_test("data/chicagoA64.txt", 605006, 100_f64, 0.95, 8);
    // test::file_fpr_test("data/chicagoA64.txt", 605006, 1000_f64, 0.95, 8);
    // test::file_fpr_test("data/chicagoA64.txt", 605006, 10000_f64, 0.95, 8);

    // test::file_fpr_test("data/chicagoB64.txt",1700741 , 100_f64, 0.95, 8);
    test::file_fpr_test("data/chicagoB64.txt", 1700741, 1000_f64, 0.95, 8);
    test::file_fpr_test("data/chicagoB64.txt", 1700741, 10000_f64, 0.95, 8);

    // test::file_fpr_test("data/sanjose64.txt", 2193052, 100_f64, 0.95, 8);
    // test::file_fpr_test("data/sanjose64.txt", 2193052, 1000_f64, 0.95, 8);
    // test::file_fpr_test("data/sanjose64.txt", 2193052, 10000_f64, 0.95, 8);
}
