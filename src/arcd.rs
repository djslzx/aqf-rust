// Arithmetic coding

// Adaptivity bits are represented as an Option<(u64, usize)>
// where the u64 holds the bits and the usize is the number of bits

const CODE_LEN: usize = 56;     // length of arithmetic code
const LG_ADAPTS: usize = 2;     // lg(adapt_rate)
const LG_EPS: usize = 4;        // -lg(eps)

pub mod extensions {
    use std::fmt;

    #[derive(Clone, Copy, PartialEq, Eq, Hash)]
    pub enum Ext {
        None,
        Some {bits: u64, len: usize},
    }

    /// Character to use when printing Ext::None
    const NONE_STR: &str = "_";

    impl Ext {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match *self {
                Ext::None => f.write_str(NONE_STR),
                Ext::Some { bits, len } => f.write_str(&format!("{:01$b}", bits, len)),
            }
        }
    }
    impl fmt::Debug for Ext {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            self.fmt(f)
        }
    }
    impl fmt::Display for Ext {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            self.fmt(f)
        }
    }
}
use extensions::*;

/// Error type for encoding failure
#[derive(Debug, Clone)]
pub struct EncodingFailure;

/// Arithmetic code
pub trait Arcd<I, O> {
    fn encode(input: [I; 64]) -> Result<O, EncodingFailure>;
    fn decode(input: O) -> [I; 64];
}

/// Computes floor(log2(x)) by counting the number of leading 0s
/// to get the position of the first leading 1
fn floor_log(x: u64) -> u32 {
    if x == 0 { 0 } // treat log(0) as 0
    else { 63 - x.leading_zeros() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floor_log() {
        for i in 0..10_000 {
            assert_eq!(floor_log(i),
                       (i as f64).log2().floor() as u32,
                       "i={}",
                       i);
        }
    }
}

/// Arithmetic code used for extension-based adaptivity
pub mod ext_arcd {
    use super::*;

    pub struct ExtensionArcd;

    impl Arcd<Ext, u64> for ExtensionArcd {
        fn encode(input: [Ext; 64]) -> Result<u64, EncodingFailure> {
            let mut low: u64 = 0;
            let mut high: u64 = !0 >> (64 - CODE_LEN);

            for ext in input.iter() {
                let range = high - low;
                // Multiply range by ~0.90624 (Pr[ext is empty])
                let mut gap = (range >> 1) + (range >> 2) + (range >> 3) + (range >> 5);
                match ext {
                    Ext::None => {
                        // If extension is empty, lower top of range
                        high = low + gap;
                    }
                    Ext::Some {bits, len} => {
                        // If extension is nonempty, raise bottom of range
                        low += gap;
                        // Multiply gap by ~0.4687
                        gap = (range >> 5) + (range >> 6);
                        // Account for probability of extension length:
                        // extension length k>0 has probability 2^{-k};
                        for _ in 0..(*len-1) {
                            low += gap;
                            gap >>= 1;
                        }
                        // Account for probability of a particular extension of length k
                        // (all equally likely -> 1/2^k)
                        gap >>= len;
                        low += bits * gap; // take bits-th 1/2^k-long piece
                        high = low + gap;
                    }
                }
                if high - low < 2 {
                    return Err(EncodingFailure);
                }
            }
            Ok(low)
        }
        fn decode(input: u64) -> [Ext; 64] {
            let mut low: u64 = 0;
            let mut high: u64 = !0 >> (64 - CODE_LEN);
            let mut out = [Ext::None; 64];

            for i in 0..64 {
                let range = high - low;
                // Multiply range by ~0.90624 (Pr[ext is empty])
                let mut gap = (range >> 1) + (range >> 2) + (range >> 3) + (range >> 5);
                out[i] = if low + gap > input {
                    high = low + gap;
                    Ext::None
                } else {
                    low += gap;
                    // Multiply gap by ~0.4687
                    gap = (range >> 5) + (range >> 6);
                    // Compute k, the length of the extension, by
                    // iteratively shrinking the gap in proportion to
                    // the probability of each length k=1, 2, ...
                    let mut len = 1;
                    while low + gap <= input {
                        low += gap;
                        gap >>= 1;
                        len += 1;
                    }
                    // Get the bits given k, the length of the extension,
                    // by dividing the interval into 2^k even pieces and
                    // determining which piece the input belongs to
                    gap >>= len; // divide gap by 2^k
                    let bits = (input - low)/gap;
                    low += bits * gap; // low = low_0 + bits * gap
                    high = low + gap;  // high = low_0 + (bits+1) * gap
                    Ext::Some{bits, len}
                };
            }
            out
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// Convert a test case string into an array of extensions
        fn str_to_ext_arr(s: &str) -> [Ext; 64] {
            let mut out = [Ext::None; 64];
            for (i, str) in s.split('|').enumerate() {
                out[i] = match str.parse::<u64>() {
                    Err(_) => Ext::None,
                    Ok(digit) => Ext::Some {
                        bits: digit,
                        len: (floor_log(digit) + 1) as usize
                    }
                };
                // print!("[str={}, ext={:?}]", str, out[i]);
            }
            // println!();
            out
        }

        // Check for each x in inputs that x = decode(encode(x))
        fn test_encode_decode_w_inputs(inputs: Vec<[Ext; 64]>) {
            for input in inputs {
                match ExtensionArcd::encode(input) {
                    // Encoding succeeds
                    Ok(code) => {
                        let decoded = ExtensionArcd::decode(code);
                        assert_eq!(
                            input, decoded,
                            "x={:?}, encode(x)={}, decode(encode(x))={:?}",
                            input, code, decoded,
                        );
                        println!("x={:?}, encode(x)={}",
                                 input.iter()
                                     .filter(|&ext|
                                         !matches!(*ext, Ext::None)
                                     )
                                     .collect::<Vec<&Ext>>(),
                                 code);
                    },
                    // Encoding fails (out of bits)
                    Err(_) => {
                        // Failed to encode but sequence encoding shouldn't overflow
                        // => error
                        panic!("Failed to encode sequence (ran out of bits): {:?}",
                               input.iter()
                                   .filter(|&ext|
                                       !matches!(*ext, Ext::None)
                                   )
                                   .collect::<Vec<&Ext>>());
                    }
                }
            }
        }

        #[test]
        fn test_encode_decode() {
            // '-': empty extension, '|': extension separator
            let inputs = [
                "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-",
                "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|0",
                "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|2|-|-|-|-|-|-|-|-|3|-|-|-|-|1|-|-|-|-|-|0",
                "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|11|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|3|-|-|-|-|-|-|-|-|-|-|13|-|-|-|-|-|-|-|-|-|-|-",
                "-|-|-|-|-|-|-|-|-|-|-|-|5|-|-|-|-|-|-|-|-|5|-|-|-|-|-|-|-|5|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-",
                "-|-|-|-|26|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-",
                "1|1|1|1|1|1|1|1|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-", // passes (8 1s)
                "1|1|1|1|1|1|1|1|1|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-", // fails (9 1s)
            ];

            // Convert inputs into ext array
            let mut ext_inputs = Vec::with_capacity(inputs.len());
            for input in inputs.iter() {
                ext_inputs.push(str_to_ext_arr(input));
            }
            test_encode_decode_w_inputs(ext_inputs);
        }

        #[test]
        fn stress_test_encode_decode() {
            // TODO: Generate and test codes programmatically
        }
    }
}

/// Arithmetic code used for selector-based adaptivity
pub mod selector_arcd {
    use super::*;

    /// Multiply x by (k-1)/k
    fn mult_k_frac(x: u64, lg_k: usize) -> u64 {
        let k = 1 << lg_k;
        // [x/k * (k-1)] + [see if [leftover bits of x] * k-1 / k > 0]
        (x >> lg_k)*(k-1) + (x & (k-1)) - (x & (k-1) != 0) as u64
    }
    /// Multiply x by a-1/a where a is adapt rate
    fn mult_n_frac(x: u64) -> u64 {
        mult_k_frac(x, LG_ADAPTS)
    }
    /// Multiply x by (eps-1)/eps where eps is error rate
    fn mult_eps_frac(x: u64) -> u64 {
        mult_k_frac(x, LG_EPS)
    }
    /// Multiply x by Pr[letter], where Pr[letter] is defined as
    /// Pr[letter] = 7/8                            if letter = 0
    ///              1/8 * eps^(letter-1) * (1-eps) if letter > 0
    fn mult_pr(x: u64, letter: u64) -> u64 {
        if letter == 0 {
            mult_n_frac(x)
        } else {
            // (eps-1)/eps * [(x * 1/8) >> (lg(eps) * (letter-1))]
            // (eps-1)/eps * [(x * 1/8) >> lg(eps^(letter-1))]
            // (eps-1)/eps * [(x * 1/8) * eps^(letter-1)]
            mult_eps_frac(x >> LG_ADAPTS >> (LG_EPS * (letter-1) as usize))
        }
    }

    pub struct SelectorArcd;

    impl Arcd<u64, u64> for SelectorArcd {
        fn encode(input: [u64; 64]) -> Result<u64, EncodingFailure> {
            let mut low: u64 = 0;
            let mut high: u64 = !0 >> (64 - CODE_LEN);

            for i in 0..64 {
                let letter = input[i];
                let range = high - low;
                let first_range = mult_n_frac(range);

                if letter == 0 {
                    high = low + first_range;
                } else {
                    let mut total = 0;
                    for j in 0..letter-1 {
                        total += range >> (LG_EPS * j as usize) >> LG_ADAPTS;
                    }
                    let top = range >> (LG_EPS * (letter - 1) as usize);
                    low += first_range + mult_eps_frac(total);
                    high = low + (mult_eps_frac(top) >> LG_ADAPTS);
                }
                // Check if out of bits
                if high - low < 2 {
                    return Err(EncodingFailure)
                }
            }
            Ok(low)
        }
        fn decode(input: u64) -> [u64; 64] {
            let mut low: u64 = 0;
            let mut high: u64 = !0 >> (64 - CODE_LEN);
            let mut out = [0; 64];

            for i in 0..64 {
                let range = high - low;
                let first_range = mult_n_frac(range);
                if low + first_range > input {
                    high = low + first_range;
                } else {
                    let mut letter = 0;
                    let mut total = 0_u64;
                    let mut last_total = 0_u64;
                    while low + first_range + mult_eps_frac(total) <= input {
                        last_total = total;
                        total += range >> LG_EPS * letter >> LG_ADAPTS;
                        letter += 1;
                    }
                    let top = range >> (LG_EPS * (letter - 1));
                    out[i] = letter as u64;
                    low += first_range + mult_eps_frac(last_total);
                    high = low + (mult_eps_frac(top) >> LG_ADAPTS);
                }
            }
            out
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::cmp::max;

        /// Measures the number of bits needed to encode arr
        fn range_size(msg: [u64; 64]) -> Result<u64, usize> {
            // Get CODE_LEN ones
            let mut range: u64 = !0 >> (64 - CODE_LEN);
            // Iteratively reduce the range by scaling it
            // using the probability of each letter
            for i in 0..64 {
                range = mult_pr(range, msg[i]);
                if range == 0 {
                    return Err(i);
                }
            }
            Ok(range)
        }

        #[test]
        fn test_encode_decode() {
            let inputs = [
                "0000000000000000000000000000000000000000000000000000000000000000",
                "0000000000000000000000000000000000000000000000000000000000000001",
                "0000000000000000000000000000000000000000000000000000000000000002",
                "0000000000000000000000000000000000000000000000000000000000000003",
                "0000000000000000000000000000000000000020000000000000000000000002",
                "0000000000000000000000000000000000000020000000000000000000000003",
                "0000000000000000000000000000000000000000000000010000000000000000",
                "0000000000000000000000000000100000000100000100010000000000000000",
                "0000000000000011000000000000100000000100000100010000000000000000",
                "0000100000000011000000000000100000000100000100010100001000000110",
                "2000000000010001000000001101000002100000000000001000000000000101",
                "0000010000000000001101020001021100000000000100000000000000000011",
                "1111111111000000000000000000000000000000000000000000000000000000",
                "0101000000000000013000000000000001010100000100000001000000010101",
                "0001010000000000000001000001101000000000100300000000010010010001",
                "3000000000000100100010001000000010110000010000000001000000000000",
                "0000001000000000000110102000102110000000000100000000000000000011",

                // Ending 200000000 -> 133333333
                "1000010000000000002100000000000010110000000000100010001120000000",

                "1111111111111111111100000000000000000000000000000000000000000000", //20 1s: range=0, out of bits
                "1111111111111111111111111111111111111111111111111111111111111111", //range=0, out of bits

                // Moby Dick
                "0010001010111100100000000000000002000010010200000001001100012000", //range=0, out of bits

                // From chicagoA64 w/ insert-query ratio of 100:1
                "0000001000010101010000100010002010000100000100002000000100000000", //range=0, out of bits

                // From chicagoA64 w/ insert-query ratio of 1000:1
                "0100000011000001000100001000000001100000001000111000000100021000", //range=0 but passed
                "0100000001000001000100001000001001100000001000111000000100021000", //range=0 but passed
                "0000100000000000021010000000000101100000100001000100011130000000", //range=0, out of bits
                "0000010000001000010100002000011000100000000201020000000009000000", //range=0, out of bits
            ];

            /// Convert a test case string into an array of integer letters
            fn str_to_int_arr(s: &str) -> [u64; 64] {
                let mut out = [0_u64; 64];
                for (i, char) in s.chars().enumerate() {
                    match char.to_digit(10) {
                        Some(d) => out[i] = d as u64,
                        None => panic!("Example has a non-digit: {}", char),
                    }
                }
                out
            }

            // Convert inputs into int array
            let mut int_inputs = Vec::with_capacity(inputs.len());
            for input in inputs.iter() {
                int_inputs.push(str_to_int_arr(input));
            }
            // Check for each x in int_inputs that x = decode(encode(x))
            for input in int_inputs {
                match SelectorArcd::encode(input) {
                    // Encoding succeeds
                    Ok(code) => {
                        let decoded = SelectorArcd::decode(code);
                        assert_eq!(
                            input, decoded,
                            "x={:?}, encode(x)={}, decode(encode(x))={:?}",
                            input, code, decoded,
                        );
                    },
                    // Encoding fails (out of bits)
                    Err(_) => {
                        // Failed to encode but sequence encoding shouldn't overflow
                        // => error
                        if let Ok(size) = range_size(input) {
                            assert!(
                                size <= max(1, (64-CODE_LEN) as u64),
                                "Ran out of bits! (input={:?}, range_size={})",
                                input, size
                            );
                        }
                    }
                }
            }
        }
    }
}
