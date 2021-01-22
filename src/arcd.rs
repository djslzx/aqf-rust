// Arithmetic coding

// Adaptivity bits are represented as an Option<(u64, usize)>
// where the u64 holds the bits and the usize is the number of bits

const CODE_LEN: usize = 56;     // length of arithmetic code
const LG_ADAPTS: usize = 2;     // lg(adapt_rate)
const LG_EPS: usize = 4;        // -lg(eps)

#[derive(Debug)]
enum AdaptBits {
    None,
    Some {bits: u64, len: usize},
}

/// Computes floor(log2(x)) by counting the number of leading 0s
/// to get the position of the first leading 1
fn floor_log(x: u64) -> u32 {
    if x == 0 { 0 } // treat log(0) as 0
    else { 63 - x.leading_zeros() }
    // this uses a compiler intrinsic, which i'm assuming is fast
}

/// Converts a letter to adapt bits
///   len = floor(log2(letter + 1))
///   bits = letter - (2^len - 1)
fn to_bits(letter: u64) -> AdaptBits {
    let len = floor_log(letter+1) as usize;
    if len == 0 {
        AdaptBits::None
    } else {
        AdaptBits::Some {
            bits: letter - ((1 << len) - 1),
            len: len as usize
        }
    }
}

/// Converts adaptivity bits to a letter
///   letter = 0                   if len=0
///          = (2^len - 1) + bits  if len>0
fn to_letter(bits: &AdaptBits) -> u64 {
    match *bits {
        AdaptBits::None => 0,
        AdaptBits::Some {bits, len} => {
            assert_ne!(len, 0, "Length-less adaptivity bits should be None");
            ((1 << len)-1) as u64 + bits
        }
    }
}

fn slow_factorial(x: u64) -> u64 {
    let mut out = 1;
    for i in 2..=x {
        out *= i;
    }
    out
}
/// Multiply x by Pr[k]
/// Pr[k] = 1/(2^k * k!e)
// FIXME: probability function is incorrect
fn mult_pr(x: u64, k: usize) -> u64 {
    let k_factorial = slow_factorial(k as u64) as f64;
    (((x >> k) >> LG_EPS) as f64/k_factorial) as u64
    // FIXME: slow b/c of factorial
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
    #[test]
    fn test_factorial() {
        let mut x = 1;
        for i in 0..20 {
            x = if i == 0 { 1 } else { x * i };
            assert_eq!(slow_factorial(i), x);
        }
    }
    #[test]
    fn test_conversions() {
        for i in 0..1000 {
            let bits = to_bits(i);
            let letter = to_letter(&bits);
            assert_eq!(i, letter,
                       "l={} -> bits={:?} -> l={}",
                       i, bits, letter);
        }
    }
    //FIXME: probability function is incorrect
    #[test]
    fn test_mult_pr() {
        fn slow_mult_pr(x: u64, k: usize) -> u64 {
            let epsilon = (1 << LG_EPS) as f64;
            ((x as f64)/(2_u64.pow(k as u32) as f64
                * slow_factorial(k as u64) as f64
                * epsilon)) as u64
        }
        for i in 0..100 {
            for j in 0..10 {
                let fast = mult_pr(i, j);
                let slow = slow_mult_pr(i, j);
                assert_eq!(fast, slow,
                           "i={}, j={}, fast={}, slow={}",
                           i, j, fast, slow);
            }
        }
    }
}

/// Arithmetic code
/// Define helper functions that encode/decode need; as long as a type
/// implements the helpers, it'll have access to encode/decode
trait Arcd {
    // TODO: implement
    fn encode(input: [u64; 64]) -> Option<u64>;
    // TODO: implement
    fn decode(input: u64) -> [u64; 64];
}

// For reference
mod old_arcd {
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
    fn encode(input: [u64; 64]) -> Option<u64> {
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
                return None
            }
        }
        Some(low)
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
                match encode(input) {
                    // Encoding succeeds
                    Some(code) => {
                        let decoded = decode(code);
                        assert_eq!(
                            input, decoded,
                            "x={:?}, encode(x)={}, decode(encode(x))={:?}",
                            input, code, decoded,
                        );
                    },
                    // Encoding fails (out of bits)
                    None => {
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