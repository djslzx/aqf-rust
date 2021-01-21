// Arithmetic coding

// Adaptivity bits are represented as an Option<(u64, usize)>
// where the u64 holds the bits and the usize is the number of bits

enum AdaptBits {
    None,
    Some {bits: u64, len: usize},
}

/// Converts a letter to adapt bits
///   len = floor(log2(letter + 1))
///   bits = letter - (2^len - 1)
fn to_bits(letter: u64) -> Option<(u64, usize)> {
    let l = letter as f32;
    let len = (l+1.0).log2().floor() as usize;
    if len == 0 {
        None
    } else {
        let bits = letter - ((1 << len) - 1);
        Some((bits, len as usize))
    }
}

/// Converts adaptivity bits to a letter
///   letter = 0                   if len=0
///          = (2^len - 1) + bits  if len>0
fn to_letter(bits: Option<(u64, usize)>) -> u64 {
    match bits {
        None => 0,
        Some((bits, len)) => {
            assert_ne!(len, 0, "Lengthless bits should be None instead");
            ((1 << len)-1) as u64 + bits
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conversions() {
        for i in 0..1000 {
            let bits = to_bits(i);
            let letter = to_letter(bits);
            assert_eq!(i, letter,
                       "l={} -> bits={:?} -> l={}",
                       i, bits, letter);
        }
    }
}

// const CODE_LEN: usize = 56;     // length of arithmetic code
// const LOG_ADAPT: usize = 2;     // lg(adapt_rate)
// const LOG_EPS: usize = 4;       // lg(eps)

// /// Multiply x by (k-1)/k
// fn mult_k_frac(x: u64, lg_k: usize) -> u64 {
//     let k = 1 << lg_k;
//     // [x/k * (k-1)] + [leftover lg_k bits of x] - [1 if leftover bits nonzero, else 0]
//     (x >> lg_k)*(k-1) + (x & (k-1)) - (x & (k-1) != 0)
// }
// /// Multiply x by a-1/a where a is adapt rate
// fn mult_n_frac(x: u64) -> u64 {
//     mult_k_frac(x, LOG_ADAPT)
// }
// /// Multiply x by (eps-1)/eps where eps is error rate
// fn mult_eps_frac(x: u64) -> u64 {
//     mult_k_frac(x, LOG_EPS)
// }
// /// Multiply x by Pr[letter], where Pr[letter] is defined as 
// /// Pr[letter] = 7/8                            if letter = 0
// ///              1/8 * eps^(letter-1) * (1-eps) if letter > 0
// fn mult_pr(x: u64, letter: u32) -> u64 {
//     if letter == 0 {
//         mult_n_frac(x)
//     } else {
//         // (eps-1)/eps * [(x * 1/8) >> (lg(eps) * (letter-1))]
//         // (eps-1)/eps * [(x * 1/8) >> lg(eps^(letter-1))]
//         // (eps-1)/eps * [(x * 1/8) * eps^(letter-1)]
//         mult_eps_frac(x >> LOG_ADAPT >> (LOG_EPS * (letter-1)))
//     }
// }

// fn encode(input: [u32; 64]) -> Option<u64> {
//     let mut low: u64 = 0;
//     let mut high: u64 = !0 >> (64 - CODE_LEN);
    
//     for i in 0..64 {
//         let letter = input[i];
//         let range = high - low;
//         let first_range = mult_n_frac(range);
        
//         if letter == 0 {
//             high = low + first_range;
//         } else {
//             let mut total = 0;
//             for j in 0..letter-1 {
//                 total += (range >> (LOG_E * j) >> LOG_ADAPT_RATE);
//             }
//             let top = range >> (LOG_E * (letter - 1));
//             low += first_range + mult_eps_frac(total);
//             high = low + (mult_eps_frac(top) >> LOG_ADAPT_RATE);
//         }
//         // Check if out of bits
//         if high - low < 2 {
//             return None
//         }
//     }
//     Some(low)
// }

// fn decode(input: u64) -> [u32; 64] {
    
// }
