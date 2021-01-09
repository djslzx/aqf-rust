// Utility functions

// 64-bit Bit Array
#[derive(Debug)]
pub struct B64(u64);

impl B64 {
    pub fn zero() -> B64 {
        B64(0)
    }
    pub fn get(&self, at: usize) -> bool {
        (self.0 & (1 << at)) == 0
    }
    pub fn set(&self, at: usize) -> B64 {
        B64(self.0 | (1 << at))
    }
    pub fn unset(&self, at: usize) -> B64 {
        B64(self.0 & !(1 << at))
    }
    pub fn set_to(&mut self, on: bool, at: usize) {
        *self = if on {
            self.set(at)
        } else {
            self.unset(at)
        };
    }
}

// Round v to nearest power of 2
// Pre: v >= 0
// https://graphics.stanford.edu/~seander/bithacks.html
pub fn nearest_pow_of_2(mut v: usize) -> usize {
    v = v-1;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v + 1
}
