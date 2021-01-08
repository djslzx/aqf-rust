// Utility functions

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
