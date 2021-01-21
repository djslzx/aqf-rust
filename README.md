# AQF in Rust
An implementation of the Adaptive Quotient Filter (AQF) in Rust.

&copy; David J. Lee, 2021.

## TODO
- Remote representation
  - [x] Add remote representation
  - [ ] Test remote rep (can't do yet until we start adapting)
- Arithmetic coding and fingerprint extensions
  - [x] Transcribe old arithmetic coding
    - [ ] Test
  - [x] Write functions translating between fingerprint extension bits and the letters used in the arithmetic code
    - [x] Test 
  - [ ] Write new arithmetic coding (optimized for new probability distribution)
  - [ ] Integrate changes into `aqf` module
    - [ ] Adapt on lookups
- Testing
  - [ ] Add integration tests (adversary, file-based, etc.)

## Implementation notes
The logic is very close to that of the C implementation, with a few notable differences:

- *The filter now grows when it runs out of space for remainders, instead of starting off
  with a static amount of extra space.* We do this by putting the filter's blocks in a 
  vector instead of an array, which could impact total space utilization but affords us
  greater flexibility and generality.
- `find_runend` is renamed as `rank_select`, as this better represents what we are doing:
  this function does not find the runend of a particular slot; it is merely a useful shorthand for 
  a blocked operation that reused multiple times (`select(Q.runends, rank(Q.occupieds, x))`), notably in insertions, queries, and the `find_first_unused`
  operation.
- `rank_select(x)` makes use of a custom Enum called `RankSelectResult`, 
  which consists of three cases (`Empty`, `Full(loc)`, and `Overflow`).
  This allows the function to more clearly signal its result:
  - `Empty` indicates that the result of `select(Q.runends, rank(Q.occupieds, x))`
    precedes `x`; i.e., the home slot for the quotient `x` is open.
  - `Full(loc)` indicates that the home slot for `x` is taken, and `loc`
    is the last slot that is taken by a run for a quotient smaller than `x`.
  - `Overflow` indicates that the search passed the last block in the filter.
- Elsewhere in the code, we replace haphazard flag values and returns by reference
  with proper return values, sometimes using `Option`s.
  Generally, `Some(y)` indicates a valid result `y` and `None` indicates an invalid result
  (e.g., running out of space or going off the end of an array).
- Loop logic is modified somewhat; Rust lacks `do-while` loops, so we make do with
  infinite loops that break when a condition is met.
- Casts between integer types are made explicit. 
- Boolean to integer casts are implemented using conditional expressions or explicit casts instead of implicit casts.
- Integer overflows are fixed. Whereas C allows them, Rust does not -- these overflows cause compiler errors. So the Rust AQF prevents these overflows from happening with additional checks. 
  This is especially important for incrementing direct/indirect offsets.
- Functions are more rigorously tested. See the `tests` module for details.
