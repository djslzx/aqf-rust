# AQF in Rust
An implementation of the Adaptive Quotient Filter (AQF) in Rust.

&copy; David J. Lee, 2021.

## TODO
- Remote representation
  - [x] Add remote representation
  - [ ] Test remote rep (can't do yet until we start adapting)
- Arithmetic coding and fingerprint extensions
  - [x] Transcribe old arithmetic coding (tested)
  - [x] Write functions translating between fingerprint extension bits and the letters used in the arithmetic code (tested)
  - [x] Write new arithmetic coding (optimized for new probability distribution)
  - [ ] Integrate changes into `aqf` module
    - [ ] adapt a remainder
      - add/extend an extension until it doesn't clash with the elt causing a false pos
- Testing
  - [ ] Add integration tests (adversary, file-based, etc.)
- Organization
  - [ ] Make RSQF a supertrait of Filter

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

# Issues from C code

## Shifting selectors/extensions during inserts
In the C implementation, we don't shift selectors during insertions.  This isn't an issue for our current test suite, because we do all inserts before doing any queries.  This means that we don't have any nonzero selectors at insert time, so inserts need not worry about shifting selectors. But to make the filter more general, we need to shift selectors/extensions during insertions (do speed this up, we can check whether extensions are 0 and skip shifting them if this is the case).

## Rebuild logic
When encoding the extensions associated with a block fails, we zero out the extensions.  This requires that we update the remote representation as well, because the remote representation maps from (quot, rem, ext) triples to (word, hash) pairs.

When clearing a block, it is trivial to find its stored remainders and their extensions, as these are contained in the block. Finding the associated quotients is nontrivial, as the quotients contained in the block (in the form of the `occupieds` bitarray) don't necessarily correspond to the stored remainders.

So we need a way to figure out the quotient for each of the remainders stored in the block `B` in a filter `Q`.  Let `i` be the index of the block, so that its first remainder `B.rems[0]` is at position `i*64` in the filter. We refer to this position as `B.start`.

We already have a function that we can use to get from quotients to their corresponding runends; this is `rank_select(x)`, a shorthand for `rank(select(Q.runends, rank(Q.occupieds, x))`.  Note that this function only gives the corresponding runend if `x` is occupied; if it is not, then the result of `rank_select(x)` is the runend for the last occupied quotient before position `x`.

Because our filter is split into blocks, we store offsets in each block to let us run `rank_select` without traversing the whole filter. For the same reason, we can't just run a mirrored version of `rank_select` that maps runends to their quotients, because this would require that we store offsets that go in the opposite direction.

Instead, we have to make do with `rank_select`; that is, we need to figure out the quotients associated with the remainders in a block given 
- a function that maps from quotients to runends,
- the block's index, and
- block metadata (each block's `occupieds`, `runends`, and `offset`).

### Scanning forward/backward
There are two approaches I've thought about:

1. Determine the quotient of the first run that intersects with the block. (This requires scanning backwards through quotients until we find the last one whose runend >= `B.start`.) Then, scan backwards through the remainders in each run and step forward to the next occupied quotient. Repeat until we hit the last quotient whose run intersects the block. (This is roughly what we do in the C code.)
2. Determine the quotient of the last run that intersects the block. (This also requires scanning backwards through quotients as in case 1.) Then, scan backwards through the remainders in each run and step backward to the previous occupied quotient. Repeat until we hit the first quotient whose run intersects the block.

I've chosen to go with the latter approach because it gives us better locality and means that we don't have to do two passes.

### The algorithm
1. Find the quotient of the last run that intersects the block. This gives us a quotient, `q`, and `end`, the end of `q`'s run.
2. Step through each remainder in the run, clearing its extension in the remote representation. Stop when we reach a runend, the start of the block, or `q`.
3. Use rank and select to jump to the previous occupied quotient and runend. If the new runend is before the start of the block, stop. Otherwise, repeat Step 2 with the new quotient and runend.
4. Clear all extensions in the block. (This is faster than clearing each extension one by one.)

In pseudocode:
```
let B = Q.blocks[block_i]
let q, end = last_intersecting_run(block_i)
let i = end
loop quots:
  loop slots:
    update remote (q, Q.rem[i], Q.ext[i]) -> (q, Q.rem[i], none)
    if i == B.start:
      break quots
    else if B.runend(i-1) or i == q:
      break slots
    else:
      i -= 1
  q, end = prev_pair(q, end)

clear B.extensions
```


#### `last_intersecting_run(block_i)`
Start searching for the last intersecting run's quotient at the block's last quotient.  If the block has no quotients, then start at `B.start` or quit depending on whether the block's offset is zero. 
- If the offset is 0 and the block has no quotients, then the block is guaranteed to be empty, so we quit early.
- If the offset is positive and the block has no quotients, then the last intersecting run belongs to a quotient from a prior block.

After determining the location from which to start the search, jump backwards through previous (quotient, runend) pairs using rank and select to find the first pair whose runend < `B.start + 63`, and backtrack to the previously encountered (quotient, runend) pair. This is the pair we want.

```
// Get start position of search
let count = popcnt(B.occs)
if count == 0:
  if B.offset == 0: 
    return None
  else B.offset > 0:
    let q = prev_q(B.start)
    let end = B.offset - 1
    return Some(q, end)

else count > 0:
  let q = select(B.occs, count-1) 
  // select indexes from 0 - the first item is the 0th
  
  // Search backwards to find last intersecting run
  let end = rank_select(q)
  if end <= B.start + 63:
    return Some(q, end)
  else:
    loop:
      let last_q = q, last_end = end       // store current (q, end)
      match prev_pair(q, end, B.start):    // use rank-select
        Some(new_q, new_end) => 
          (q, end) = (new_q, new_end)
        None => return (last_q, last_end)
      if end < B.start + 63:
        return Some(last_q, last_end)
```

#### `prev_pair(q, end, bound)`
Use rank and select to get the previous (quot, runend) pair given (`q`, `end`). As an optimization, exit early (return None) if the previous pair's runend is before `bound` -- in this case, the caller does not need the pair.

```
let B_i = end/64
let B = blocks[B_i]
let count = bitrank(B.ends, end%64)
if count == 1:
  loop:
    B_i -= 1
    B = blocks[B_i]
    if B_i * 64 < bound:
      return None
    if B.ends != 0: 
      break

let new_end = B_i * 64 + highest_set_bit(B.ends)  
match prev_q(q) {
  Some(new_q) => return Some(new_q, new_end),
  None => panic!
```

#### `prev_q(q)`
Gets the last occupied quotient before `q`, returning `Some(loc)` if one exists and `None` if `q` is the first occupied quotient in the filter.

```
let block_i = q/64
let B = blocks[block_i]
let n = if q%64 > 0 { bitrank(B.occs, q%64-1) } else { 0 }
if n == 0:
  while block_i > 0:
    block_i -= 1
    B = blocks[block_i]
    if B.occs != 0:
      return block_i * 64 + highest_set_bit(B.occs) 
```

#### `highest_set_bit(bits)`
Gets the index of the MSB set to 1.
```rust
fn highest_set_bit(bits: u64) -> usize {
  assert_ne!(bits, 0);
  63 - bits.leading_zeros()
}
```

### Finding the quotients associated with the remainders stored in a block
Let's consider finding the quotients whose runs intersect the `i`-th block `B` in a filter `Q` where `B.start = i * 64`. Note that `rank_select(x)` is a shorthand for .

Idea: start with the last possible quotient that could point to a remainder in the block and work backwards.
