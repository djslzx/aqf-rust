# AQF in Rust
An implementation of the Adaptive Quotient Filter (AQF) in Rust.

&copy; David J. Lee, 2021.

## Transitioning from C to Rust
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

## Offsets

In an RSQF, each block stores an offset, which lets us avoid searching the entire filter when looking for a run; instead, we can start at the block containing the quotient of the word we are interested in and use offsets to account for runs from earlier blocks.

The offset `O_i` of the absolute index `i` in filter `Q` is defined as
```
O_i = max(select(Q.runends, rank(Q.occupieds, i)) - i, 0)
```
This measures the distance from `i` to the `k`-th runend in `Q`, where `k` is the number of runs with quotients in the interval `[0, i]`.  Therefore, when `i` is occupied, `O_i` will be the distance from `i` to its runend; when `i` is not occupied, `O_i` will be the distance from `i` to the runend of the last quotient before `i`. In the second case, the distance may be negative (if the run before `i` ends before `i`); instead of storing this distance, we store `O_i = 0`, using `max(x, 0)` to clamp.

### Blocking

Because we only store offsets for the first slot in each block, for a block `B`,
- if `B[0]` is occupied, then `B.offset` will be the distance from `B.start` to the runend for `B[0]`;
- if `B[0]` is unoccupied, then `B.offset` will be either `0` or the distance from `B.start` to the runend for the last run from a block before `B`. 

To compute the offset of a slot `j` where `i = j - (j mod 64)`, i.e., `i` is the position of the start of the block that `j` belongs to, Pandey et al. suggest the following procedure:
```
d = rank(Q.occupieds[i+1, j], j-i-1)
O_j = select(Q.runends[i+O_i+1, end], d)
```
`d` represents the number of occupied quotients in the interval `[i+1, j]`; that is, between `i` and `j` and excluding `i`. `O_j` is computed as the `d`-th runend in `Q` after `i+0_i`.

### Insertions
Insertions will modify offsets in different ways depending on whether an inserted element's home slot is empty and whether the inserted elements extend existing runs or create new ones.

> #### We only increment non-negative offsets
> At insert time, we only increment those offsets that are already non-negative.  A negative offset will either remain negative or become zero, but because zero offsets and negative offsets are represented the same way (as zeros), we need not track this transition explicitly.

> #### Owned and unowned offsets
> If a block `B`'s offset points to the runend for `B[0]`, then we call `B.offset` _owned_.  If, instead, `B`'s offset points to the runend of the last runend in a block before `B`, then we say that `B.offset` is _unowned_.
>
> We can distinguish owned from unowned offsets by checking whether `B[0]` is occupied.
> If it is occupied, then `B.offset` is owned; otherwise, it is unowned.

Let `x` be the element we are inserting, with `q = quot(x)`. 

If the home slot `q` is empty, then there is no runend at `q`, so there can't be a non-negative offset pointing to it (because non-negative offsets must point to runends).  Therefore, we don't need to increment offsets in this case.

If the home slot `q` is taken, then we need to lookup the appropriate location to insert into.  First, we find the runend `r = rank_select(q)`.  

If `q` is occupied, then `r` marks the end of `q`'s run.  Inserting `x` means extending this run with `rem(x)`.  We first make room for `rem(x)` by finding `u`, the first unused slot after `r`, and shifting forward (by 1) the remainders and runends in the interval `[r+1, u-1]`, updating offsets with targets in this interval. Now that slot `r+1` is open, we insert `rem(x)` at `r+1`: we shift the runend bit at `r` to `r+1` and increment any non-negative offsets formerly pointing to `r`.

If `q` is unoccupied, then `r` marks the end of the last run before `q`. Here, inserting `x` means adding a new run.  As in the previous case, we make room for `rem(x)` at `r+1` by shifting runends and remainders and updating offsets.  Next, we handle offsets previously pointing to `r` in the following manner:
  - If the offset is owned, then don't increment it. Call the offset `o` and let `B` be the block it belongs to.  `o` being owned means that `r` is the index of the end of the run starting at `B.start`.  Adding a new run at `r+1` should not affect `B`'s offset.
  - If the offset is unowned, then:
    - If `q < B.start` for `B` the block containing `o`, then increment.  This condition must be met for `B.offset` to be affected by `q`'s runend.
    - If `q = B.start`, then increment.  Adding a new run for `B[0]` should set `B.offset` to `r+1`.  Because it was previously pointing to `r`, this means that incrementing `B.offset` will set it to the correct value.
    - If `q > B.start`, don't increment.  A run for a quotient after `B` should not affect `B.offset`.

## Handling fingerprint collisions

### Remote representation

Given elements `A` and `B` that have the same quotient `q` and remainder `r`, a non-adaptive quotient filter need only store one copy of the pair `(q, r)`.  

In the AQF, we store both `A` and `B` in the remote representation, so that the fingerprints for both can adapt independent of each other. This leaves us with a choice of what to store in the local filter:

1. "One-to-many": We store one copy of the pair `(q, r)` in the filter, which maps to `A` and `B` in the remote representation.  When there is a collision at `(q, r)` for a query `C`, we adapt either one or all of the elements associated with it, which requires that we perform an insertion to handle the case where the extensions for `A` and `B` (relative to `C`) diverge.  

2. "Many-to-many": We store two copies of the pair `(q, r)` in the filter. This means that we need not perform additional inserts when a query collides with `(q, r)`.

We go with the latter approach because it doesn't require additional insertions at adapt time.

### Decision points

There are a few key points in the code where we need to choose policies to deal with fingerprint collisions consistently.

At query time, we need policies to determine the following:
 - _Whether a query fingerprint matches a stored fingerprint._ Our current policy is to find the first stored fingerprint whose quotient, remainder, and extension match that of the query fingerprint. This means that if stored fingerprints `f(a)` and `f(b)` share the same quotient and remainder, but `f(a)` has an extension and `f(b)` doesn't, such that `f(a) = q:r:e` and `f(b) = q:r:_`, then when querying with `a`, we may find `q:r:_` instead of `q:r:e`. 

 - _Which remote elements are associated with a local element._ Once we have found a stored fingerprint `f(s)` that matches a query fingerprint `f(q)` of query element `q`, we need to determine whether `q` is in the filter (i.e., whether the result is a true or false positive) to determine whether we need to adapt.  This requires retrieving `S(f(s))`, the set of stored elements associated with `f(s)`, and seeing whether any `s` in `S(f(s))` matches `q`. Therefore, we need to define `S`, which maps from a stored fingerprint to its associated stored elements.

At adapt time, we need to determine the following:
 - _Which remote element to update._ If `q` does not match any element in `S(f(s))`, then we need to adapt.  But if there are multiple elements in `S(f(s))`, which of these elements should we extend the fingerprints of?
   - One: extend the fingerprint of a single element plucked arbitrarily from `S(f(s))`.
     - Pros: fast, rebuilds less often
     - Cons: under-adapts
   - All: extend the fingerprints for all `s` in `S(f(s))`
     - Pros: fixes collisions more thoroughly
     - Cons: slow;  over-adapts; only distinguishes false matches from query, but not with each other 
 - _Which remote elements to rebuild._ If we run out of space in a block's encoding, we clear all the extensions of local fingerprints in the block. This necessitates that we also clear the extensions of the remote elements, which in turn requires that we have a mapping from local fingerprints to remote elements.  

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

1. Determine the quotient of the first run that intersects with the block. (This requires scanning backwards through quotients until we find the last one whose runend >= `B.start`.) Then, walk backwards through the remainders in each run, but step forward to the next occupied quotient when done going through a given run. Repeat until we hit the last quotient whose run intersects the block. (This is roughly what we do in the C code.)
2. Determine the quotient of the last run that intersects the block. (This requires scanning backwards through quotients until we find the first one whose run doesn't overshoot the end of the block.) Then, walk backwards through the remainders in each run, and when done with a run, step backward to the previous occupied quotient. Repeat until we hit the first quotient whose run intersects the block.

I've chosen to go with the latter approach because it gives us better locality and means that we don't have to do two passes.

### The algorithm
1. Find the quotient (`q`) and runend (`end`) of the last run that intersects the block.
2. Step through each remainder in the run, clearing its extension in the remote representation if the remainder is in the block of interest. Stop when we reach a runend, the start of the block, or `q`.
3. Use rank and select to jump to the previous occupied quotient and runend. If the new runend is before the start of the block, stop. Otherwise, repeat Step 2 with the new quotient and runend.
4. Clear all extensions in the block. (This is faster than clearing each extension one by one.)

In pseudocode, where `Q` is the filter, `B` is the block of interest, and `block_i` is the index of the block: 
```
let B = Q.blocks[block_i]
let (q, i) = last_intersecting_run(block_i)

// skip backward until we hit the block 
while i >= block_start + 64:
  i -= 1
  if Q.runend(i):
    q = prev_quot(q)
    // q is a quotient that is in block B or any earlier block.
    // therefore, if it is not intersecting B, that means that it
    // must have been pushed over by another run; 
    // so, if there's no runend between the end of B and q's runend,
    // then the run for q overlaps with B.
  loop block:
    loop run:
      update remote (q, Q.rem[i], Q.ext[i]) -> (q, Q.rem[i], none)
      if i == B.start:
        break block
      else if i == q or Q.runend(i-1):
        break run
      else
        i -= 1
    (q, i) = prev_pair(q, i)

clear B.extensions
```

### Helper functions

#### `last_intersecting_run(block_i)`
Determines the quotient and runend for the last run that intersects the `block_i`-th block.

Start searching for the last intersecting run's quotient at `B`'s last quotient.  If `B` has no quotients, then start at `B.start := block_i * 64` or quit depending on whether `B`'s offset is zero. 
- If the offset is 0 and the block has no quotients, then the block is guaranteed to be empty, so we quit early.
- If the offset is positive and the block has no quotients, then the last intersecting run belongs to a quotient from a prior block.

If `B` has occupied quotients, then we start at the last occupied quotient. (We can do this by using either `leading_zeros` or `bitselect(B.occs, popcnt(B.occs)-1)`.)

After determining the location from which to start the search, jump backwards through previous (quotient, runend) pairs using rank and select to find the first pair whose runend < `B.start + 63`, and backtrack to the previously encountered (quotient, runend) pair. This is the pair we want.

```
// (1) Get start position of search
let n = popcnt(B.occs)
if n == 0:
  if B.offset == 0 and B[0] is not a runend: 
    return None
  else:
    let q = prev_q(B.start)
else:
  let q = highest_set_bit(B.occs) // get last occupied quotient 
  //      ^ defined below, same as select(B.occs, count-1) in this context
  
// (2) Search backwards to find last intersecting run
let end = rank_select(q)
if end <= B.start + 63:
  return Some(q, end)
else:
  loop:
    let last_q = q, last_end = end     // store current (q, end)
    if prev_pair(q, end) gives None:    // prev_pair uses rank-select and 
                                        // returns None if no valid result found
      return last_q, last_end
    else:
      q, end = prev_pair(q, end)
    
    if end < B.start + 63:
      return last_q, last_end
```

#### `prev_pair(q, end, bound)`
Uses rank and select to get the previous (quot, runend) pair given (`q`, `end`). As an optimization, exits early (returns None) if the previous pair's runend is before `bound`.

This can be implemented straightforwardly by defining the helpers `prev_end(end, bound)` and `prev_q(q)`.

#### `prev_end(end, bound_i)`
```
let block_i = end/64
let b = Q.blocks[block_i]
let n = if end%64 == 0 then 0, else rank(b.runends, end%64 - 1)
if n == 0:
  for i in [block_i-1, bound_i]:
    b = Q.blocks[i]
    if b.runends > 0:
      return i*64 + highest_set_bit(b.runends)
  return None
else:
  return block_i * 64 + select(b.runends, n-1)
```

#### `prev_q(q)`
Gets the last occupied quotient (`q'`) before `q`, returning `Some(q')` if one exists and `None` if `q` is the first occupied quotient in the filter.

```
let block_i = q/64
let b = Q.blocks[block_i]
let n = if q%64 == 0 then 0, else rank(b.occupieds, q%64 - 1)
if n == 0:
  for i in [block_i-1, 0]:
    let b = Q.blocks[i]
    if b.occupieds > 0:
      return i*64 + highest_set_bit(b.occupieds)
  return None
else:
  return i*64 + select(b.occupieds, n-1)
```

#### `highest_set_bit(bits)`
Gets the index of the MSB set to 1.
```rust
fn highest_set_bit(bits: u64) -> usize {
  assert_ne!(bits, 0);
  63 - bits.leading_zeros()
}
```