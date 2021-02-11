## TODO
- Testing
    - [ ] Add unit tests for `raw_insert` and `raw_query`, paying special attention to offset handling
    - [ ] Add integration tests (adversary, file-based, etc.)
- Organization
    - [ ] Make RSQF a super-trait of Filter

## To think about
- Only store the hash in remote and ditch the word (comparing word to elt is potentially slow, faster to compare hashes)
