# metaRE news

## Unreleased (chore/modernize-package)

### Bug fixes

- **Fix incorrect motif matching on Windows / any LLP64 platform.** Several bit
  masks were built by shifting a 32-bit literal (`1L`, `1`, `(unsigned long)1`)
  by up to 60–62 positions. On LLP64 platforms (Windows), `long`/`int` are
  32-bit, so these shifts are undefined behaviour and produced corrupted masks.
  Symptoms:
  - `enumeratePatterns()` on long, degenerate (N-containing) patterns whose
    length mod 16 was ≥ 8 matched almost every sequence (~96%) instead of the
    true frequency, because the last 64-bit cell of `IUPACMotif::includes()` was
    mis-masked (`src/Motifs/IUPACMotif.cpp`).
  - `Pattern` stored its k-mer bitset in `unsigned long` words while indexing
    bits 0–63 (`MASK_SHIFT=6`, `MASK_FILTER=63`); on Windows bits 32–63 aliased
    onto 0–31, causing false-positive core k-mers in the spaced-dyad / repeat
    counters (`src/Pattern/Pattern.{h,cpp}`).
  - Reverse-complement / last-cell masks in the motif builders had the same
    32-bit-shift defect (`src/Motifs/IUPACMotifBuilder.cpp`,
    `src/Motifs/CompactMotifBuilder.cpp`).

  All are fixed by performing the shifts in 64-bit (`(cell)1`, `(uint64_t)1`,
  `~(cell)0`). Added regression tests in
  `tests/testthat/test-enumerateMotifs-longpattern.R` covering pattern lengths
  that straddle the 16-nt cell boundary. After the fix, **forward-strand**
  matching by `enumeratePatterns()` and `enumerateDyadsWithCore()` agrees
  exactly with an independent regex count on real promoters.

### Known remaining issue

- **Reverse-complement matching of long motifs (`rc = TRUE`).** For
  non-palindromic motifs longer than 16 nt (spanning more than one IUPAC cell),
  `enumeratePatterns(..., rc = TRUE)` still misses occurrences on the reverse
  strand — the multi-cell reverse-complement window is not assembled correctly.
  Palindromic patterns and motifs ≤ 16 nt are unaffected. The whole-cell-shift
  UB in `fillComplementBuf` was fixed, but a further defect remains in the
  reverse-complement bookkeeping and is not yet located. **Workaround:** run
  `enumeratePatterns(seq, P, rc = FALSE)` and, separately, on the reverse
  complement of the sequences, then take the union (this uses only the
  verified-correct forward matcher). The regression test for this case is
  `skip()`ped and references this note.
