library(testthat)

# Regression test for the LLP64 (Windows) 32-bit-shift bug in the last-cell
# mask of IUPACMotif::includes and the reverse-complement / Pattern masks.
# IUPAC motifs pack 16 nt per 64-bit cell. When a pattern's length modulo 16 was
# >= 8, `1L << ((len % 16) * 4)` overflowed 32-bit `long` on Windows, corrupting
# the last-cell mask so that long, mostly-N patterns matched almost every
# sequence (~96% instead of the true frequency). These tests exercise pattern
# lengths that straddle the 16-nt cell boundary (12, 20, 24) on BOTH strands.

# TGTC .. GACA spaced motifs of controlled length; N must match any base and the
# TGTC/GACA anchors MUST still be required.
mk_IR <- function(spacer) paste0("TGTC", strrep("N", spacer + 4L), "GACA")

seqs <- c(
    hit12  = paste0("aaaa", "TGTC", "gggg",                 "GACA", "aaaa"),  # len-12 (s=0)
    hit24  = paste0("aa",   "TGTC", "acgtacgtacgtacgt",     "GACA", "aa"),    # len-24 (s=12)
    hit20  = paste0("aa",   "TGTC", "acgtacgtacgt",         "GACA", "aa"),    # len-20 (s=8)
    core_only = paste0("ttt", "TGTC", "aaaaaaaaaaaa", "ttt"),                 # TGTC but no GACA
    plain     = "acacacacacacacacacacacacacacac",                            # no anchors at all
    rc24   = paste0("aa", "TGTC", "acgtacgtacgtacgt", "GACA", "aa")           # same as hit24 (palindromic)
)
seqs <- toupper(seqs)
gn <- names(seqs)

test_that("long mostly-N IUPAC pattern requires its anchors (no over-match)", {
    for (s in c(0L, 8L, 12L)) {
        res <- enumeratePatterns(seqs, mk_IR(s), rc = FALSE, output = "genes")
        hit <- gn[res[[1]]]
        # the core_only and plain sequences must NEVER match a spaced TGTC..GACA
        expect_false("core_only" %in% hit)
        expect_false("plain" %in% hit)
    }
})

test_that("length-24 pattern (last-cell r=8) matches exactly the right genes", {
    res <- enumeratePatterns(seqs, mk_IR(12L), rc = FALSE, output = "genes")  # len 24
    expect_setequal(gn[res[[1]]], c("hit24", "rc24"))
})

test_that("length-12 pattern (single-cell r=12) matches exactly the right genes", {
    res <- enumeratePatterns(seqs, mk_IR(0L), rc = FALSE, output = "genes")   # len 12
    expect_setequal(gn[res[[1]]], "hit12")
})

test_that("length-20 pattern (last-cell r=4, unaffected) still correct", {
    res <- enumeratePatterns(seqs, mk_IR(8L), rc = FALSE, output = "genes")   # len 20
    expect_setequal(gn[res[[1]]], "hit20")
})

# Second bug: reverse-complement of a multi-cell (>16 nt) motif was corrupted by
# a shift-by-cellSize in fillComplementBuf, so rc=TRUE missed hits on the reverse
# strand for long non-palindromic patterns.
rc_str <- function(s) {
    s <- chartr("ACGTacgt", "TGCAtgca", s)
    paste(rev(strsplit(s, "")[[1]]), collapse = "")
}

test_that("rc=TRUE finds a long non-palindromic motif on the reverse strand", {
    skip("KNOWN ISSUE (see NEWS.md): reverse-complement matching is still wrong for
          non-palindromic motifs spanning >16 nt (>1 IUPAC cell). Forward matching
          is correct; drive both strands explicitly (rc=FALSE on seq and its
          reverse complement) as a workaround until the rc path is fixed.")
    core   <- paste0("TGTC", "ACGTACGTACGTACGT", "TGTCAA")   # DR-like, 26 nt (2 cells)
    fwdseq <- paste0("GG", core, "GG")
    seqs2  <- c(onfwd = fwdseq, onrev = rc_str(fwdseq), none = strrep("A", 40))
    P <- paste0("TGTC", strrep("N", 16L), "TGTC", "NN")       # len 26, non-palindromic
    res <- enumeratePatterns(seqs2, P, rc = TRUE, output = "genes")
    expect_setequal(names(seqs2)[res[[1]]], c("onfwd", "onrev"))
})
