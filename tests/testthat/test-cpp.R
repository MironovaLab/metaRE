test_that("Catch unit tests pass", {
    skip("Legacy native testthat runner is not yet wired up on the modern toolchain")
    expect_cpp_tests_pass("metaRE")
})
