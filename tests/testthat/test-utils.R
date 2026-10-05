test_that("rel100 normalizes each sample to 100 percent", {
  m <- matrix(c(10, 5, 1, 2), nrow = 2,
              dimnames = list(c("a", "b"), c("s1", "s2")))
  r <- phyloseqExplorer:::rel100(m)
  expect_equal(as.numeric(colSums(r)), c(100, 100), tolerance = 1e-6)
  expect_equal(r["a", "s1"], 10 / 15 * 100, tolerance = 1e-6)
})

test_that("natural_key zero-pads numbers for natural sorting", {
  expect_true(phyloseqExplorer:::natural_key("ASV10") > phyloseqExplorer:::natural_key("ASV2"))
  expect_equal(phyloseqExplorer:::natural_key("abc"), "abc")
})

test_that("stars thresholds are correct", {
  expect_equal(phyloseqExplorer:::stars(c(0.0005, 0.005, 0.03, 0.2, NA)),
               c("***", "**", "*", "ns", ""))
})

test_that("taxon_cols always greys Others", {
  c2 <- phyloseqExplorer:::taxon_cols(c("A", "Others"))
  expect_equal(unname(c2["Others"]), "#B0B0B0")
  expect_equal(names(c2), c("A", "Others"))
})

test_that("stars_blank drops ns", {
  expect_equal(phyloseqExplorer:::stars_blank(c(0.2, 0.009)), c("", "**"))
})

test_that("make_palette recycles and extends", {
  p30 <- phyloseqExplorer:::make_palette(5, "contrast")
  expect_length(p30, 5)
  p_long <- phyloseqExplorer:::make_palette(100, "set1_60")
  expect_length(p_long, 100)
  expect_false(any(duplicated(p_long[1:60])))
})
