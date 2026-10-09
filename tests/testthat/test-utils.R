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

test_that("cooc_stats derives r and p from the same matrix", {
  # four dominant 18S species, real counts, uneven library sizes. Converting to
  # relative abundance reorders samples, so the counts-based and composition-based
  # Spearman disagree here - for t1 ~ t3 they differ in sign (+0.3 vs -0.6).
  cnt <- structure(c(48943L, 21602L, 2191L, 345L, 30156L, 15970L, 2697L,
                     196L, 55709L, 13936L, 21541L, 107L, 20818L, 1789L, 1587L,
                     149L, 10784L, 8676L, 3341L, 4828L), dim = 4:5,
                   dimnames = list(c("t1", "t2", "t3", "t4"), paste0("s", 1:5)))
  cs  <- phyloseqExplorer:::cooc_stats(cnt, "spearman")
  rel <- phyloseqExplorer:::rel100(cnt)

  # r is the composition-based coefficient ...
  expect_equal(cs$r["t1", "t3"], cor(rel["t1", ], rel["t3", ], method = "spearman"),
               tolerance = 1e-9)
  # ... and it really does differ from the counts-based one, so this test has teeth
  expect_false(isTRUE(all.equal(cs$r["t1", "t3"],
                                cor(cnt["t1", ], cnt["t3", ], method = "spearman"))))

  # every p must be the p-value of the r reported beside it
  n <- ncol(cnt)
  for (pr in list(c("t1", "t3"), c("t1", "t2"), c("t2", "t4"), c("t3", "t4"))) {
    r <- cs$r[pr[1], pr[2]]
    expect_equal(cs$p[pr[1], pr[2]],
                 2 * pt(abs(r) / sqrt((1 - r^2) / (n - 2)), n - 2, lower.tail = FALSE),
                 tolerance = 1e-8)
  }
})

test_that("cooc_stats BH correction uses only the upper triangle", {
  set.seed(1)
  cnt <- matrix(rpois(30, 100), nrow = 6, dimnames = list(letters[1:6], paste0("s", 1:5)))
  cs <- phyloseqExplorer:::cooc_stats(cnt, "spearman")
  k <- sum(!is.na(cs$p))
  expect_equal(k, choose(6, 2))
  expect_equal(cs$p_adj[!is.na(cs$p)], p.adjust(cs$p[!is.na(cs$p)], "BH"), tolerance = 1e-12)
})

test_that("node_tax_lookup resolves labels at the requested rank", {
  tt <- cbind(Kingdom = c("Bacteria", "Bacteria"), Genus = c("HIMB11", "HIMB11"))
  rownames(tt) <- c("ASV1", "ASV7")
  out <- phyloseqExplorer:::node_tax_lookup(tt, "Genus", "HIMB11")
  expect_equal(out$Taxon, "HIMB11")
  expect_equal(out$Kingdom, "Bacteria")
  # a label with no match at that rank yields NA rather than a wrong lineage
  miss <- phyloseqExplorer:::node_tax_lookup(tt, "Genus", "Karenia")
  expect_true(is.na(miss$Kingdom))
})

test_that("node_tax_lookup does not leak dataset-1 lineage onto a colliding label", {
  # dataset 1 (16S), node labels are genera; "ASV1" here is a 16S ASV id
  tt1 <- cbind(Kingdom = "Bacteria", Genus = "HIMB11"); rownames(tt1) <- "ASV1"
  # dataset 2 (18S) contributes a node literally named "ASV1"
  tt2 <- cbind(Domain = "Eukaryota", Species = "Karenia_selliformis"); rownames(tt2) <- "ASV1"
  d1 <- phyloseqExplorer:::node_tax_lookup(tt1, "Genus", character(0))
  d2 <- phyloseqExplorer:::node_tax_lookup(tt2, "(ASV)", "ASV1")
  both <- phyloseqExplorer:::rbind_fill(d1, d2)
  expect_equal(both$Species, "Karenia_selliformis")
  expect_false("Bacteria" %in% unlist(both))
})

test_that("rbind_fill unions columns and fills gaps with NA", {
  a <- data.frame(Taxon = "x", Domain = "Eukaryota", stringsAsFactors = FALSE)
  b <- data.frame(Taxon = "y", Kingdom = "Bacteria", stringsAsFactors = FALSE)
  out <- phyloseqExplorer:::rbind_fill(a, b)
  expect_equal(sort(names(out)), c("Domain", "Kingdom", "Taxon"))
  expect_true(is.na(out$Kingdom[out$Taxon == "x"]))
  expect_true(is.na(out$Domain[out$Taxon == "y"]))
})
