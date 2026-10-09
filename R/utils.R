# ---------------------------------------------------------------------
# Internal helpers: palettes, encoding repair, data access, exporters
# ---------------------------------------------------------------------
has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

options(shiny.maxRequestSize = 500 * 1024^2)
suppressWarnings(try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE))

# ---------------------------------------------------------------------
# PALETTES (taken from your script)
# ---------------------------------------------------------------------
PALETTE_60 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
  "#1b4b32", "#fec4df", "#1cebf5", "#a38a00", "#840d35", "#aec4f1",
  "#00948b", "#524000", "#0037a5", "#f1cfa7", "#fd3e85", "#497100",
  "#ef6efe", "#c68d82", "#b0e0bd", "#006170", "#9e8bb6", "#4e3a5f",
  "#00b8f6", "#7fee1d", "#965d62", "#896aff", "#613726", "#f4af02",
  "#bdc16a", "#7f885c", "#0a7d4d", "#7300b3", "#9b0000", "#686a94",
  "#816a3b", "#c40970", "#00c19d", "#b2976f", "#02f298", "#b97400",
  "#ff6d4e", "#dd89b8", "#820063", "#dbb5fb", "#60afba", "#859e00",
  "#054e00", "#555a31", "#fabeaf", "#f600be", "#6d99fa", "#0f889f",
  "#9bdefb", "#cd0040", "#87ad88", "#005a98", "#374ef9")

PALETTE_30 <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFD92F","#A65628",
  "#F781BF","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#1B9E77","#D95F02",
  "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#1F78B4","#B2DF8A",
  "#33A02C","#FB9A99","#CAB2D6","#6A3D9A","#FFFF99","#B3B3B3","#FDBF6F",
  "#FF7F00","#FBB4AE")

make_palette <- function(n, which = "set1_60") {
  if (n <= 0) return(character(0))
  base <- switch(which,
                 set1_60  = PALETTE_60,
                 contrast = PALETTE_30,
                 NULL)
  if (is.null(base)) {
    return(if (which == "viridis") scales::viridis_pal()(n) else scales::hue_pal()(n))
  }
  if (n > length(base)) return(colorRampPalette(base)(n))
  base[seq_len(n)]
}

# colors for taxa levels: "Others" is always #B0B0B0 (as in your script)
taxon_cols <- function(lv, which = "set1_60") {
  lv <- as.character(lv)
  reg <- lv[!grepl("^Others?( |$|\\()", lv)]
  cols <- setNames(make_palette(length(reg), which), reg)
  oth <- setdiff(lv, reg)
  if (length(oth)) cols <- c(cols, setNames(rep("#B0B0B0", length(oth)), oth))
  cols[lv]
}

group_cols <- function(lv) setNames(make_palette(length(lv), "set1_60"), lv)

# ---------------------------------------------------------------------
# ENCODING REPAIR (ANSI_X3.4-1968 -> UTF-8 problems, e.g. "Temperature (degC)")
# ---------------------------------------------------------------------
enc_state <- new.env(); enc_state$n <- 0

fix_chr <- function(x) {
  if (is.factor(x)) { levels(x) <- fix_chr(levels(x)); return(x) }
  if (!is.character(x)) return(x)
  nonascii <- !is.na(x) & grepl("[^\\x01-\\x7F]", x, perl = TRUE, useBytes = TRUE)
  need <- nonascii & (Encoding(x) != "UTF-8")
  if (any(need)) {
    enc_state$n <- enc_state$n + sum(need)
    valid <- is.na(x) | validUTF8(x)
    v <- need & valid
    if (any(v)) Encoding(x[v]) <- "UTF-8"
    inv <- need & !valid
    if (any(inv)) x[inv] <- iconv(x[inv], "latin1", "UTF-8", sub = "byte")
  }
  x
}

fix_ps <- function(p) {
  enc_state$n <- 0
  sample_names(p) <- fix_chr(sample_names(p))
  taxa_names(p)   <- fix_chr(taxa_names(p))
  sdf <- sd_df(p)
  names(sdf) <- fix_chr(names(sdf))
  sdf[] <- lapply(sdf, fix_chr)
  rownames(sdf) <- fix_chr(rownames(sdf))
  sample_data(p) <- sample_data(sdf)
  tt <- as(tax_table(p), "matrix"); dn <- dimnames(tt)
  tt <- matrix(fix_chr(as.vector(tt)), nrow = nrow(tt),
               dimnames = list(fix_chr(dn[[1]]), fix_chr(dn[[2]])))
  tt[is.na(tt) | trimws(tt) == ""] <- "Unknown"      # as in your script
  tax_table(p) <- tax_table(tt)
  p
}

load_phyloseq <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj)) {
    env <- new.env()
    loaded <- tryCatch(load(path, envir = env), error = function(e) character(0))
    hits <- Filter(function(n) inherits(env[[n]], "phyloseq"), loaded)
    if (length(hits)) obj <- env[[hits[1]]]
  }
  if (is.null(obj) || !inherits(obj, "phyloseq")) stop("This file does not contain a phyloseq object.")
  if (is.null(tax_table(obj, errorIfNULL = FALSE))) stop("The phyloseq object has no taxonomy table.")
  if (is.null(sample_data(obj, errorIfNULL = FALSE)))
    sample_data(obj) <- sample_data(data.frame(Sample_ID = sample_names(obj), row.names = sample_names(obj)))
  obj
}

# ---------------------------------------------------------------------
# METADATA helpers (upload table + "group from sample-name pattern")
# ---------------------------------------------------------------------
read_meta <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) return(openxlsx::read.xlsx(path))
  sep <- if (ext == "csv") "," else "\t"
  utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
}

merge_meta <- function(p, mt) {
  sd0 <- sd_df(p)
  id <- fix_chr(as.character(mt[[1]]))
  mt <- mt[, -1, drop = FALSE]
  names(mt) <- fix_chr(names(mt)); mt[] <- lapply(mt, fix_chr)
  idx <- match(rownames(sd0), id)
  add <- mt[idx, , drop = FALSE]; rownames(add) <- rownames(sd0)
  keep0 <- setdiff(names(sd0), names(add))
  new <- cbind(sd0[, keep0, drop = FALSE], add)
  sample_data(p) <- sample_data(new)
  list(ps = p, n = sum(!is.na(idx)))
}

add_rule_var <- function(p, var, rules_text) {
  lines <- trimws(unlist(strsplit(rules_text, "\n")))
  lines <- lines[grepl("=", lines)]
  if (!length(lines) || !nzchar(trimws(var))) return(p)
  pat <- trimws(sub("=[^=]*$", "", lines)); lab <- trimws(sub("^.*=", "", lines))
  nm <- sample_names(p); out <- rep("Unknown", length(nm))
  for (i in rev(seq_along(pat))) out[grepl(pat[i], nm)] <- lab[i]   # first rule wins
  sdf <- sd_df(p); sdf[[trimws(var)]] <- out
  sample_data(p) <- sample_data(sdf)
  p
}

# ---------------------------------------------------------------------
# DATA helpers
# ---------------------------------------------------------------------
otu_mat <- function(p) {
  m <- as(otu_table(p), "matrix")
  if (!taxa_are_rows(p)) m <- t(m)
  m
}
tax_mat <- function(p) as(tax_table(p), "matrix")

agg_matrix <- function(p, rank) {            # taxa(rank) x samples, counts
  otu <- otu_mat(p)
  lab <- tax_mat(p)[, rank][match(rownames(otu), taxa_names(p))]
  rowsum(otu, group = lab, reorder = FALSE)
}
rel100 <- function(m) sweep(m, 2, pmax(colSums(m), 1e-12), "/") * 100

natural_key <- function(x) {
  vapply(as.character(x), function(s) {
    m <- gregexpr("[0-9]+", s)
    if (m[[1]][1] == -1) return(s)
    nums <- regmatches(s, m)[[1]]
    regmatches(s, m) <- list(formatC(as.numeric(nums), width = 12, flag = "0", format = "f", digits = 0))
    s
  }, character(1), USE.NAMES = FALSE)
}
nat_order <- function(...) order(...,  method = "radix")

# NOTE: as(sample_data, "data.frame") mangles names like "Temperature (degC)" -> keep them intact
sd_df <- function(p) {
  d <- data.frame(sample_data(p), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(d) <- sample_names(p)
  d
}
grp_of <- function(p, var) {
  nm <- sample_names(p)
  if (is.null(var) || var %in% c("None", "")) return(NULL)
  if (var == "Sample") return(setNames(nm, nm))
  g <- as.character(sd_df(p)[[var]]); g[is.na(g)] <- "NA"
  setNames(g, nm)
}

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
stars_blank <- function(p) sub("^ns$", "", stars(p))

is_num <- function(x) {
  if (is.numeric(x)) return(TRUE)
  y <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(y)) >= 0.8 && sum(!is.na(y)) >= 3 && length(unique(y[!is.na(y)])) > 2
}

# ---------------------------------------------------------------------
# CO-OCCURRENCE STATISTICS
# ---------------------------------------------------------------------
# Correlation coefficients and p-values for a taxa x samples count matrix.
#
# Both statistics come from the SAME matrix. Counts are converted to per-sample
# relative abundances first, because that is the scale the network is read on.
# Deriving r from relative abundances but p from raw counts puts two different
# analyses in one output row: when library sizes are uneven, normalising changes
# a taxon's rank order across samples, so the two statistics stop agreeing (e.g.
# a reported r of 0.7 carrying the p-value that belongs to r = 1).
cooc_stats <- function(cnt, method = "spearman") {
  rel <- rel100(cnt)
  rr  <- cor(t(rel), method = method)
  n   <- nrow(rr)
  pm  <- matrix(NA_real_, n, n, dimnames = dimnames(rr))
  if (n >= 2) for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    ct <- suppressWarnings(cor.test(rel[i, ], rel[j, ], method = method, exact = FALSE))
    pm[i, j] <- ct$p.value
  }
  # p.adjust drops NA before setting n, so the upper triangle is corrected alone
  padj <- matrix(p.adjust(as.vector(pm), "BH"), n, n, dimnames = dimnames(rr))
  list(r = rr, p = pm, p_adj = padj)
}

# ---------------------------------------------------------------------
# NETWORK NODE TAXONOMY
# ---------------------------------------------------------------------
# Co-occurrence node names are rank labels ("Alteromonas"), not tax_table
# rownames, so joining a node table to the taxonomy by rowname matches nothing
# at any rank above ASV. Worse, in a combined two-dataset network a dataset-2
# node whose label collides with a dataset-1 ASV id ("ASV1") silently picks up
# the other dataset's lineage. Resolve each label against its own taxonomy
# table at its own rank, taking the most frequent value per column.
node_tax_lookup <- function(tt, rank, labels) {
  labels <- unique(as.character(labels))
  if (is.null(tt) || !length(labels) || !ncol(tt)) return(NULL)
  by_asv <- is.null(rank) || identical(rank, "(ASV)") || !rank %in% colnames(tt)
  consensus <- function(lb) {
    rows <- if (by_asv) which(rownames(tt) == lb) else which(tt[, rank] == lb)
    if (!length(rows)) return(setNames(rep(NA_character_, ncol(tt)), colnames(tt)))
    apply(tt[rows, , drop = FALSE], 2, function(v) {
      v <- v[!is.na(v)]
      if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else NA_character_
    })
  }
  m <- vapply(labels, consensus, character(ncol(tt)))
  m <- matrix(m, nrow = ncol(tt), dimnames = list(colnames(tt), labels))
  data.frame(Taxon = labels, t(m), check.names = FALSE,
             stringsAsFactors = FALSE, row.names = NULL)
}

# row-bind two frames that may not share columns (missing cells become NA)
rbind_fill <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  cols <- union(names(a), names(b))
  for (cl in setdiff(cols, names(a))) a[[cl]] <- NA_character_
  for (cl in setdiff(cols, names(b))) b[[cl]] <- NA_character_
  rbind(a[cols], b[cols])
}

# exporters ----------------------------------------------------------
save_gg <- function(g, file, fmt, w, h, dpi) {
  if (fmt == "tiff") {
    grDevices::tiff(file, width = w * dpi, height = h * dpi, units = "px", res = dpi, compression = "lzw")
    on.exit(grDevices::dev.off())
    print(g)
  } else {
    dev <- if (fmt == "svg" && has("svglite")) svglite::svglite else fmt
    ggplot2::ggsave(file, g, device = dev, width = w, height = h, units = "in", dpi = dpi)
  }
}
save_base <- function(file, fmt, w, h, dpi, draw) {
  if (fmt == "svg") { if (has("svglite")) svglite::svglite(file, width = w, height = h) else grDevices::svg(file, width = w, height = h) }
  else if (fmt == "png") grDevices::png(file, width = w * dpi, height = h * dpi, res = dpi)
  else if (fmt == "tiff") grDevices::tiff(file, width = w * dpi, height = h * dpi, units = "px", res = dpi, compression = "lzw")
  else grDevices::pdf(file, width = w, height = h)
  on.exit(grDevices::dev.off())
  draw()
}
make_xlsx <- function(file, sheets) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    d <- sheets[[nm]]; if (is.null(d)) next
    sn <- substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, sn)
    openxlsx::writeData(wb, sn, d)
    openxlsx::setColWidths(wb, sn, cols = seq_len(ncol(d)), widths = "auto")
    openxlsx::freezePane(wb, sn, firstRow = TRUE)
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

dl_row <- function(fig, xl) div(downloadButton(fig, "Figure", class = "btn-sm"),
                                downloadButton(xl, "Excel tables", class = "btn-sm"))

# footer: which packages this tab uses (rendered dynamically in the server)
pkg_line <- function(pkgs) {
  pkgs <- as.character(pkgs)
  div(style = "margin:16px 0 2px 0; padding-top:5px; border-top:1px solid #eee; color:#999; font-size:11px;",
      HTML(paste0("Packages used in this tab: <code>", paste(pkgs, collapse = "</code> &middot; <code>"), "</code>")))
}
