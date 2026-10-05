# ---------------------------------------------------------------------
# Application server logic (internal)
# ---------------------------------------------------------------------

#' Application server (internal)
#' @noRd
app_server <- function(input, output, session) {
  # ---- encoding warnings: tell the user, then continue --------------
  seen <- new.env()
  notify_enc <- function(msg) {
    msg <- iconv(msg, "UTF-8", "UTF-8", sub = "byte"); key <- substr(msg, 1, 80)
    if (!isTRUE(seen[[key]])) {
      seen[[key]] <- TRUE
      showNotification(paste("Encoding problem ignored, continuing:", substr(msg, 1, 200)), type = "warning", duration = 12)
    }
  }
  enc_safe <- function(expr) {
    pat <- "cannot be translated|invalid in this locale|invalid multibyte|input string"
    withCallingHandlers(
      tryCatch(expr, error = function(e) {
        if (grepl(pat, conditionMessage(e))) { notify_enc(conditionMessage(e)); NULL } else stop(e)
      }),
      warning = function(w) {
        if (grepl(pat, conditionMessage(w))) { notify_enc(conditionMessage(w)); invokeRestart("muffleWarning") }
      })
  }
  show_gg <- function(g) enc_safe(withCallingHandlers(print(g),
    message = function(m) if (grepl("Too few points", conditionMessage(m))) invokeRestart("muffleMessage"),
    warning = function(w) if (grepl("^Removed [0-9]+ rows", conditionMessage(w))) invokeRestart("muffleWarning")))

  # ---- loading -------------------------------------------------------
  base <- reactiveVal(NULL)
  set_ps <- function(obj) {
    obj <- fix_ps(obj)
    if (enc_state$n > 0)
      showNotification(sprintf(paste0("Encoding issue found and fixed in %d text value(s) (e.g. a name like ",
                                      "'Temperature (\u00b0C)' flagged ANSI_X3.4-1968 but valid UTF-8). ",
                                      "Fixed automatically, nothing for you to do."), enc_state$n),
                       type = "warning", duration = 15)
    base(obj)
  }
  observeEvent(input$file, {
    res <- tryCatch(load_phyloseq(input$file$datapath), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL); NULL })
    if (!is.null(res)) set_ps(res)
  })
  observeEvent(input$demo, {
    e <- new.env(); data("GlobalPatterns", package = "phyloseq", envir = e); set_ps(e$GlobalPatterns)
  })

  rules <- reactiveVal(NULL)
  observeEvent(input$apply_rules, rules(list(name = input$rule_name, text = input$rule_text)))

  raw <- reactive({
    p <- base(); req(p)
    if (!is.null(input$meta_file)) {
      mt <- tryCatch(read_meta(input$meta_file$datapath, input$meta_file$name), error = function(e) {
        showNotification(paste("Could not read metadata:", conditionMessage(e)), type = "error"); NULL })
      if (!is.null(mt) && ncol(mt) >= 2) {
        r <- merge_meta(p, mt)
        if (r$n == 0) showNotification("No sample names in the metadata match the phyloseq object.", type = "error")
        p <- r$ps
      }
    }
    r <- rules()
    if (!is.null(r)) p <- add_rule_var(p, r$name, r$text)
    p
  })

  keep_sel <- function(id, choices, default) {
    cur <- isolate(input[[id]])
    if (is.null(cur) || !all(cur %in% choices)) cur <- default
    cur
  }
  observeEvent(raw(), {
    p <- raw(); ranks <- rank_names(p); sv <- sample_variables(p)
    numv <- sv[vapply(sv, function(v) is_num(sd_df(p)[[v]]), logical(1))]
    def <- ranks[grepl("phylum|division", ranks, ignore.case = TRUE)][1]; if (is.na(def)) def <- ranks[min(2, length(ranks))]
    deftm <- ranks[grepl("^class$", ranks, ignore.case = TRUE)][1]; if (is.na(deftm)) deftm <- def
    grp <- unique(c("Sample", sv)); grpn <- c("None", grp)
    updateCheckboxGroupInput(session, "samples", choices = sample_names(p),
                             selected = if (is.null(isolate(input$samples))) sample_names(p) else intersect(isolate(input$samples), sample_names(p)))
    updateSelectInput(session, "filt_rank", choices = c("None", ranks))
    for (id in c("rank", "b_taxrank")) updateSelectInput(session, id, choices = ranks, selected = keep_sel(id, ranks, def))
    updateSelectInput(session, "tm_rank", choices = ranks, selected = keep_sel("tm_rank", ranks, deftm))
    updateSelectInput(session, "hm_level", choices = c("ASV", ranks), selected = keep_sel("hm_level", c("ASV", ranks), "ASV"))
    updateSelectInput(session, "env_level", choices = c("ASV", ranks), selected = keep_sel("env_level", c("ASV", ranks), ranks[min(5, length(ranks))]))
    lr <- c("(none)", ranks)
    updateSelectInput(session, "hm_l1", choices = lr, selected = keep_sel("hm_l1", lr, lr[min(length(lr), 7)]))
    updateSelectInput(session, "hm_l2", choices = lr, selected = keep_sel("hm_l2", lr, lr[min(length(lr), 6)]))
    lab_choices <- c("ASV only" = "asv", "ASV_species" = "asv_sp", "ASV_species_class" = "asv_sp_cl",
                     "Custom (ranks below)" = "custom")
    updateSelectInput(session, "hm_lab", choices = lab_choices,
                      selected = keep_sel("hm_lab", unname(lab_choices), "asv_sp_cl"))
    for (id in c("facet", "tm_group", "hm_group", "b_group2", "env_group", "env_shape", "pat_group"))
      updateSelectInput(session, id, choices = grpn, selected = keep_sel(id, grpn, "None"))
    for (id in c("a_group", "b_group"))
      updateSelectInput(session, id, choices = grp, selected = keep_sel(id, grp, if (length(sv)) sv[1] else "Sample"))
    updateSelectizeInput(session, "env_vars", choices = numv, selected = keep_sel("env_vars", numv, head(numv, 4)))
  })

  observeEvent(input$filt_rank, {
    p <- raw(); req(p)
    ch <- if (input$filt_rank == "None") character(0) else sort(unique(tax_mat(p)[, input$filt_rank]))
    updateSelectizeInput(session, "filt_vals", choices = ch, selected = character(0))
  })

  # ---- shared filtered object ---------------------------------------
  ps_f <- reactive({
    p <- raw()
    validate(need(!is.null(p), "Upload a phyloseq file (or load the demo) to start."))
    req(input$samples)
    keep <- intersect(input$samples, sample_names(p))
    validate(need(length(keep) > 0, "Select at least one sample."))
    p <- prune_samples(keep, p)
    p <- prune_samples(sample_sums(p) >= input$min_depth, p)
    validate(need(nsamples(p) > 0, "No samples pass the read-depth filter."))
    p <- prune_taxa(taxa_sums(p) >= input$min_abund, p)
    p <- filter_taxa(p, function(x) sum(x > 0) >= input$min_prev, prune = TRUE)
    if (!identical(input$filt_rank, "None") && length(input$filt_vals)) {
      p <- prune_taxa(taxa_names(p)[tax_mat(p)[, input$filt_rank] %in% input$filt_vals], p)
    }
    validate(need(ntaxa(p) > 0, "No taxa pass the filters."))
    p
  })

  ps_rare <- reactive({
    p <- ps_f()
    if (!isTRUE(input$a_rare)) return(p)
    set.seed(42)
    suppressMessages(suppressWarnings(
      rarefy_even_depth(p, sample.size = min(sample_sums(p)), rngseed = 42, replace = FALSE,
                        trimOTUs = FALSE, verbose = FALSE)))
  })

  # ---- generic download builders ------------------------------------
  dl_gg <- function(id, getplot, base_name)
    output[[id]] <- downloadHandler(
      filename = function() paste0(base_name, ".", input$fmt),
      content = function(file) save_gg(getplot(), file, input$fmt, input$ex_w, input$ex_h, input$ex_dpi))
  dl_base <- function(id, drawer, base_name)
    output[[id]] <- downloadHandler(
      filename = function() paste0(base_name, ".", input$fmt),
      content = function(file) save_base(file, input$fmt, input$ex_w, input$ex_h, input$ex_dpi, drawer()))
  dl_xl <- function(id, getsheets, base_name)
    output[[id]] <- downloadHandler(filename = function() paste0(base_name, ".xlsx"),
                                    content = function(file) make_xlsx(file, getsheets()))
  dt <- function(df, ...) datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10), ...)

  # =================================================================
  # OVERVIEW
  # =================================================================
  output$ps_print <- renderPrint(enc_safe(ps_f()))

  seq_stats <- reactive({
    p <- ps_f(); lib <- sample_sums(p)
    rs <- refseq(p, errorIfNULL = FALSE)
    data.frame(Metric = c("Samples", "Total reads", "ASVs", if (!is.null(rs)) "Mean ASV length (bp)",
                          "Min library", "Median library", "Max library"),
               Value = c(nsamples(p), format(sum(lib), big.mark = ","), ntaxa(p),
                         if (!is.null(rs)) sprintf("%.1f", mean(nchar(as.character(rs)))),
                         format(min(lib), big.mark = ","), format(median(lib), big.mark = ","),
                         format(max(lib), big.mark = ",")), stringsAsFactors = FALSE)
  })
  rank_census <- reactive({
    tt <- tax_mat(ps_f())
    data.frame(Rank = colnames(tt),
               N_unique_taxa = vapply(colnames(tt), function(r) length(setdiff(unique(tt[, r]), "Unknown")), numeric(1)),
               row.names = NULL)
  })
  output$seq_stats   <- renderDT(enc_safe(datatable(seq_stats(), rownames = FALSE, options = list(dom = "t"))))
  output$rank_census <- renderDT(enc_safe(datatable(rank_census(), rownames = FALSE, options = list(dom = "t", pageLength = 20))))
  output$meta_tbl    <- renderDT(enc_safe(dt(sd_df(ps_f()))))

  output$depth_plot <- renderPlot(show_gg({
    s <- sort(sample_sums(ps_f()))
    df <- data.frame(Sample = factor(names(s), levels = names(s)), Reads = as.numeric(s))
    ggplot(df, aes(Sample, Reads)) + geom_col(fill = "#377EB8", width = 0.85) +
      labs(x = NULL, y = "Reads") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
  }))

  rare_gg <- reactive({
    p <- ps_f(); X <- t(otu_mat(p)); X <- X[, colSums(X) > 0, drop = FALSE]
    tot <- rowSums(X)
    out <- do.call(rbind, lapply(rownames(X), function(s) {
      d <- unique(round(seq(1, tot[s], length.out = min(100, max(2, tot[s])))))
      data.frame(Name = s, Depth = d, Richness = suppressWarnings(as.numeric(vegan::rarefy(round(X[s, ]), d))))
    }))
    out$Name <- factor(out$Name, levels = sort(unique(out$Name)))
    g <- ggplot(out, aes(Depth, Richness, color = Name)) + geom_line(linewidth = 0.5) +
      labs(x = "Sample Size", y = "Species Richness", color = "Description") + theme_grey(base_size = 11) +
      theme(panel.grid.minor = element_blank(), legend.text = element_text(size = 7),
            legend.title = element_text(size = 9, face = "bold"), legend.key.size = unit(0.4, "cm")) +
      guides(color = guide_legend(ncol = if (nlevels(out$Name) > 20) 2 else 1))
    if (isTRUE(input$rc_clip)) g <- g + coord_cartesian(xlim = c(0, min(tot)))
    g
  })
  output$rare_plot <- renderPlot(show_gg(rare_gg()))
  dl_gg("dl_rare", rare_gg, "rarefaction_curves")
  dl_xl("xl_over", function() list(Sequence_stats = seq_stats(), Rank_census = rank_census(),
                                   Metadata = data.frame(Sample = sample_names(ps_f()), sd_df(ps_f()), check.names = FALSE)),
        "overview_tables")

  # =================================================================
  # COMPOSITION (stacked bars)
  # =================================================================
  comp <- reactive({
    req(input$rank); p <- ps_f()
    cnt  <- agg_matrix(p, input$rank)          # counts: clade x samples
    aggm <- rel100(cnt)                        # % per sample (used for the bars)
    tot_reads <- sort(rowSums(cnt), decreasing = TRUE)
    tab <- data.frame(Clade = names(tot_reads), Reads = as.numeric(tot_reads),
                      stringsAsFactors = FALSE)
    tab$Pct    <- tab$Reads / sum(tab$Reads) * 100
    tab$CumPct <- cumsum(tab$Pct)
    n_keep <- if (input$comp_mode == "cum") {
      i <- which(tab$CumPct >= input$comp_cum)[1]; if (is.na(i)) nrow(tab) else i
    } else min(req(input$topn), nrow(tab))
    top <- tab$Clade[seq_len(n_keep)]
    M <- aggm[top, , drop = FALSE]
    if (nrow(aggm) > n_keep) M <- rbind(M, Others = colSums(aggm[!rownames(aggm) %in% top, , drop = FALSE]))
    nm <- colnames(M)
    g <- grp_of(p, input$facet)
    ord <- if (is.null(g)) seq_along(nm) else nat_order(g[nm], natural_key(nm))
    lv <- nm[ord]
    df <- data.frame(Sample = factor(rep(nm, each = nrow(M)), levels = lv),
                     Rank = factor(rep(rownames(M), ncol(M)), levels = rownames(M)),
                     Abundance = as.vector(M), stringsAsFactors = FALSE)
    if (!is.null(g)) df$Group <- factor(g[as.character(df$Sample)], levels = unique(g[lv]))
    tab$Shown <- ifelse(seq_len(nrow(tab)) <= n_keep, "shown", "in Others")
    tab$ASVs      <- as.integer(table(tax_mat(p)[, input$rank])[tab$Clade])
    tab$Pct_ASVs  <- round(tab$ASVs / sum(tab$ASVs) * 100, 2)
    list(df = df, tab = tab, n_keep = n_keep, shown = tab$CumPct[n_keep], aggm = aggm)
  })

  bar_gg <- reactive({
    cc <- comp(); df <- cc$df
    cols <- taxon_cols(levels(df$Rank), input$pal)
    g <- ggplot(df, aes(Sample, Abundance, fill = Rank)) + geom_col(width = 0.85) +
      scale_fill_manual(values = cols) +
      labs(y = "Relative abundance (%)", fill = input$rank,
           title = sprintf("%d clades = %.2f%% of total reads (rank: %s)", cc$n_keep, cc$shown, input$rank)) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = if (nlevels(df$Sample) > 40) 5 else 8),
            legend.title = element_text(face = "bold"), legend.text = element_text(size = 7),
            strip.text = element_text(face = "bold", size = 11), strip.background = element_rect(fill = "grey90"),
            plot.title = element_text(size = 11, face = "bold")) +
      guides(fill = guide_legend(reverse = FALSE))
    if ("Group" %in% names(df)) g <- g + facet_grid(~ Group, scales = "free_x", space = "free_x")
    g
  })
  output$bar_plot <- renderPlot(show_gg(bar_gg()))
  dl_gg("dl_bar", bar_gg, "composition_barplot")

  output$cum_tbl <- renderDT(enc_safe({
    t <- comp()$tab
    t$Pct <- round(t$Pct, 2); t$CumPct <- round(t$CumPct, 2)
    names(t)[names(t) == "ASVs"]     <- "n_ASVs"
    names(t)[names(t) == "Pct_ASVs"] <- "ASV_%"
    names(t)[names(t) == "Reads"]    <- "n_reads"
    names(t)[names(t) == "Pct"]      <- "Reads_%"
    names(t)[names(t) == "CumPct"]   <- "Cum_reads_%"
    dt(t[, c("Clade", "n_ASVs", "ASV_%", "n_reads", "Reads_%", "Cum_reads_%", "Shown")])
  }))

  grp_summary <- reactive({
    cc <- comp(); p <- ps_f(); m <- cc$aggm            # taxa x samples (% per sample)
    tab <- data.frame(Taxon = rownames(m), ASVs = as.integer(table(tax_mat(p)[, input$rank])[rownames(m)]),
                      Overall_pct = round(rowMeans(m), 2), stringsAsFactors = FALSE)
    g <- grp_of(p, input$facet)
    if (!is.null(g)) for (lv in unique(g[colnames(m)]))
      tab[[paste0(gsub("[^A-Za-z0-9]+", "_", lv), "_pct")]] <- round(rowMeans(m[, g[colnames(m)] == lv, drop = FALSE]), 2)
    tab[order(-tab$Overall_pct), ]
  })
  output$grp_tbl <- renderDT(enc_safe(dt(grp_summary())))
  dl_xl("xl_bar", function() {
    cc <- comp()
    list(Clade_percentages = cc$tab, Taxon_summary = grp_summary(),
         Per_sample_pct = data.frame(Taxon = rownames(cc$aggm), round(cc$aggm, 4), check.names = FALSE, row.names = NULL))
  }, "composition_tables")

  # =================================================================
  # TREEMAP
  # =================================================================
  tm_data <- reactive({
    req(input$tm_rank); p <- ps_f()
    m <- rel100(agg_matrix(p, input$tm_rank)); g <- grp_of(p, input$tm_group)
    groups <- if (is.null(g)) list(All = colnames(m)) else split(colnames(m), g[colnames(m)])
    thr <- input$tm_thr
    dfs <- lapply(groups, function(s) {
      pct <- rowMeans(m[, s, drop = FALSE]); d <- data.frame(Class = names(pct), Pct = as.numeric(pct), stringsAsFactors = FALSE)
      keep <- d$Pct > thr; oth <- sum(d$Pct[!keep]); d <- d[keep, , drop = FALSE]
      olab <- sprintf("Other (<%s%%)", format(thr))
      if (oth > 0.01) d <- rbind(d, data.frame(Class = olab, Pct = oth))
      d$Label <- paste0(d$Class, "\n", sprintf("%.2f%%", d$Pct))
      d[order(-d$Pct), ]
    })
    allc <- unique(unlist(lapply(dfs, `[[`, "Class")))
    cols <- taxon_cols(sub("^Other \\(.*", "Others", allc), input$pal); names(cols) <- allc
    cols[grepl("^Other \\(", allc)] <- "#BFBFBF"
    dfs <- lapply(dfs, function(d) { d$Color <- unname(cols[d$Class]); d })
    list(dfs = dfs, thr = thr)
  })
  draw_tm <- function() {
    dfs <- tm_data()$dfs; n <- length(dfs); nc <- min(3, n); nr <- ceiling(n / nc)
    grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, nc)))
    for (i in seq_len(n)) {
      grid::pushViewport(grid::viewport(layout.pos.row = ceiling(i / nc), layout.pos.col = (i - 1) %% nc + 1))
      treemap::treemap(dfs[[i]], index = "Label", vSize = "Pct", vColor = "Color", type = "color",
                       title = names(dfs)[i], vp = grid::viewport(width = 0.98, height = 0.98),
                       fontsize.labels = 14, fontface.labels = "bold", inflate.labels = FALSE,
                       border.col = "black", border.lwds = 1.5)
      grid::popViewport()
    }
    grid::popViewport()
  }
  if (has("treemap")) output$tm_plot <- renderPlot(enc_safe(draw_tm()))
  dl_base("dl_tm", function() draw_tm, "treemap")
  dl_xl("xl_tm", function() {
    dfs <- tm_data()$dfs
    long <- do.call(rbind, lapply(names(dfs), function(n) data.frame(Group = n, Class = dfs[[n]]$Class, Pct = round(dfs[[n]]$Pct, 2))))
    wide <- Reduce(function(x, y) merge(x, y, by = "Class", all = TRUE),
                   lapply(names(dfs), function(n) { d <- dfs[[n]][, c("Class", "Pct")]; names(d)[2] <- gsub("[^A-Za-z0-9]+", "_", n); d }))
    wide$Total_mean_pct <- round(rowMeans(wide[, -1, drop = FALSE], na.rm = TRUE), 2)
    wide[, -1] <- lapply(wide[, -1, drop = FALSE], round, 2)
    list(Class_by_group = wide[order(-wide$Total_mean_pct), ], Treemap_values = long)
  }, "treemap_values")

  # =================================================================
  # HEATMAP
  # =================================================================
  hm_calc <- reactive({
    req(input$hm_level, input$hm_n); p <- ps_f(); tt <- tax_mat(p)
    if (input$hm_level == "ASV") {
      cnt <- otu_mat(p); rel <- rel100(cnt)
      top <- names(sort(rowSums(cnt), decreasing = TRUE))[seq_len(min(input$hm_n, nrow(cnt)))]
      mat <- rel[top, , drop = FALSE]
      hm_lab <- if (is.null(input$hm_lab)) "asv_sp_cl" else input$hm_lab
      lab_ranks <- switch(hm_lab,
                          asv       = character(0),
                          asv_sp    = intersect("Species", colnames(tt)),
                          asv_sp_cl = intersect(c("Species", "Class"), colnames(tt)),
                          custom    = intersect(c(input$hm_l1, input$hm_l2), colnames(tt)))
      lab_ranks <- setdiff(lab_ranks, "(none)")
      lab <- vapply(top, function(a) {
        parts <- c(a)
        for (r in lab_ranks) {
          v <- tt[a, r]; if (is.na(v) || trimws(v) == "") v <- "Unknown"
          parts <- c(parts, gsub(" ", "_", v))
        }
        paste(parts, collapse = "_")
      }, character(1))
      rownames(mat) <- make.unique(unname(lab)); tax_info <- data.frame(ID = top, tt[top, , drop = FALSE], check.names = FALSE)
      total <- rowSums(cnt)[top]
    } else {
      cnt <- agg_matrix(p, input$hm_level); rel <- rel100(cnt)
      top <- names(sort(rowSums(cnt), decreasing = TRUE))[seq_len(min(input$hm_n, nrow(cnt)))]
      mat <- rel[top, , drop = FALSE]; tax_info <- data.frame(Taxon = top); total <- rowSums(cnt)[top]
    }
    z <- t(scale(t(mat))); z[is.na(z)] <- 0
    cap <- input$hm_cap; z[z > cap] <- cap; z[z < -cap] <- -cap
    g <- grp_of(p, input$hm_group); nm <- colnames(z)
    if (input$hm_cols == "blocks") {
      ord <- if (is.null(g)) order(natural_key(nm)) else nat_order(g[nm], natural_key(nm))
      z <- z[, ord, drop = FALSE]; mat <- mat[, ord, drop = FALSE]
    }
    ann <- NULL; ann_col <- NULL; gaps <- NULL
    if (!is.null(g)) {
      ann <- data.frame(g[colnames(z)], row.names = colnames(z)); names(ann) <- input$hm_group
      lv <- unique(ann[[1]]); ann_col <- setNames(list(group_cols(lv)), input$hm_group)
      if (input$hm_cols == "blocks") gaps <- head(cumsum(table(factor(ann[[1]], levels = lv))), -1)
    }
    list(z = z, mat = mat, ann = ann, ann_col = ann_col, gaps = gaps, tax_info = tax_info, total = total)
  })
  hm_obj <- reactive({
    h <- hm_calc(); clust_cols <- input$hm_cols == "clust" && ncol(h$z) > 1
    pheatmap::pheatmap(h$z, cluster_rows = nrow(h$z) > 1, clustering_method = "ward.D2",
                       cluster_cols = clust_cols, clustering_method_cols = "ward.D2", gaps_col = h$gaps,
                       color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                       breaks = seq(-input$hm_cap, input$hm_cap, length.out = 101),
                       main = sprintf("Top %d (%s) - row Z-scores", nrow(h$z), input$hm_level),
                       fontsize_row = 8, fontsize_col = 9, angle_col = 90, border_color = NA,
                       annotation_col = h$ann, annotation_colors = h$ann_col,
                       annotation_names_col = FALSE, silent = TRUE)
  })
  draw_hm <- function() { grid::grid.newpage(); grid::grid.draw(hm_obj()$gtable) }
  if (has("pheatmap")) output$hm_plot <- renderPlot(enc_safe(draw_hm()))
  dl_base("dl_hm", function() draw_hm, "heatmap")
  dl_xl("xl_hm", function() {
    h <- hm_calc(); o <- hm_obj()
    ro <- if (!is.null(o$tree_row)) o$tree_row$order else seq_len(nrow(h$z))
    co <- if (!is.null(o$tree_col)) o$tree_col$order else seq_len(ncol(h$z))
    z <- h$z[ro, co, drop = FALSE]; m <- h$mat[ro, co, drop = FALSE]
    list(Heatmap_order = data.frame(Row_in_figure = seq_len(nrow(z)), Label = rownames(z), h$tax_info[ro, , drop = FALSE],
                                    Total_reads = as.numeric(h$total[ro]), Mean_pct = round(rowMeans(m), 3), check.names = FALSE, row.names = NULL),
         Zscores = data.frame(Label = rownames(z), round(z, 3), check.names = FALSE, row.names = NULL),
         Relative_pct = data.frame(Label = rownames(m), round(m, 4), check.names = FALSE, row.names = NULL))
  }, "heatmap_tables")

  # =================================================================
  # ALPHA DIVERSITY
  # =================================================================
  observe({
    p <- raw(); req(p)
    tr <- !is.null(phy_tree(p, errorIfNULL = FALSE)) && has("picante")
    ch <- c("Observed", "Simpson", "Chao1", "Shannon", if (tr) "PD")
    updateCheckboxGroupInput(session, "a_meas", choices = ch, selected = intersect(isolate(input$a_meas), ch))
  })

  alpha_df <- reactive({
    p <- ps_rare(); X <- t(otu_mat(p)); X <- round(X)
    out <- data.frame(Sample = rownames(X), Observed = rowSums(X > 0),
                      Simpson = vegan::diversity(X, "simpson"), Shannon = vegan::diversity(X, "shannon"),
                      Chao1 = suppressWarnings(as.numeric(t(vegan::estimateR(X))[, "S.chao1"])), stringsAsFactors = FALSE)
    tr <- phy_tree(p, errorIfNULL = FALSE)
    if (!is.null(tr) && has("picante")) {
      out$PD <- tryCatch(picante::pd(X, tr, include.root = ape::is.rooted(tr))$PD, error = function(e) NA_real_)
    }
    out
  })
  alpha_long <- reactive({
    req(input$a_meas, input$a_group); a <- alpha_df(); p <- ps_f(); g <- grp_of(p, input$a_group)
    idx <- intersect(c("Observed", "Simpson", "Chao1", "Shannon", "PD"), input$a_meas); idx <- intersect(idx, names(a))
    d <- do.call(rbind, lapply(idx, function(m) data.frame(Sample = a$Sample, Group = g[a$Sample], Index = m, value = a[[m]], stringsAsFactors = FALSE)))
    d$Index <- factor(d$Index, levels = idx)
    d[is.finite(d$value), ]
  })
  alpha_stats <- reactive({
    d <- alpha_long(); validate(need(input$a_group != "Sample" && length(unique(d$Group)) >= 2,
                                     "Pick a metadata variable with at least 2 groups to run the statistics."))
    d$Group <- factor(d$Group)
    norm <- list(); varr <- list(); glob <- list(); post <- list(); meth <- c()
    for (idx in levels(d$Index)) {
      di <- d[d$Index == idx, ]; if (nrow(di) < 3) next
      sh <- do.call(rbind, lapply(split(di, di$Group), function(g) {
        n <- nrow(g); ok <- n >= 3 && n <= 5000 && length(unique(g$value)) > 1
        sw <- if (ok) shapiro.test(g$value) else NULL
        data.frame(Index = idx, Group = as.character(g$Group[1]), n = n, W = if (ok) unname(sw$statistic) else NA_real_,
                   p_value = if (ok) sw$p.value else NA_real_, stringsAsFactors = FALSE)
      }))
      sh$Interpretation <- ifelse(sh$n < 3, "Skipped (n < 3)", ifelse(is.na(sh$W), "Skipped (constant or invalid)",
                                  ifelse(sh$p_value >= 0.05, "Normal (fail to reject H0)", "NOT normal (reject H0)")))
      norm[[idx]] <- sh
      normal_ok <- all(!is.na(sh$p_value) & sh$p_value >= 0.05)
      bt <- tryCatch(bartlett.test(value ~ Group, data = di), error = function(e) NULL)
      var_ok <- !is.null(bt) && bt$p.value >= 0.05
      if (!is.null(bt)) varr[[idx]] <- data.frame(Index = idx, Statistic = unname(bt$statistic), df = unname(bt$parameter), p_value = bt$p.value,
                                                  Interpretation = ifelse(var_ok, "Equal variances (fail to reject H0)", "UNEQUAL variances (reject H0)"))
      if (normal_ok && var_ok) {
        fit <- aov(value ~ Group, data = di); r <- summary(fit)[[1]]
        glob[[idx]] <- data.frame(Index = idx, Test = "ANOVA (F-test)", Statistic = round(r[1, "F value"], 4), df1 = r[1, "Df"], df2 = r[2, "Df"],
                                  p_value = signif(r[1, "Pr(>F)"], 4), Interpretation = ifelse(r[1, "Pr(>F)"] < 0.05, "Significant", "Not significant"))
        tk <- as.data.frame(TukeyHSD(fit)$Group)
        post[[idx]] <- data.frame(Index = idx, Comparison = rownames(tk), Difference = tk$diff, Lower_CI = tk$lwr, Upper_CI = tk$upr,
                                  p_adj = tk$`p adj`, Method = "Tukey HSD", Significance = stars(tk$`p adj`))
        meth[idx] <- "p"
      } else {
        kw <- kruskal.test(value ~ Group, data = di)
        glob[[idx]] <- data.frame(Index = idx, Test = "Kruskal-Wallis", Statistic = round(unname(kw$statistic), 4), df1 = unname(kw$parameter), df2 = NA,
                                  p_value = signif(kw$p.value, 4), Interpretation = ifelse(kw$p.value < 0.05, "Significant", "Not significant"))
        pw <- tryCatch(pairwise.wilcox.test(di$value, di$Group, p.adjust.method = "BH", exact = FALSE)$p.value, error = function(e) NULL)
        if (!is.null(pw)) {
          pl <- as.data.frame(as.table(pw)); pl <- pl[!is.na(pl$Freq), ]
          if (nrow(pl)) post[[idx]] <- data.frame(Index = idx, Comparison = paste(pl$Var2, "vs", pl$Var1), Difference = NA, Lower_CI = NA, Upper_CI = NA,
                                                   p_adj = pl$Freq, Method = "Pairwise Wilcoxon (BH)", Significance = stars(pl$Freq))
        }
        meth[idx] <- "np"
      }
    }
    bind <- function(l) if (length(l)) do.call(rbind, c(l, make.row.names = FALSE)) else NULL
    list(norm = bind(norm), var = bind(varr), glob = bind(glob), post = bind(post), meth = meth)
  })
  output$a_note <- renderUI({
    if (!has("ggpubr")) div(class = "alert alert-info", "Install 'ggpubr' to draw p-value brackets on the boxplots: install.packages('ggpubr').")
  })

  a_dot_gg <- reactive({
    d <- alpha_long(); g <- d$Group; ord <- unique(d$Sample[order(g, natural_key(d$Sample))])
    d$Sample <- factor(d$Sample, levels = rev(ord))
    cols <- group_cols(sort(unique(d$Group)))
    ggplot(d, aes(value, Sample, color = Group)) + geom_point(size = 2.3) + scale_color_manual(values = cols) +
      facet_wrap(~ Index, scales = "free_x", nrow = 1) + labs(x = NULL, y = NULL, color = input$a_group) +
      theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 6), strip.background = element_rect(fill = "grey90"),
                                       strip.text = element_text(face = "bold"))
  })
  a_box_gg <- reactive({
    d <- alpha_long(); st <- alpha_stats(); d$Group <- factor(d$Group); lv <- levels(d$Group)
    cols <- group_cols(lv)
    g <- ggplot(d, aes(Group, value, fill = Group)) + geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.4) +
      geom_jitter(width = 0.12, size = 1.6, color = "black", alpha = 0.7) + scale_fill_manual(values = cols) +
      facet_wrap(~ Index, scales = "free_y", nrow = 1) + labs(x = NULL, y = NULL, fill = input$a_group) +
      theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "grey90"),
                                       strip.text = element_text(face = "bold"), legend.position = "none")
    if (has("ggpubr")) {
      use_np <- switch(input$a_test, p = FALSE, np = TRUE, auto = any(st$meth == "np"))
      gm <- if (use_np) "kruskal.test" else "anova"; pm <- if (use_np) "wilcox.test" else "t.test"
      g <- g + ggpubr::stat_compare_means(method = gm, label.y.npc = "top", size = 3, vjust = -0.3)
      if (length(lv) >= 2 && length(lv) <= 5)
        g <- g + ggpubr::stat_compare_means(comparisons = combn(lv, 2, simplify = FALSE), method = pm, label = "p.signif", hide.ns = TRUE)
      g <- g + scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
    }
    g
  })
  output$a_dot <- renderPlot(show_gg(a_dot_gg()))
  output$a_box <- renderPlot(show_gg(a_box_gg()))
  dl_gg("dl_alpha", a_box_gg, "alpha_boxplots")
  output$a_vals <- renderDT(enc_safe({ a <- alpha_df(); a[-1] <- lapply(a[-1], round, 3)
    dt(data.frame(a, Group = if (input$a_group == "Sample") a$Sample else grp_of(ps_f(), input$a_group)[a$Sample])) }))
  for (nm in c("norm", "var", "glob", "post")) local({ n <- nm
    output[[paste0("a_", n)]] <- renderDT(enc_safe({ x <- alpha_stats()[[n]]; validate(need(!is.null(x), "No result.")); dt(x) })) })
  dl_xl("xl_alpha", function() { st <- alpha_stats()
    list(Alpha_values = alpha_df(), Normality = st$norm, Variance_Bartlett = st$var, Global_test = st$glob, Posthoc = st$post) }, "alpha_diversity_stats")

  # =================================================================
  # BETA DIVERSITY
  # =================================================================
  beta <- reactive({
    p <- ps_f(); validate(need(nsamples(p) >= 3, "Ordination needs at least 3 samples."))
    otu <- otu_mat(p); tot <- sort(rowSums(otu), decreasing = TRUE)
    n_ret <- which(cumsum(tot) / sum(tot) * 100 >= input$b_cum)[1]; if (is.na(n_ret)) n_ret <- length(tot)
    keep <- names(tot)[seq_len(n_ret)]; Y <- t(otu[keep, , drop = FALSE])
    Y <- switch(input$b_trans, sqrt = sqrt(Y), rel = Y / pmax(rowSums(Y), 1e-12), hell = vegan::decostand(Y, "hellinger"), Y)
    set.seed(42)
    d <- vegan::vegdist(Y, method = input$b_dist, binary = input$b_dist == "jaccard")
    hc <- hclust(d, "average"); pc <- cmdscale(d, k = 2, eig = TRUE)
    ev <- pc$eig; var_exp <- 100 * ev[ev > 0] / sum(ev[ev > 0])
    sp <- as.data.frame(vegan::wascores(pc$points, Y)); names(sp) <- c("Axis1", "Axis2"); sp$ASV <- rownames(sp)
    list(d = d, hc = hc, pc = pc, var = var_exp, Y = Y, sp = sp, n_ret = n_ret, n_all = length(tot),
         pct_ret = cumsum(tot)[n_ret] / sum(tot) * 100, p = p)
  })
  b_meta <- reactive({
    b <- beta(); p <- b$p; nm <- rownames(b$pc$points)
    data.frame(Sample = nm, F1 = grp_of(p, input$b_group)[nm],
               F2 = if (input$b_group2 == "None") NA_character_ else grp_of(p, input$b_group2)[nm], stringsAsFactors = FALSE)
  })
  b_cols <- reactive(group_cols(sort(unique(b_meta()$F1))))

  b_dend_gg <- reactive({
    b <- beta(); m <- b_meta(); validate(need(has("ggdendro"), "Install 'ggdendro' for the dendrogram: install.packages('ggdendro')"))
    dd <- ggdendro::dendro_data(b$hc, type = "rectangle"); lab <- dd$labels
    lab$F1 <- m$F1[match(lab$label, m$Sample)]; ym <- max(ggdendro::segment(dd)$y)
    ggplot() + geom_segment(data = ggdendro::segment(dd), aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.3) +
      geom_tile(data = lab, aes(x = x, y = -0.04 * ym, fill = F1), height = 0.06 * ym, width = 0.9, alpha = 0.9) +
      geom_text(data = lab, aes(x = x, y = 0, label = label, color = F1), angle = 90, hjust = 1.05, vjust = 0.5, size = 2.4, show.legend = FALSE) +
      scale_fill_manual(values = b_cols(), name = input$b_group) + scale_color_manual(values = b_cols()) +
      coord_cartesian(ylim = c(-0.12 * ym, ym * 1.05), clip = "off") +
      labs(x = NULL, y = paste(input$b_dist, "dissimilarity")) + theme_classic() +
      theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(b = 60))
  })
  output$b_dend <- renderPlot(enc_safe({
    if (has("ggdendro")) print(b_dend_gg()) else plot(beta()$hc, main = "Hierarchical clustering", xlab = "", sub = "")
  }))

  b_pcoa_gg <- reactive({
    b <- beta(); m <- b_meta(); df <- data.frame(PCo1 = b$pc$points[, 1], PCo2 = b$pc$points[, 2], m)
    shp <- !is.na(df$F2[1])
    g <- ggplot(df, aes(PCo1, PCo2, color = F1, shape = if (shp) F2 else NULL)) + geom_point(size = 3, alpha = 0.85, stroke = 1)
    if (isTRUE(input$b_ell)) g <- g + stat_ellipse(aes(x = PCo1, y = PCo2, color = F1), inherit.aes = FALSE, type = "norm", linetype = 2, linewidth = 0.8, show.legend = FALSE)
    g <- g + scale_color_manual(values = b_cols()) +
      labs(x = sprintf("PCo1 (%.1f%%)", b$var[1]), y = sprintf("PCo2 (%.1f%%)", b$var[2]), color = input$b_group,
           shape = if (shp) input$b_group2 else NULL) +
      theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), aspect.ratio = 1)
    if (shp) g <- g + scale_shape_manual(values = rep(c(16, 17, 15, 18, 3, 4, 8, 7, 9, 10), 5)[seq_along(unique(df$F2))])
    g
  })
  output$b_pcoa <- renderPlot(show_gg(b_pcoa_gg()))

  b_bi_gg <- reactive({
    b <- beta(); m <- b_meta(); req(input$b_taxrank)
    ss <- data.frame(Axis1 = b$pc$points[, 1], Axis2 = b$pc$points[, 2], m)
    sp <- b$sp; tt <- tax_mat(b$p); sp$Taxon <- tt[sp$ASV, input$b_taxrank]
    top <- names(sort(table(sp$Taxon), decreasing = TRUE))[seq_len(min(12, length(unique(sp$Taxon))))]
    sp$Tax_plot <- ifelse(sp$Taxon %in% top, sp$Taxon, "Others")
    sp$Tax_plot <- factor(sp$Tax_plot, levels = c(setdiff(top, "Others"), "Others")[c(setdiff(top, "Others"), "Others") %in% sp$Tax_plot])
    xl <- range(c(ss$Axis1, sp$Axis1)) * 1.05; yl <- range(c(ss$Axis2, sp$Axis2)) * 1.05
    pa <- ggplot(ss, aes(Axis1, Axis2)) + geom_text(aes(label = Sample, color = F1), size = 2.5, fontface = "bold", show.legend = FALSE) +
      coord_cartesian(xlim = xl, ylim = yl) + scale_color_manual(values = b_cols()) +
      labs(title = "Samples", x = sprintf("Axis 1 [%.1f%%]", b$var[1]), y = sprintf("Axis 2 [%.1f%%]", b$var[2])) +
      theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
    pb <- ggplot(sp, aes(Axis1, Axis2, color = Tax_plot)) + geom_point(size = 1.8, alpha = 0.7) +
      coord_cartesian(xlim = xl, ylim = yl) + scale_color_manual(values = taxon_cols(levels(sp$Tax_plot), input$pal)) +
      labs(title = "Taxa", x = sprintf("Axis 1 [%.1f%%]", b$var[1]), y = sprintf("Axis 2 [%.1f%%]", b$var[2]), color = input$b_taxrank) +
      theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank(),
                                       legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))
    if (has("patchwork")) patchwork::wrap_plots(pa, pb, nrow = 1)
    else if (has("ggpubr")) ggpubr::ggarrange(pa, pb, ncol = 2) else pb
  })
  output$b_bi <- renderPlot(show_gg(b_bi_gg()))
  dl_gg("dl_beta", b_bi_gg, "beta_biplot")

  b_stats <- reactive({
    b <- beta(); m <- b_meta(); n <- input$b_perm; set.seed(42)
    d <- b$d; facs <- list(F1 = input$b_group); if (!is.na(m$F2[1])) facs$F2 <- input$b_group2
    df <- data.frame(F1 = factor(m$F1), F2 = factor(m$F2))
    perm <- NULL; disp <- NULL
    for (f in names(facs)) {
      r <- tryCatch(vegan::adonis2(d ~ g, data = data.frame(g = df[[f]]), permutations = n, by = "margin"), error = function(e) NULL)
      if (!is.null(r)) perm <- rbind(perm, data.frame(Factor = facs[[f]], Df = r$Df[1], F = r$F[1], R2 = r$R2[1], p = r$`Pr(>F)`[1]))
      bd <- tryCatch(vegan::permutest(vegan::betadisper(d, df[[f]]), permutations = n), error = function(e) NULL)
      if (!is.null(bd)) disp <- rbind(disp, data.frame(Factor = facs[[f]], F = bd$tab$F[1], p = bd$tab$`Pr(>F)`[1]))
    }
    if (length(facs) == 2) {
      r <- tryCatch(vegan::adonis2(d ~ F1 * F2, data = df, permutations = n, by = "terms"), error = function(e) NULL)
      if (!is.null(r) && "F1:F2" %in% rownames(r)) perm <- rbind(perm, data.frame(Factor = paste(facs$F1, "x", facs$F2),
                                                  Df = r["F1:F2", "Df"], F = r["F1:F2", "F"], R2 = r["F1:F2", "R2"], p = r["F1:F2", "Pr(>F)"]))
    }
    if (!is.null(perm)) perm$Significance <- stars(perm$p)
    if (!is.null(disp)) disp$Significance <- ifelse(disp$p < 0.05, "Dispersion differs", "Homogeneous dispersion")
    pair <- NULL; lv <- levels(df$F1); dm <- as.matrix(d)
    if (length(lv) >= 2) for (i in seq_len(length(lv) - 1)) for (j in (i + 1):length(lv)) {
      keep <- df$F1 %in% lv[c(i, j)]; if (sum(keep) < 3) next
      r <- tryCatch(vegan::adonis2(as.dist(dm[keep, keep]) ~ g, data = data.frame(g = droplevels(df$F1[keep])), permutations = n), error = function(e) NULL)
      if (!is.null(r)) pair <- rbind(pair, data.frame(Comparison = paste(lv[i], "vs", lv[j]), R2 = r$R2[1], F = r$F[1], p = r$`Pr(>F)`[1]))
    }
    if (!is.null(pair)) { pair$p_adj_BH <- p.adjust(pair$p, "BH"); pair$Significance <- stars(pair$p) }
    meth <- data.frame(Parameter = c("Total ASVs", "ASVs retained (cumulative cutoff)", "Percentage of reads retained", "Transformation",
                                     "Distance", "Ordination", "PCoA axis 1 variance", "PCoA axis 2 variance", "Permutations"),
                       Value = c(b$n_all, b$n_ret, sprintf("%.2f%%", b$pct_ret), input$b_trans, input$b_dist, "PCoA",
                                 sprintf("%.2f%%", b$var[1]), sprintf("%.2f%%", b$var[2]), n))
    list(perm = perm, disp = disp, pair = pair, meth = meth)
  })
  for (nm in c("perm", "disp", "pair", "meth")) local({ n <- nm
    output[[paste0("b_", n, "_t")]] <- renderDT(enc_safe({ x <- b_stats()[[n]]
      validate(need(!is.null(x), "Not available (need >= 2 groups and enough samples).")); x[] <- lapply(x, function(v) if (is.numeric(v)) signif(v, 4) else v); dt(x) })) })
  dl_xl("xl_beta", function() { s <- b_stats(); b <- beta(); m <- b_meta()
    list(Method_Summary = s$meth, PERMANOVA = s$perm, Betadisper = s$disp, Pairwise_PERMANOVA = s$pair,
         Sample_scores = data.frame(m, PCo1 = b$pc$points[, 1], PCo2 = b$pc$points[, 2], row.names = NULL),
         Taxa_scores = data.frame(b$sp, tax_mat(b$p)[b$sp$ASV, , drop = FALSE], check.names = FALSE, row.names = NULL)) }, "beta_diversity_stats")

  # =================================================================
  # ENVIRONMENT: RDA + correlations
  # =================================================================
  two_cols <- function(m) { m <- as.matrix(m); if (ncol(m) < 2) m <- cbind(m, 0); m[, 1:2, drop = FALSE] }

  rda_run <- reactive({
    req(input$env_vars, input$env_level); p <- ps_f(); sd <- sd_df(p)
    env <- sd[, input$env_vars, drop = FALSE]; env[] <- lapply(env, function(x) suppressWarnings(as.numeric(as.character(x))))
    lab_map <- setNames(names(env), make.names(names(env), unique = TRUE)); names(env) <- names(lab_map)
    ok <- stats::complete.cases(env); validate(need(sum(ok) >= 4, "Need >= 4 samples with complete values for the selected variables."))
    p <- prune_samples(rownames(env)[ok], p); env <- env[ok, , drop = FALSE]
    env <- env[, vapply(env, function(x) sd(x) > 0, logical(1)), drop = FALSE]; validate(need(ncol(env) >= 1, "Selected variables are constant."))
    m <- if (input$env_level == "ASV") otu_mat(p) else agg_matrix(p, input$env_level)
    m <- m[rowSums(m) > 0, , drop = FALSE]; m <- m[order(rowSums(m), decreasing = TRUE), , drop = FALSE]
    Yfull <- vegan::decostand(t(m)[rownames(env), , drop = FALSE], "hellinger")
    top <- rownames(m)[seq_len(min(input$env_topn, nrow(m)))]; Y <- Yfull[, top, drop = FALSE]
    mean_pct <- colMeans(rel100(m)[top, rownames(env), drop = FALSE] |> t())
    envs <- as.data.frame(scale(env))
    cand <- vapply(names(envs), function(v) { x <- envs[[v]]; tryCatch(vegan::RsquareAdj(vegan::rda(Y ~ x))$adj.r.squared, error = function(e) NA_real_) }, numeric(1))
    cand_df <- data.frame(Variable = unname(lab_map[names(cand)]), Marginal_adj_R2 = round(cand, 4), row.names = NULL)
    dropped <- data.frame(Variable = character(), Reason = character(), stringsAsFactors = FALSE)
    repeat {
      if (ncol(envs) < 2) break; cm <- abs(cor(envs)); diag(cm) <- 0; if (max(cm) < input$env_r) break
      ij <- which(cm == max(cm), arr.ind = TRUE)[1, ]; pr <- colnames(cm)[ij]; dv <- pr[which.max(colMeans(cm)[pr])]
      dropped <- rbind(dropped, data.frame(Variable = lab_map[dv], Reason = sprintf("|r| >= %.2f with %s", input$env_r, lab_map[setdiff(pr, dv)[1]])))
      envs[[dv]] <- NULL
    }
    while (ncol(envs) > nrow(envs) - 2 && ncol(envs) > 1) {
      dv <- names(which.min(cand[names(envs)])); dropped <- rbind(dropped, data.frame(Variable = lab_map[dv], Reason = "more variables than degrees of freedom allow")); envs[[dv]] <- NULL
    }
    vif_removed <- character()
    repeat {
      mod <- vegan::rda(Y ~ ., data = envs)
      v <- tryCatch(vegan::vif.cca(mod), error = function(e) NULL)
      if (!isTRUE(input$env_vif) || is.null(v) || ncol(envs) <= 1) break
      v[is.na(v)] <- Inf; if (max(v) <= 10) break
      dv <- names(which.max(v)); vif_removed <- c(vif_removed, lab_map[dv]); envs[[dv]] <- NULL
    }
    set.seed(42); perm <- input$env_perm
    glob <- tryCatch(anova(mod, permutations = perm), error = function(e) NULL)
    marg <- tryCatch(anova(mod, by = "margin", permutations = perm), error = function(e) NULL)
    vf <- tryCatch(vegan::vif.cca(mod), error = function(e) setNames(rep(NA_real_, ncol(envs)), names(envs)))
    ra <- vegan::RsquareAdj(mod); eigs <- c(mod$CCA$eig, mod$CA$eig); pct <- 100 * eigs[1:2] / mod$tot.chi
    sc <- scores(mod, display = c("sites", "species", "bp"), choices = 1:2, scaling = 2)
    sites <- two_cols(sc$sites); spec <- two_cols(sc$species); bp <- two_cols(sc$biplot); axn <- colnames(scores(mod, display = "sites", choices = 1:2))[1:2]
    colnames(sites) <- colnames(spec) <- colnames(bp) <- axn; rownames(bp) <- unname(lab_map[rownames(sc$biplot)])
    # correlations taxon ~ env (final variables, Hellinger abundance)
    envf <- env[, names(envs), drop = FALSE]; long <- NULL
    for (tx in colnames(Y)) for (vn in names(envf)) {
      ct <- tryCatch(cor.test(Y[, tx], envf[[vn]]), error = function(e) NULL)
      long <- rbind(long, data.frame(Taxon = tx, Variable = unname(lab_map[vn]), r = if (is.null(ct)) NA else unname(ct$estimate),
                                     p_value = if (is.null(ct)) NA else ct$p.value, stringsAsFactors = FALSE))
    }
    long$p_adj_BH <- p.adjust(long$p_value, "BH"); long$Stars <- stars_blank(long$p_adj_BH)
    stats_df <- data.frame(Variables_final = paste(unname(lab_map[names(envs)]), collapse = ", "),
                           Removed_by_collinearity = paste(dropped$Variable, collapse = ", "), Removed_by_VIF = paste(vif_removed, collapse = ", "),
                           R2 = round(unname(ra$r.squared), 4), Adj_R2 = round(unname(ra$adj.r.squared), 4),
                           Global_p = if (!is.null(glob)) glob$`Pr(>F)`[1] else NA, Axis1_pct = round(pct[1], 2), Axis2_pct = round(pct[2], 2),
                           Axis2_type = if (length(mod$CCA$eig) >= 2) "constrained (RDA2)" else "first unconstrained (PC1)",
                           row.names = NULL)
    list(mod = mod, sites = sites, spec = spec, bp = bp, axn = axn, pct = pct, mean_pct = mean_pct, long = long, lab_map = lab_map,
         stats = stats_df, glob = glob, marg = marg, vif = vf, cand = cand_df, dropped = dropped, p = p, level = input$env_level)
  })

  rda_gg <- reactive({
    r <- rda_run(); p <- r$p; ss <- as.data.frame(r$sites); ss$Sample <- rownames(ss)
    gv <- if (input$env_group == "None") NULL else grp_of(p, input$env_group); sv <- if (input$env_shape == "None") NULL else grp_of(p, input$env_shape)
    ss$G <- if (is.null(gv)) "samples" else gv[ss$Sample]; ss$S <- if (is.null(sv)) "samples" else sv[ss$Sample]
    A <- r$axn; names(ss)[1:2] <- c("X", "Y")
    lim <- max(abs(c(ss$X, ss$Y))) * 1.15
    mb <- 0.85 * max(abs(c(ss$X, ss$Y))) / max(abs(r$bp), 1e-9); bp <- as.data.frame(r$bp * mb); names(bp) <- c("X", "Y"); bp$lab <- rownames(bp)
    mt <- 0.9 * max(abs(c(ss$X, ss$Y))) / max(abs(r$spec), 1e-9); sp <- as.data.frame(r$spec * mt); names(sp) <- c("X", "Y"); sp$lab <- rownames(sp)
    sp$mean <- r$mean_pct[sp$lab]
    if (isTRUE(input$env_sig)) {
      L <- r$long; sig <- L[!is.na(L$p_value) & L$p_value <= 0.05, ]
      rmax <- tapply(abs(sig$r), sig$Taxon, max); keep <- names(sort(rmax, decreasing = TRUE))[seq_len(min(10, length(rmax)))]
      sp_plot <- sp[sp$lab %in% keep, , drop = FALSE]; sp_lab <- sp_plot
    } else {
      sp_plot <- sp; d0 <- sqrt(sp$X^2 + sp$Y^2); sp_lab <- sp[rank(-d0) <= (if (r$level == "ASV") 20 else 10), , drop = FALSE]
    }
    rep_t <- function(...) if (has("ggrepel")) ggrepel::geom_text_repel(...) else geom_text(..., vjust = -0.7)
    gl <- sort(unique(ss$G)); sl <- sort(unique(ss$S))
    g <- ggplot(ss, aes(X, Y)) + geom_hline(yintercept = 0, linetype = 2, color = "grey60") + geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
      geom_segment(data = bp, aes(x = 0, y = 0, xend = X, yend = Y), arrow = arrow(length = unit(0.10, "cm"), type = "closed"), color = "darkblue", linewidth = 0.4) +
      rep_t(data = bp, aes(X, Y, label = lab), color = "darkblue", size = 3.4) +
      geom_point(data = sp_plot, aes(X, Y), shape = 4, size = 1.4, color = "brown3", stroke = 0.5) +
      rep_t(data = sp_lab, aes(X, Y, label = lab), color = "brown3", size = 2.6) +
      geom_point(aes(color = G, shape = S), size = 3) +
      scale_color_manual(values = if (is.null(gv)) c(samples = "black") else group_cols(gl)) +
      scale_shape_manual(values = if (is.null(sv)) c(samples = 16) else rep(c(16, 17, 15, 18, 3, 4, 8, 7, 9, 10), 5)[seq_along(sl)]) +
      labs(x = sprintf("%s (%.1f%%)", A[1], r$pct[1]), y = sprintf("%s (%.1f%%)", A[2], r$pct[2]),
           color = if (is.null(gv)) NULL else input$env_group, shape = if (is.null(sv)) NULL else input$env_shape,
           title = sprintf("RDA - %s level, top %d taxa | adj R2 = %.3f | global p = %s", r$level, nrow(sp), r$stats$Adj_R2,
                           format.pval(r$stats$Global_p, digits = 3))) +
      theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11))
    if (is.null(gv)) g <- g + guides(color = "none"); if (is.null(sv)) g <- g + guides(shape = "none")
    g
  })
  output$rda_plot <- renderPlot(show_gg(rda_gg()))
  dl_gg("dl_rda", rda_gg, "rda_biplot")

  cor_gg <- reactive({
    L <- rda_run()$long; L$Taxon <- factor(L$Taxon, levels = rev(unique(L$Taxon)))
    ggplot(L, aes(Variable, Taxon, fill = r)) + geom_tile(color = "white") + geom_text(aes(label = Stars), size = 4) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Pearson r") +
      labs(x = NULL, y = NULL, caption = "* p_adj<=0.05  ** <=0.01  *** <=0.001 (BH)") + theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  output$cor_plot <- renderPlot(show_gg(cor_gg()))

  tbl_or_msg <- function(x) { validate(need(!is.null(x), "Not available.")); dt(x) }
  output$e_model <- renderDT(enc_safe(tbl_or_msg(rda_run()$stats)))
  output$e_glob  <- renderDT(enc_safe({ x <- rda_run()$glob; tbl_or_msg(if (is.null(x)) NULL else data.frame(Term = rownames(x), x, check.names = FALSE)) }))
  output$e_marg  <- renderDT(enc_safe({ x <- rda_run()$marg; tbl_or_msg(if (is.null(x)) NULL else data.frame(Term = rownames(x), x, check.names = FALSE)) }))
  output$e_vif   <- renderDT(enc_safe({ v <- rda_run()$vif; tbl_or_msg(data.frame(Variable = unname(rda_run()$lab_map[names(v)]), VIF = round(as.numeric(v), 3))) }))
  output$e_cand  <- renderDT(enc_safe(tbl_or_msg(rda_run()$cand)))
  cor_mats <- reactive({
    r <- rda_run(); L <- r$long; tt <- if (r$level == "ASV") tax_mat(r$p) else NULL
    mk <- function(col, f = identity) { w <- reshape(L[, c("Taxon", "Variable", col)], idvar = "Taxon", timevar = "Variable", direction = "wide")
      names(w) <- sub(paste0("^", col, "\\."), "", names(w)); w[-1] <- lapply(w[-1], f); rownames(w) <- NULL; w }
    rr <- mk("r", function(x) round(x, 2)); pp <- mk("p_value", function(x) round(x, 4)); pa <- mk("p_adj_BH", function(x) round(x, 4))
    rs <- rr; for (v in names(rr)[-1]) rs[[v]] <- paste0(rr[[v]], stars_blank(pa[[v]]))
    add <- function(w) if (is.null(tt)) w else cbind(w, tt[w$Taxon, , drop = FALSE])
    list(r = add(rs), p = add(pp), padj = add(pa))
  })
  output$e_cor <- renderDT(enc_safe(dt(cor_mats()$r)))
  dl_xl("xl_rda", function() { r <- rda_run(); cm <- cor_mats(); wrap <- function(x) if (is.null(x)) NULL else data.frame(Term = rownames(x), x, check.names = FALSE)
    list(Model_summary = r$stats, ANOVA_global = wrap(r$glob), ANOVA_margin = wrap(r$marg),
         VIF = data.frame(Variable = unname(r$lab_map[names(r$vif)]), VIF = round(as.numeric(r$vif), 3)),
         Correlation_screen = if (nrow(r$dropped)) r$dropped else data.frame(Note = "no variable removed"),
         Candidate_marginal_R2 = r$cand, Site_scores = data.frame(Sample = rownames(r$sites), r$sites, row.names = NULL),
         Taxa_scores = data.frame(Taxon = rownames(r$spec), r$spec, Mean_pct = round(r$mean_pct[rownames(r$spec)], 3), row.names = NULL),
         Env_arrows = data.frame(Variable = rownames(r$bp), r$bp, row.names = NULL),
         Corr_r_stars = cm$r, Corr_p_raw = cm$p, Corr_p_BH = cm$padj) }, "rda_results")

  # =================================================================
  # TAXA TABLES
  # =================================================================
  tax_df <- reactive({
    p <- ps_f(); cnt <- otu_mat(p); tt <- as.data.frame(tax_mat(p), stringsAsFactors = FALSE)
    tot <- rowSums(cnt)[taxa_names(p)]
    df <- data.frame(ASV = paste0("ASV", seq_len(ntaxa(p))), ID = taxa_names(p), tt, Total_reads = as.numeric(tot),
                     Rel_abundance_pct = round(100 * as.numeric(tot) / sum(tot), 3), Prevalence = rowSums(cnt > 0)[taxa_names(p)],
                     check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    rs <- refseq(p, errorIfNULL = FALSE); if (!is.null(rs)) df$Sequence <- as.character(rs)[taxa_names(p)]
    df[order(-df$Total_reads), ]
  })
  output$tax_tbl <- renderDT(enc_safe({
    df <- tax_df(); if ("Sequence" %in% names(df)) df$Sequence <- ifelse(nchar(df$Sequence) > 40, paste0(substr(df$Sequence, 1, 40), "..."), df$Sequence)
    datatable(df, rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15))
  }))
  output$dl_tax <- downloadHandler(filename = function() "taxa_table.csv",
                                   content = function(file) write.csv(tax_df(), file, row.names = FALSE, fileEncoding = "UTF-8"))

  pat_res <- reactive({
    pats <- trimws(unlist(strsplit(input$pat_text, "\n"))); pats <- pats[nzchar(pats)]
    validate(need(length(pats) > 0, "Enter one or more patterns (regular expressions), one per line."))
    p <- ps_f(); rel <- t(rel100(otu_mat(p))); tt <- tax_mat(p)[colnames(rel), , drop = FALSE]
    g <- grp_of(p, input$pat_group)
    do.call(rbind, lapply(pats, function(pt) {
      hit <- apply(tt, 1, function(r) any(grepl(pt, r, ignore.case = TRUE))); if (!any(hit)) return(NULL)
      ps <- rowSums(rel[, hit, drop = FALSE])
      d <- data.frame(Pattern = pt, ASVs = sum(hit), Overall_pct = round(mean(ps), 2), stringsAsFactors = FALSE)
      if (!is.null(g)) for (lv in unique(g[names(ps)])) d[[paste0(gsub("[^A-Za-z0-9]+", "_", lv), "_pct")]] <- round(mean(ps[g[names(ps)] == lv]), 2)
      d$Top_sample <- names(ps)[which.max(ps)]; d$Top_sample_pct <- round(max(ps), 2); d
    }))
  })
  output$pat_tbl <- renderDT(enc_safe({ x <- pat_res(); validate(need(!is.null(x), "No taxonomy matches those patterns.")); dt(x) }))
  output$dl_pat <- downloadHandler(filename = function() "group_search.csv",
                                   content = function(file) write.csv(pat_res(), file, row.names = FALSE, fileEncoding = "UTF-8"))

  output$dl_ps <- downloadHandler(filename = function() "filtered.phyloseq.rds", content = function(file) saveRDS(ps_f(), file))
}
