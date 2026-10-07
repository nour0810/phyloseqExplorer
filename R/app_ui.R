# ---------------------------------------------------------------------
# Application UI (wrapped as a function so it works inside the package)
# ---------------------------------------------------------------------

#' Application UI (internal)
#' @noRd
app_ui <- function(request) {
fluidPage(
  tags$head(tags$style(HTML(
    "#title-bar { position: fixed; top: 0; left: 0; right: 0; z-index: 1000; background: #fff; padding: 6px 15px 4px 15px; border-bottom: 1px solid #ddd; } body { padding-top: 82px; }"
  ))),
  div(id = "title-bar",
      titlePanel(
        div(
          "phyloseq Explorer",
          div(style = "font-size:16px; color:#444; font-weight:normal; margin-top:3px;",
              "\u00A9 2026 Mathlouthi Nourelhouda (TN)  ",
              span(style = "font-size:16px;", "\U0001F1F9\U0001F1F3"))
        )
      )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("file", "phyloseq object (.phyloseq / .rds / .RData)"),
      actionButton("demo", "Load demo (GlobalPatterns)", class = "btn-sm"),
      tags$details(
        tags$summary(strong("Metadata (optional)")),
        fileInput("meta_file", "Upload metadata table (csv / tsv / xlsx; 1st column = sample names)"),
        textInput("rule_name", "New group variable from sample names", "Group"),
        textAreaInput("rule_text", "One rule per line:  pattern = label  (first match wins)",
                      rows = 4, placeholder = "B200 = Large\nB60 = Medium\nMCL = Small"),
        actionButton("apply_rules", "Apply rules", class = "btn-sm")
      ),
      tags$details(
        open = NA,
        tags$summary(strong("Filters (all tabs)")),
        checkboxGroupInput("samples", "Samples", choices = NULL),
        numericInput("min_depth", "Min reads per sample", 1, min = 0),
        numericInput("min_abund", "Min total reads per taxon", 1, min = 0),
        numericInput("min_prev", "Min samples a taxon appears in", 1, min = 1),
        selectInput("filt_rank", "Keep only taxa where rank...", choices = "None"),
        selectizeInput("filt_vals", "...is one of (empty = no filter)", choices = NULL, multiple = TRUE),
        selectInput("filt_rank2", "Remove taxa where rank...", choices = "None"),
        selectizeInput("filt_vals2", "...is one of (excluded)", choices = NULL, multiple = TRUE)
      ),
      tags$details(
        tags$summary(strong("Palette & export")),
        selectInput("pal", "Taxon palette",
                    c("Set1-seeded max-min (60)" = "set1_60", "Contrasting (30)" = "contrast",
                      "ggplot hue" = "hue", "Viridis" = "viridis")),
        radioButtons("fmt", "Figure format", c("svg", "png", "pdf", "tiff"), inline = TRUE),
        radioButtons("padj_mode", "P-values in tests", c("BH-adjusted" = "bh", "Raw p" = "raw"), inline = TRUE),
        numericInput("ex_w", "Width (in)", 12, min = 3), numericInput("ex_h", "Height (in)", 6, min = 3),
        numericInput("ex_dpi", "DPI (png)", 300, min = 72)
      ),
      hr(), downloadButton("dl_ps", "Filtered phyloseq (.rds)", class = "btn-sm")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        # ---------------- Overview
        tabPanel("Overview", br(),
          verbatimTextOutput("ps_print"),
          fluidRow(column(6, h5("Sequence statistics"), DTOutput("seq_stats")),
                   column(6, h5("Unique taxa per rank (excl. Unknown)"), DTOutput("rank_census"))),
          h5("Library sizes"), plotOutput("depth_plot", height = "280px"),
          h5("Rarefaction curves"),
          checkboxInput("rc_clip", "Clip x axis at the smallest library", FALSE),
          plotOutput("rare_plot", height = "380px"),
          h5("Sample metadata"), DTOutput("meta_tbl"), br(),
          dl_row("dl_rare", "xl_over")),
        # ---------------- Composition
        tabPanel("Composition", br(),
          fluidRow(
            column(3, selectInput("rank", "Taxonomic rank", NULL)),
            column(3, radioButtons("comp_mode", "Show", c("Cumulative % cutoff" = "cum", "Top N" = "top"), inline = TRUE)),
            column(2, numericInput("comp_cum", "Cutoff (%)", 99, 50, 100)),
            column(2, numericInput("topn", "Top N", 20, 1, 80)),
            column(2, selectInput("facet", "Facet by", "None"))),
          plotOutput("bar_plot", height = "560px"),
          dl_row("dl_bar", "xl_bar"),
          fluidRow(column(5, h5("Clades shown (printed table of your script)"), DTOutput("cum_tbl")),
                   column(7, h5("Mean % per group"), DTOutput("grp_tbl")))),
        # ---------------- Treemap
        tabPanel("Treemap", br(),
          fluidRow(column(3, selectInput("tm_rank", "Rank", NULL)),
                   column(3, selectInput("tm_group", "Panels (group)", "None")),
                   column(3, numericInput("tm_thr", "Merge taxa below (%)", 1, 0, 20, step = 0.5))),
          if (has("treemap")) plotOutput("tm_plot", height = "520px") else
            div(class = "alert alert-warning", "Install the 'treemap' package: install.packages('treemap')"),
          dl_row("dl_tm", "xl_tm")),
        # ---------------- Heatmap
        tabPanel("Heatmap", br(),
          fluidRow(column(2, selectInput("hm_level", "Rows", "ASV")),
                   column(2, numericInput("hm_n", "Top N", 40, 3, 150)),
                   column(2, selectInput("hm_group", "Column groups", "None")),
                   column(2, numericInput("hm_cap", "Z-score cap", 1.5, 0.5, 5, step = 0.5)),
                   column(2, selectInput("hm_lab", "Row label",
                                         c("ASV only" = "asv", "ASV_species" = "asv_sp",
                                           "ASV_species_class" = "asv_sp_cl", "Custom (ranks below)" = "custom")))),
          fluidRow(column(4, radioButtons("hm_cols", "Columns",
                                          c("Fixed in group blocks" = "blocks", "Cluster (Ward.D2)" = "clust"), inline = TRUE)),
                   column(2, selectInput("hm_l1", "Custom rank 1", "(none)")),
                   column(2, selectInput("hm_l2", "Custom rank 2", "(none)"))),
          if (has("pheatmap")) plotOutput("hm_plot", height = "700px") else
            div(class = "alert alert-warning", "Install 'pheatmap': install.packages('pheatmap')"),
          dl_row("dl_hm", "xl_hm")),
        # ---------------- Alpha
        tabPanel("Alpha diversity", br(),
          fluidRow(column(3, selectInput("a_group", "Group by", "Sample")),
                   column(4, checkboxGroupInput("a_meas", "Indices", c("Observed", "Simpson", "Chao1", "Shannon"),
                                                selected = c("Observed", "Simpson", "Chao1", "Shannon"), inline = TRUE)),
                   column(3, checkboxInput("a_rare", "Rarefy to min depth (seed 42)", TRUE)),
                   column(2, radioButtons("a_test", "Figure test", c("Auto" = "auto", "ANOVA/t" = "p", "Kruskal/Wilcoxon" = "np")))),
          uiOutput("a_note"),
          h5("Per-sample values"), plotOutput("a_dot", height = "420px"),
          h5("Group comparison"), plotOutput("a_box", height = "420px"),
          dl_row("dl_alpha", "xl_alpha"),
          tabsetPanel(
            tabPanel("Values", DTOutput("a_vals")), tabPanel("Normality", DTOutput("a_norm")),
            tabPanel("Variance", DTOutput("a_var")), tabPanel("Global test", DTOutput("a_glob")),
            tabPanel("Post-hoc", DTOutput("a_post")))),
        # ---------------- Beta
        tabPanel("Beta diversity", br(),
          fluidRow(column(2, selectInput("b_group", "Color / test factor", "Sample")),
                   column(2, selectInput("b_group2", "Shape / 2nd factor", "None")),
                   column(2, numericInput("b_cum", "ASV cumulative cutoff (%)", 95, 50, 100)),
                   column(2, selectInput("b_trans", "Transform", c("sqrt" = "sqrt", "none" = "none", "relative" = "rel", "Hellinger" = "hell"))),
                   column(2, selectInput("b_dist", "Distance", c("bray", "jaccard", "euclidean"))),
                   column(2, numericInput("b_perm", "Permutations", 999, 99, 9999, step = 100))),
          fluidRow(column(3, checkboxInput("b_ell", "Ellipses", TRUE)),
                   column(3, selectInput("b_taxrank", "Biplot taxa colored by", NULL))),
          h5("A - Dendrogram"), plotOutput("b_dend", height = "380px"),
          fluidRow(column(5, h5("B - PCoA"), plotOutput("b_pcoa", height = "430px")),
                   column(7, h5("C - PCoA biplot (samples | taxa)"), plotOutput("b_bi", height = "430px"))),
          dl_row("dl_beta", "xl_beta"),
          tabsetPanel(
            tabPanel("PERMANOVA", DTOutput("b_perm_t")), tabPanel("Betadisper", DTOutput("b_disp_t")),
            tabPanel("Pairwise PERMANOVA", DTOutput("b_pair_t")), tabPanel("Method", DTOutput("b_meth_t")))),
        # ---------------- Network
        tabPanel("Network", br(),
          if (has("igraph")) tagList(
            fluidRow(column(2, selectInput("net_mode", "Mode", c("Sample similarity (A)" = "sample", "Taxon co-occurrence (B)" = "cooc")))),
            fluidRow(
              column(3, fileInput("net_file2", "Second dataset (optional: .rds/.RData)")),
              column(2, selectInput("net_rank2", "2nd: rank", NULL)),
              column(4, selectizeInput("net_tax2", "2nd: taxa to include (pick from list)", NULL, multiple = TRUE)),
              column(2, numericInput("net_min2", "2nd: min reads", 100, 0, 1e9)),
              column(1, checkboxInput("net_shape2", "Triangles = 2nd", TRUE))),
            fluidRow(
              column(2, selectInput("net_dist", "Distance (A)", c("bray", "jaccard", "euclidean", "manhattan", "canberra"))),
              column(2, sliderInput("net_maxd", "Edge: max dist (A)", 0.05, 0.95, 0.4, 0.05)),
              column(2, checkboxInput("net_isol", "Keep isolated samples (A)", TRUE)),
              column(2, selectInput("net_col", "Node color (A)", "Sample")),
              column(2, selectInput("net_shp", "Node shape (A)", "None")),
              column(2, radioButtons("net_nsize", "Node size (A)", c("Reads" = "reads", "Uniform" = "uni"), inline = TRUE))),
            fluidRow(
              column(2, textInput("net_lay", "Layout (fr/kk/circle/graphopt/mds/grid)", "graphopt")),
              column(2, numericInput("net_seed", "Layout seed", 42, 1, 1e6)),
              column(2, checkboxInput("net_lbl", "Label samples (A)", TRUE)),
              column(2, numericInput("net_lbls", "Label size (A)", 3, 1, 8)),
              column(2, numericInput("net_esz", "Edge width (A)", 0.6, 0.1, 3, step = 0.1)),
              column(2, numericInput("net_alpha", "Edge alpha (A)", 0.7, 0.1, 1, step = 0.1))),
            tags$hr(),
            fluidRow(
              column(2, selectInput("net_rank", "Node rank (B)", NULL)),
              column(2, numericInput("net_topn", "Top N taxa (B)", 60, 10, 300)),
              column(2, numericInput("net_mina", "Min rel. % (B)", 0.01, 0, 1, step = 0.01)),
              column(2, numericInput("net_minp", "Min prev. % (B)", 10, 0, 100)),
              column(2, selectInput("net_corm", "Correlation (B)", c("spearman", "pearson", "kendall"))),
              column(2, sliderInput("net_r", "Min |r| (B)", 0.2, 0.95, 0.6, 0.05))),
            fluidRow(
              column(2, numericInput("net_minpct", "Min total reads %", 0, 0, 100, step = 0.01)),
              column(2, numericInput("net_minreads", "Min abundance (reads)", 0, 0, 1e9)),
              column(2, selectInput("net_puse", "Filter edges by p", c("BH-adjusted" = "bh", "Raw p" = "raw"))),
              column(6, helpText("Min total reads % / Min abundance: total-contribution filters (before Top N). 'Filter edges by p' decides which p-value the cutoff slider below applies to (overrides the sidebar default for this tab)."))),
            fluidRow(
              column(2, sliderInput("net_fdr", "FDR cutoff (B)", 0.001, 0.2, 0.05, 0.001)),
              column(2, selectInput("net_edge", "Edges (B)", c("Both signs" = "both", "Positive" = "pos", "Negative" = "neg"))),
              column(2, numericInput("net_deg", "Min degree (B)", 1, 0, 20)),
              column(2, radioButtons("net_nsizeb", "Node size (B)", c("Reads" = "reads", "Degree" = "deg"), inline = TRUE)),
              column(2, numericInput("net_hub", "Label top N hubs (B)", 10, 0, 50)),
              column(2, selectInput("net_colrank", "Node color rank (B)", NULL))),
            tags$hr(),
            tags$small(strong("Figure style (B)")),
            fluidRow(
              column(2, radioButtons("net_nfill", "Node fill", c("By rank" = "rank", "Single" = "one"), inline = TRUE)),
              column(2, textInput("net_fcol", "Node hex", "#F0653A")),
              column(2, selectInput("net_ecol", "Edge color", c("By sign" = "sign", "Custom" = "one"))),
              column(2, textInput("net_ecolhex", "Edge hex", "#5BA8A0")),
              column(2, selectInput("net_lfmt", "Label", c("Name" = "name", "Name - rank" = "both", "Top ASV - name" = "asv"))),
              column(2, checkboxInput("net_leg", "Legend", TRUE))),
            fluidRow(
              column(2, sliderInput("net_nscale", "Node size x", 0.5, 5, 1, 0.1)),
              column(2, sliderInput("net_nstroke", "Node border", 0, 2, 0.6, 0.1)),
              column(2, numericInput("net_lbls_b", "Label size", 3.2, 1, 8)),
              column(2, sliderInput("net_ewb", "Edge width x", 0.2, 3, 1, 0.1)),
              column(2, sliderInput("net_eab", "Edge alpha", 0.1, 1, 0.6, 0.05)),
              column(2, br())),
            fluidRow(
              column(2, sliderInput("net_curv", "Edge curve", 0, 0.5, 0.15, 0.05)),
              column(2, checkboxInput("net_bold", "Bold labels", TRUE)),
              column(2, checkboxInput("net_eq", "Equal aspect", TRUE)),
              column(2, textInput("net_title", "Custom title", "")),
              column(6, helpText("Custom title replaces the stats line, e.g.: Identifying nodes with degree > 3 at the genus level (all connections were positive)"))),
            plotOutput("net_plot", height = "620px"),
            dl_row("dl_net", "xl_net"),
            tabsetPanel(
              tabPanel("Edges", DTOutput("net_edge_t")), tabPanel("Nodes", DTOutput("net_node_t"))))
          else div(class = "alert alert-warning", "Install 'igraph' for network analysis: install.packages('igraph'")),
        # ---------------- Environment
        tabPanel("Environment (RDA)", br(),
          fluidRow(column(4, selectizeInput("env_vars", "Numeric environmental variables (from metadata)", NULL, multiple = TRUE)),
                   column(2, selectInput("env_level", "Taxa level", "ASV")),
                   column(2, numericInput("env_topn", "Top N taxa", 20, 3, 300)),
                   column(2, numericInput("env_r", "Collinearity |r| >=", 0.9, 0.5, 1, step = 0.05)),
                   column(2, checkboxInput("env_vif", "VIF guard (<= 10)", TRUE))),
          fluidRow(column(3, selectInput("env_group", "Color sites by", "None")),
                   column(3, selectInput("env_shape", "Shape sites by", "None")),
                   column(3, checkboxInput("env_sig", "Show only taxa with raw p <= 0.05 (max 10 by |r|)", FALSE)),
                   column(3, numericInput("env_perm", "Permutations", 999, 99, 9999, step = 100))),
          plotOutput("rda_plot", height = "560px"),
          h5("Taxon ~ environment correlations (r, stars from BH-adjusted p)"),
          plotOutput("cor_plot", height = "420px"),
          dl_row("dl_rda", "xl_rda"),
          tabsetPanel(
            tabPanel("Model", DTOutput("e_model")), tabPanel("ANOVA global", DTOutput("e_glob")),
            tabPanel("ANOVA margin", DTOutput("e_marg")), tabPanel("VIF", DTOutput("e_vif")),
            tabPanel("Candidate R2", DTOutput("e_cand")), tabPanel("Correlations (r + stars)", DTOutput("e_cor")),
            tabPanel("Correlations (p-values)", DTOutput("e_cor_p")))),
        # ---------------- Tables
        tabPanel("Taxa tables", br(),
          h5("ASV table (taxonomy, reads, prevalence, sequence)"),
          downloadButton("dl_tax", "CSV", class = "btn-sm"), DTOutput("tax_tbl"), hr(),
          h5("Group search: % of reads for any taxon pattern (regex, one per line)"),
          fluidRow(column(4, textAreaInput("pat_text", NULL, rows = 6,
                                           placeholder = "Alveolata\nCopepoda\nSyndiniales.?I([^I]|$)")),
                   column(3, selectInput("pat_group", "Summarise by", "None")),
                   column(2, br(), downloadButton("dl_pat", "CSV", class = "btn-sm"))),
          DTOutput("pat_tbl"))
      )
    )
  )
)
}
