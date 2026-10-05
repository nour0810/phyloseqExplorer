# phyloseqExplorer

An interactive R/Shiny application for **microbiome & metabarcoding data stored as
[phyloseq](https://bioconductor.org/packages/phyloseq) objects** — from raw QC to
publication-ready figures and statistics, fully offline.

## Features

| Tab | Contents |
|---|---|
| **Overview** | Object summary, sequence statistics, library sizes, smooth rarefaction curves, sample metadata |
| **Composition** | Stacked bars (cumulative % cutoff or Top-N with "Others"), treemaps, per-clade table with **n ASVs, ASV %, n reads, reads %, cumulative reads %** |
| **Heatmap** | Z-score heatmap (Ward.D2 clustering or fixed group blocks) with row labels `ASV` / `ASV_species` / `ASV_species_class` / custom ranks |
| **Alpha diversity** | Observed, Simpson, Chao1, Shannon (+ Faith's PD when a tree exists); Shapiro-Wilk + Bartlett assumption checks; ANOVA/Tukey or Kruskal-Wallis/pairwise Wilcoxon (BH) with significance brackets |
| **Beta diversity** | Bray-Curtis (or Jaccard/Euclidean) dendrogram, PCoA + biplot; PERMANOVA, betadisper, pairwise PERMANOVA (BH); method summary table |
| **Environment (RDA)** | Constrained ordination with collinearity screen (\|r\| threshold) and VIF guard; taxon–environment correlation heatmap with BH-adjusted significance |
| **Taxa tables** | Full ASV table (taxonomy, reads, prevalence, sequence) + regex group search |

**Exports:** every figure as SVG / PNG / PDF / **TIFF (LZW, 300–600 dpi)**; every tab
as a formatted Excel workbook (frozen headers, auto widths); and the filtered phyloseq
object as `.rds`.

**Design goals:** works with *any* phyloseq object (rank/group selectors auto-populate),
repairs encoding issues in sample/taxonomy names automatically, needs zero server —
your data never leaves your machine.

## Installation

```r
# from GitHub (once published):
remotes::install_github("nour0810/phyloseqExplorer")
```

For local development, clone the repo and use:

```r
devtools::load_all()   # or: install.packages(".", repos = NULL, type = "source")
```

## Usage

```r
phyloseqExplorer::run_app()
```

Click **Load demo (GlobalPatterns)** to explore with the built-in phyloseq dataset,
or upload your own `.rds` / `.RData` phyloseq object. Optional metadata (csv/tsv/xlsx)
and group-from-sample-name rules can be added in the sidebar.

Optional packages unlock extra features (the app tells you what to install if missing):
`ggpubr` (p-value brackets), `ggdendro` (ggplot dendrogram), `treemap`, `pheatmap`,
`ggrepel` (label repulsion), `svglite` (better SVG), `patchwork` (biplot panels),
`picante` (Faith's PD).

## Methods & reproducibility notes

* All random procedures use `set.seed(42)`; record your settings (filters, cutoffs,
  permutations) in your lab notes or methods section when publishing.
* Alpha diversity is computed after rarefying all samples to the minimum library
  depth — used here *only* to compare richness indices on equal footing (see
  McMurdie & Holmes 2014, *PLoS Comput Biol*, for why rarefying is inadmissible in
  count models). Unrarefied counts are retained everywhere else.
* Alpha group tests are assumption-checked per index: Shapiro–Wilk (normality) and
  Bartlett (variance); parametric (ANOVA + Tukey HSD) or non-parametric
  (Kruskal–Wallis + pairwise Wilcoxon, BH) accordingly.
* Beta diversity: PERMANOVA (`vegan::adonis2`) with a homogeneity-of-dispersion
  check (`betadisper`) and BH-adjusted pairwise tests.
* RDA uses Hellinger-transformed abundances with automated collinearity screening
  and a VIF ≤ 10 guard.

## Screenshots

<!-- TODO: add one screenshot per tab, e.g. stored in man/figures/ -->
<!-- ![Composition](man/figures/composition.png) -->

## Known limitations (roadmap)

* Very large datasets (hundreds of samples) recompute all tabs on every filter
  change; caching (`shiny::bindCache`) is planned.
* The phylogeny panel, differential abundance and indicator-species analyses are
  planned for a future release.

## License

MIT. Citation: forthcoming (Zenodo DOI will be minted on the first release).
