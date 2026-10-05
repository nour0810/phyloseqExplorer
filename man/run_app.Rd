\name{run_app}
\alias{run_app}
\title{Launch the phyloseq Explorer app}
\usage{
run_app(launch.browser = TRUE, ...)
}
\arguments{
  \item{launch.browser}{Open the app in the browser? Default \code{TRUE}.}
  \item{...}{Additional arguments passed to \code{\link[shiny]{runApp}}.}
}
\value{
Called for its side effect; returns the \code{\link[shiny]{runApp}} result.
}
\description{
Opens the interactive Shiny application for exploring any \code{phyloseq}
object: sequence/QC overview, rarefaction curves, composition plots and
tables, alpha diversity with assumption-checked statistics, beta diversity
(PCoA, dendrogram, biplot, PERMANOVA, betadisper), constrained ordination
(RDA) with collinearity and VIF guards, taxon-environment correlations and
taxonomy tables. Figures export as SVG, PNG, PDF or TIFF; every tab exports
Excel workbooks. All processing happens locally.
}
\examples{
\dontrun{
phyloseqExplorer::run_app()
}
}
