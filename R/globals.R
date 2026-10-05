# Quiet "no visible binding" notes for column names used inside aes()/data.frames.
# These are data-frame columns, not real globals.
utils::globalVariables(c(
  "Abundance", "Axis1", "Axis2", "Depth", "F1", "F2", "G", "Group", "Name",
  "PCo1", "PCo2", "Rank", "Reads", "Richness", "S", "Sample", "Stars",
  "Tax_plot", "Taxon", "Variable", "label", "value", "x", "xend", "y", "yend"
))
