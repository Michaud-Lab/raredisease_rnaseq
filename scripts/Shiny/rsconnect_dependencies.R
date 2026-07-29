# Not sourced by the app. Packages below are loaded dynamically at runtime via
# load_install_library() in rnaseq_helper_functions.R, so rsconnect's static
# dependency scan (used by rsconnect::writeManifest()) can't see them from a
# character vector. Listing them here as literal library() calls makes them
# show up in manifest.json's package list.
if (FALSE) {
  library(reactable)
  library(plotly)
  library(shinyjs)
  library(igvShiny)
  library(GenomicAlignments)
  library(ggtranscript)
  library(patchwork)
  library(Hmisc)
  library(ggrepel)
  library(rtracklayer)
  library(GenomeInfoDb)
}
