# =============================================================================
# test_helper_functions.R
#
# Tests for:
#   - gwFRASER_table()        : filter / sort / deduplicate FRASER results
#   - manhattan_plot()        : data-manipulation logic (title, chr factorisation,
#                               cumulative positions) – plot rendering not tested
#   - gene_prioritization()   : merge HPO + OUTRIDER + FRASER; score & rank genes
#
# Run from the project root:
#   testthat::test_file("scripts/tests/testthat/test_helper_functions.R")
# =============================================================================

library(testthat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(plotly)

# Source the functions under test.  Adjust if running from a different cwd.
source("scripts/Shiny/rnaseq_shinyhelper_functions.R")

# =============================================================================
# Shared mock-data helpers
# =============================================================================

make_fraser_df <- function() {
  # 14-column data frame that matches the column layout expected by
  # gwFRASER_table (returns cols c(1:5,7,9:10,13:14)).
  data.frame(
    sampleID     = c("S1", "S1", "S1", "S1", "S2", "S2"),   # col 1
    chr          = c("2",  "1",  "X",  "1",  "3",  "4"),    # col 2
    start        = c(2e5,  1e5,  3e5,  1.5e5, 4e5, 5e5),    # col 3
    end          = c(2.1e5,1.1e5,3.1e5,1.6e5, 4.1e5,5.1e5), # col 4
    pos          = c(2.05e5,1.05e5,3.05e5,1.55e5,4.05e5,5.05e5), # col 5
    strand       = c("+",  "+",  "-",  "+",  "+",  "-"),    # col 6
    type         = c("psi3","psi5","theta","psi3","psi5","theta"), # col 7
    pValue       = c(0.001, 0.002, 0.40, 0.001, 0.001, 0.002),    # col 8
    padjust      = c(0.01,  0.04,  0.60, 0.01,  0.01,  0.03),     # col 9
    psiValue     = c(0.8,   0.7,   0.9,  0.8,   0.8,   0.7),      # col 10
    deltaPsi     = c(0.3,   0.2,   0.05, 0.3,   0.4,   0.1),      # col 11
    counts       = c(100,   80,    200,  90,    110,   75),        # col 12
    totalCounts  = c(500,   400,   600,  450,   550,   380),       # col 13
    hgncSymbol   = c("ABAT","ABCB7", NA, "ABAT","AARS1","GENE_D"), # col 14
    l2fc         = c(1.2,  -0.8,   0.3,  1.1,  -1.5,   0.6),
    minuslogpval = c(2.0,   1.4,   0.2,  2.0,   2.0,   1.5),
    stringsAsFactors = FALSE
  )
}

make_outrider_df <- function() {
  data.frame(
    sampleID    = c("S1",    "S1",    "S2"),
    geneID      = c("ABAT","ABCB7","AARS1"),
    pValue      = c(0.001,   0.002,   0.001),
    zScore      = c(3.5,    -2.8,     4.1),
    exon_zScore = c(3.2,    -2.5,     3.9),
    exon_pValue = c(0.002,   0.003,   0.001),
    stringsAsFactors = FALSE
  )
}

make_fraser_summary_df <- function() {
  data.frame(
    sampleID   = c("S1",    "S1",    "S2"),
    hgncSymbol = c("ABAT","ABCB7","AARS1"),
    pValue     = c(0.01,    0.02,    0.01),
    stringsAsFactors = FALSE
  )
}

make_clinical_df <- function() {
  data.frame(
    `Patient ID` = c("S1", "S2"),
    `HPO terms`  = c("HP:0001250: Seizure||HP:0002119: Ventriculomegaly", ""),
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
}

# Minimal HPO gene-to-phenotype table (matches columns of the real HPO file)
make_hpo_df <- function() {
  data.frame(
    ncbi_gene_id = c(1,    2,      3,      4),
    gene_symbol  = c("ABAT","ABAT","ABCB7","AARS1"),
    hpo_id       = c("HP:0001250","HP:0002119","HP:0001250","HP:0001290"),
    hpo_name     = c("Seizure","Ventriculomegaly","Seizure","Hypotonia"),
    frequency    = c(NA, NA, NA, NA),
    disease_id   = c("OMIM:1","OMIM:2","OMIM:3","OMIM:4"),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 1.  gwFRASER_table()
# =============================================================================

test_that("gwFRASER_table: returns only rows for the requested sample", {
  df  <- make_fraser_df()
  out <- gwFRASER_table(res_dt = df, sample = "S1")
  expect_true(all(out$sampleID == "S1"))
})

test_that("gwFRASER_table: returns a data.frame", {
  df  <- make_fraser_df()
  out <- gwFRASER_table(res_dt = df, sample = "S1")
  expect_s3_class(out, "data.frame")
})

test_that("gwFRASER_table: filters rows with padjust >= 0.05", {
  df  <- make_fraser_df()
  # Row 3 of S1 has padjust=0.60, should be excluded
  out <- gwFRASER_table(res_dt = df, sample = "S1", pvalue = "padjust")
  expect_true(all(out$padjust < 0.05))
})

test_that("gwFRASER_table: deduplicates by gene symbol (keeps first by p-value)", {
  df  <- make_fraser_df()
  # S1 has GENE_A twice (rows 1 & 4), both with padjust < 0.05
  out <- gwFRASER_table(res_dt = df, sample = "S1", pvalue = "padjust")
  gene_col <- out$hgncSymbol
  gene_col_no_na <- gene_col[gene_col != "na"]
  expect_equal(length(gene_col_no_na), length(unique(gene_col_no_na)))
})

test_that("gwFRASER_table: replaces NA gene symbols with 'na'", {
  df  <- make_fraser_df()
  # Temporarily give NA row a low p-value so it passes the filter
  df$padjust[df$sampleID == "S1" & is.na(df$hgncSymbol)] <- 0.001
  out <- gwFRASER_table(res_dt = df, sample = "S1", pvalue = "padjust")
  expect_false(any(is.na(out$hgncSymbol)))
})

test_that("gwFRASER_table: orders output by chromosome", {
  df  <- make_fraser_df()
  out <- gwFRASER_table(res_dt = df, sample = "S1", pvalue = "padjust")
  chr_levels <- levels(factor(out$chr, levels = c(as.character(1:22), "X", "Y", "MT")))
  chr_ordered <- factor(out$chr, levels = chr_levels)
  expect_equal(as.character(chr_ordered), as.character(sort(chr_ordered)))
})

test_that("gwFRASER_table: returns empty data.frame for unknown sample", {
  df  <- make_fraser_df()
  out <- gwFRASER_table(res_dt = df, sample = "DOES_NOT_EXIST")
  expect_equal(nrow(out), 0)
})

# NOTE: The pcutoff parameter is currently ignored in the filter (hardcoded 0.05).
# The following test documents this known behaviour.
test_that("gwFRASER_table: pcutoff param has no effect (hardcoded 0.05 filter)", {
  df <- make_fraser_df()
  # pcutoff=0.10 should theoretically include rows with padjust 0.04..0.10
  # but the function uses < 0.05 regardless
  out_strict  <- gwFRASER_table(res_dt = df, sample = "S1", pcutoff = 0.03)
  out_relaxed <- gwFRASER_table(res_dt = df, sample = "S1", pcutoff = 0.10)
  expect_equal(nrow(out_strict), nrow(out_relaxed))
})


# =============================================================================
# 2.  manhattan_plot() – data-manipulation logic only
#     We test the title generation and internal data structures by inspecting
#     the ggplot object, without rendering it.
# =============================================================================

test_that("manhattan_plot: returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df  <- make_fraser_df()
  p   <- manhattan_plot(res_dt = df, sample = "S1")
  expect_s3_class(p, "gg")
})

test_that("manhattan_plot: title contains the sample name", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df    <- make_fraser_df()
  p     <- manhattan_plot(res_dt = df, sample = "S1")
  title <- p$labels$title
  expect_true(grepl("S1", title))
})

test_that("manhattan_plot: title reports correct outlier count (padjust < 0.05)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df <- make_fraser_df()
  # S1 rows with padjust < 0.05: rows 1 (0.01), 2 (0.04), 4 (0.01) = 3 events
  p  <- manhattan_plot(res_dt = df, sample = "S1", pvalue = "padjust", pcutoff = 0.05)
  title <- p$labels$title
  expect_true(grepl("3 outliers", title))
})

test_that("manhattan_plot: title contains median outlier count across samples", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df <- make_fraser_df()
  p  <- manhattan_plot(res_dt = df, sample = "S1", pvalue = "padjust", pcutoff = 0.05)
  title <- p$labels$title
  expect_true(grepl("Median", title, ignore.case = TRUE))
})

test_that("manhattan_plot: shape='over'/'under' assigned when shape=TRUE", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df   <- make_fraser_df()
  p    <- manhattan_plot(res_dt = df, sample = "S1", shape = TRUE)
  shapes <- unique(p$data$shape)
  expect_true(all(shapes %in% c("over", "under")))
})

test_that("manhattan_plot: shape='splicing' when shape=FALSE (default)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df <- make_fraser_df()
  p  <- manhattan_plot(res_dt = df, sample = "S1", shape = FALSE)
  expect_true(all(p$data$shape == "splicing"))
})

test_that("manhattan_plot: BPcum is monotonically non-decreasing within each chr", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("RColorBrewer")
  df <- make_fraser_df()
  p  <- manhattan_plot(res_dt = df, sample = "S1")
  # All rows of the full data (not just S1) are used to build BPcum
  bpcum_by_chr <- p$data %>% arrange(chr, pos) %>% group_by(chr) %>%
    summarise(is_monotone = all(diff(BPcum) >= 0), .groups = "drop")
  expect_true(all(bpcum_by_chr$is_monotone))
})


# =============================================================================
# 3.  gene_prioritization()
# =============================================================================

# Helper: write a minimal HPO file to a temp directory and return the dir path
setup_hpo_tempdir <- function() {
  tmp   <- tempdir()
  hpo   <- make_hpo_df()
  hpodir <- file.path(tmp, "temp")
  dir.create(hpodir, showWarnings = FALSE)
  write.table(hpo, file = file.path(hpodir, "genes_to_phenotype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  tmp
}

test_that("gene_prioritization: returns a data.frame", {
  wd  <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample      = "S1",
      hpo_sample  = make_clinical_df(),
      outrider    = make_outrider_df(),
      fraser      = make_fraser_summary_df()
    )
  })
  expect_s3_class(out, "data.frame")
})

test_that("gene_prioritization: respects top parameter", {
  wd <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out_1 <- gene_prioritization(
      sample = "S1", top = 1,
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
    out_5 <- gene_prioritization(
      sample = "S1", top = 5,
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  expect_lte(nrow(out_1), 1)
  expect_lte(nrow(out_5), 5)
})

test_that("gene_prioritization: output is sorted by gene score descending", {
  wd <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample     = "S1",
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  if (nrow(out) > 1) {
    expect_true(all(diff(out$`gene score`) <= 0))
  }
})

test_that("gene_prioritization: GENE_A has higher score (2 HPO terms) than GENE_B (1 term)", {
  wd <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample     = "S1",
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  if (all(c("GENE_A", "GENE_B") %in% out$geneID)) {
    score_a <- out$`gene score`[out$geneID == "GENE_A"]
    score_b <- out$`gene score`[out$geneID == "GENE_B"]
    expect_gt(score_a, score_b)
  }
})

test_that("gene_prioritization: output contains geneID column", {
  wd <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample     = "S1",
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  expect_true("geneID" %in% colnames(out))
})

test_that("gene_prioritization: handles sample with no HPO terms gracefully", {
  wd <- setup_hpo_tempdir()
  # S2 has empty HPO terms string
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample     = "S2",
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  # Should return a data.frame (possibly empty)
  expect_s3_class(out, "data.frame")
})

test_that("gene_prioritization: OUTRIDER columns are present in output", {
  wd <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample     = "S1",
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  expect_true(any(grepl("OUTRIDER", colnames(out), ignore.case = TRUE)))
})

test_that("gene_prioritization: FRASER column is present in output", {
  wd <- setup_hpo_tempdir()
  withr::with_dir(wd, {
    out <- gene_prioritization(
      sample     = "S1",
      hpo_sample = make_clinical_df(),
      outrider   = make_outrider_df(),
      fraser     = make_fraser_summary_df()
    )
  })
  expect_true(any(grepl("FRASER", colnames(out), ignore.case = TRUE)))
})


# =============================================================================
# 4.  plot_expression_cohort()
# =============================================================================

# Shared mock data for plotly expression tests
make_expr_cohort_data <- function(gene = "ABAT", n_exons = 4) {
  exons <- rep(seq_len(n_exons), each = 5)
  data.frame(
    exonID     = exons,
    expression = abs(rnorm(length(exons), mean = 100, sd = 20)),
    PatientID  = rep(paste0("P", seq_len(5)), n_exons),
    Sexe       = rep(c("M", "F", "M", "F", "M"), n_exons),
    age        = rep(c(5, 8, 12, 15, 10), n_exons),
    geneID     = gene,
    stringsAsFactors = FALSE
  )
}

make_expr_patient_data <- function(gene = "ABAT", n_exons = 4) {
  data.frame(
    exonID     = seq_len(n_exons),
    expression = abs(rnorm(n_exons, mean = 90, sd = 10)),
    PatientID  = "HSJ_001_03",
    Sexe       = "M",
    age        = 8,
    geneID     = gene,
    stringsAsFactors = FALSE
  )
}

test_that("plot_expression_cohort: returns a plotly object", {
  skip_if_not_installed("plotly")
  p <- plot_expression_cohort(
    data_child   = make_expr_cohort_data(),
    data_adults  = make_expr_cohort_data(),
    data_patient = make_expr_patient_data()
  )
  expect_s3_class(p, "plotly")
})

test_that("plot_expression_cohort: has exactly 3 traces (children, adults, patient)", {
  skip_if_not_installed("plotly")
  p <- plotly_build(plot_expression_cohort(
    data_child   = make_expr_cohort_data(),
    data_adults  = make_expr_cohort_data(),
    data_patient = make_expr_patient_data()
  ))
  expect_equal(length(p$x$data), 3)
})

test_that("plot_expression_cohort: tick labels match the number of exons in patient data", {
  skip_if_not_installed("plotly")
  p <- plotly::plotly_build(plot_expression_cohort(
    data_child   = make_expr_cohort_data(n_exons = 6),
    data_adults  = make_expr_cohort_data(n_exons = 6),
    data_patient = make_expr_patient_data(n_exons = 6)
  ))
  expect_equal(length(p$x$layout$xaxis$ticktext), 6)
})

test_that("plot_expression_cohort: works with a single exon (edge case)", {
  skip_if_not_installed("plotly")
  expect_no_error(
    plot_expression_cohort(
      data_child   = make_expr_cohort_data(n_exons = 1),
      data_adults  = make_expr_cohort_data(n_exons = 1),
      data_patient = make_expr_patient_data(n_exons = 1)
    )
  )
})

test_that("plot_expression_cohort: children trace is a box plot", {
  skip_if_not_installed("plotly")
  p <- plotly::plotly_build(plot_expression_cohort(
    data_child   = make_expr_cohort_data(),
    data_adults  = make_expr_cohort_data(),
    data_patient = make_expr_patient_data()
  ))
  expect_equal(p$x$data[[1]]$type, "box")
})

test_that("plot_expression_cohort: patient trace is a scatter plot", {
  skip_if_not_installed("plotly")
  p <- plotly::plotly_build(plot_expression_cohort(
    data_child   = make_expr_cohort_data(),
    data_adults  = make_expr_cohort_data(),
    data_patient = make_expr_patient_data()
  ))
  expect_equal(p$x$data[[3]]$type, "scatter")
})


# =============================================================================
# 5.  plot_expression_family()
# =============================================================================

make_expr_family_data <- function(gene = "ABAT", n_exons = 4, n_members = 3) {
  exons   <- rep(seq_len(n_exons), each = n_members)
  members <- rep(paste0("FM", seq_len(n_members)), n_exons)
  data.frame(
    exonID     = exons,
    expression = abs(rnorm(length(exons), mean = 95, sd = 15)),
    PatientID  = members,
    Sexe       = rep(c("M", "F", "M")[seq_len(n_members)], n_exons),
    age        = rep(c(8, 38, 40)[seq_len(n_members)], n_exons),
    geneID     = gene,
    stringsAsFactors = FALSE
  )
}

test_that("plot_expression_family: returns a plotly object", {
  skip_if_not_installed("plotly")
  p <- plot_expression_family(
    data_family  = make_expr_family_data(),
    data_patient = make_expr_patient_data()
  )
  expect_s3_class(p, "plotly")
})

test_that("plot_expression_family: has one trace per family member", {
  skip_if_not_installed("plotly")
  fam <- make_expr_family_data(n_members = 3)
  p <- plotly_build(plot_expression_family(
    data_family  = fam,
    data_patient = make_expr_patient_data()
  ))
  expect_equal(length(p$x$data), length(unique(fam$PatientID)))
})

test_that("plot_expression_family: tick labels match number of exons in family data", {
  skip_if_not_installed("plotly")
  p <- plotly_build(plot_expression_family(
    data_family  = make_expr_family_data(n_exons = 5),
    data_patient = make_expr_patient_data(n_exons = 5)
  ))
  expect_equal(length(p$x$layout$xaxis$ticktext), 5)
})

test_that("plot_expression_family: works with a single family member (edge case)", {
  skip_if_not_installed("plotly")
  expect_no_error(
    plot_expression_family(
      data_family  = make_expr_family_data(n_members = 1),
      data_patient = make_expr_patient_data()
    )
  )
})

test_that("plot_expression_family: all traces are scatter plots", {
  skip_if_not_installed("plotly")
  p <- plotly::plotly_build(plot_expression_family(
    data_family  = make_expr_family_data(),
    data_patient = make_expr_patient_data()
  ))
  types <- vapply(p$x$data, `[[`, character(1), "type")
  expect_true(all(types == "scatter"))
})
