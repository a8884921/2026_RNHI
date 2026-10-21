# Regression tests for the RNHI calculator.
# Run from the application folder:  Rscript tests/test_rnhi.R
#
# tests/reference_compendium_rnhi.csv: RNHI of the 3,310 AML transcriptomes of the bulk compendium with a
#   detected lineage, computed once on the complete compendium (all samples sharing the same detection limits).
# tests/reference_expression_proportions.csv: cell-state proportions of the simulated expression
#   example estimated with MuSiC::music.iter (MuSiC 1.0.0, default settings).

args <- commandArgs(FALSE)
file_arg <- sub('^--file=', '', args[grep('^--file=', args)])
if (length(file_arg)) setwd(normalizePath(file.path(dirname(file_arg[[1]]), '..')))
suppressPackageStartupMessages(library(nnls))
source('rnhi_core.R', local = TRUE)

check <- function(ok, label) { if (!isTRUE(ok)) stop('FAILED: ', label, call. = FALSE); cat('PASS ', label, '\n', sep = '') }
fails_with <- function(expr, pattern) { m <- tryCatch({ expr; '' }, error = function(e) conditionMessage(e)); grepl(pattern, m) }
write_tmp <- function(x, ext = '.csv', ...) { f <- tempfile(fileext = ext); utils::write.table(x, f, sep = if (ext == '.tsv') '\t' else ',', quote = FALSE, ...); f }
tol <- 1e-12

compendium <- read_input('data/bulk_compendium_composition.csv')
parsed <- parse_table(compendium)
health <- parsed$x[parsed$healthy, , drop = FALSE]
dl <- detection_limits(closure(parsed$x[!parsed$healthy, , drop = FALSE]))
reference <- read.csv('tests/reference_compendium_rnhi.csv', stringsAsFactors = FALSE)
check(nrow(health) == 45L && nrow(reference) == 3310L, 'reference data: 45 healthy donors and 3,310 AML transcriptomes')

# 1-3. Bulk compendium: complete set and single samples
full <- calculate_rnhi(compendium, health, query_dl = dl)$data
idx <- match(reference$sample_id, full$sample_id)
check(!anyNA(idx) && max(abs(full$RNHI[idx] - reference$RNHI)) < tol, 'complete compendium reproduces the reference values')
est <- calculate_rnhi(compendium, health)$data
check(max(abs(est$RNHI[idx] - reference$RNHI)) < tol, 'estimated detection limits equal the fixed limits on the complete compendium')
aml <- compendium[compendium$healthy_reference == 0, , drop = FALSE]
single <- vapply(seq_len(nrow(aml)), function(i) calculate_rnhi(aml[i, , drop = FALSE], health, query_dl = dl)$data$RNHI, numeric(1))
check(max(abs(single - reference$RNHI[match(aml$sample_id, reference$sample_id)])) < tol,
      'single-sample scoring with fixed detection limits reproduces the reference values')

# 4. Error messages
check(fails_with(calculate_rnhi(aml[1, , drop = FALSE], health), '^Query samples: no detection limit can be estimated'),
      'estimated limits on a single sample with zeros give an informative error')
zero <- aml[1:2, ]; zero[1, LINEAGES] <- 0
check(fails_with(calculate_rnhi(zero, health, query_dl = dl), zero$sample_id[1]), 'a sample with no detected lineage is named in the error')

# 5. Bulk expression input
basis <- readRDS('data/joint_reference_basis.rds')
expr <- read_input('data/simulated_expression_demo.csv')
dc <- deconvolve_expression(expr, basis)
check(nrow(dc$data) == 8L && all(is.finite(as.matrix(dc$data[, LINEAGES]))) && all(grepl('^converged', dc$qc$status)),
      'expression matrix deconvolved for all samples')
ref_p <- read.csv('tests/reference_expression_proportions.csv', check.names = FALSE)
X <- parse_expression_matrix(expr)$matrix; common <- intersect(rownames(basis$D), rownames(X))
dmax <- 0
for (i in seq_len(ncol(X))) {
  y <- X[common, i]; pos <- y > 0
  p <- music_wnnls(y[pos] / sum(y), basis$D[common, ][pos, ], basis$S[colnames(basis$D)], basis$Sigma[common, ][pos, ])$p
  dmax <- max(dmax, max(abs(p - as.numeric(ref_p[ref_p$sample_id == colnames(X)[i], colnames(basis$D)]))))
}
check(dmax < 1e-10, 'MuSiC estimator reproduces MuSiC::music.iter proportions')
one <- deconvolve_expression(expr[, 1:2], basis)
check(max(abs(as.numeric(one$data[1, LINEAGES]) - as.numeric(dc$data[1, LINEAGES]))) < tol,
      'deconvolution of one profile does not depend on the other uploaded profiles')
blank_header <- read_input(write_tmp(setNames(expr[, 1:3], c('', names(expr)[2:3])), col.names = TRUE, row.names = FALSE))
check(nrow(deconvolve_expression(blank_header, basis)$data) == 2L, 'expression matrix with a blank first header is read')
log_expr <- expr; log_expr[, -1] <- log2(as.matrix(expr[, -1]) / 1000 + 1)
check(fails_with(deconvolve_expression(log_expr, basis), 'log-transformed'), 'log-transformed expression values are rejected')
dup <- expr[, c(1, 2, 2)]; names(dup) <- c('gene', 'S1', 'S1')
check(fails_with(deconvolve_expression(dup, basis), 'unique names'), 'duplicated sample columns are rejected')

# 6. Input handling
ids <- read_input(write_tmp(data.frame(sample_id = c('007', '008'), myeprog = 1:2, mono = 1, gran = 1, mega = 1,
                                       earlyery = 1, lateery = 1, b = 1, nk = 1, dc = 1), row.names = FALSE), '.csv')
pt <- parse_table(ids)
check(identical(pt$ids, c('007', '008')) && is_expression_matrix(ids) == FALSE, 'sample IDs kept as text; lineage names matched without regard to case')
tsv <- read_input(write_tmp(compendium[1:5, ], '.tsv', row.names = FALSE), 'x.tsv')
check(nrow(parse_table(tsv)$x) == 5L, 'TSV input is read')

# 7. scRNA-seq and Visium
sc_meta <- data.frame(orig.ident = rep(c('S01', 'S02'), each = length(LINEAGES)), celltype = rep(LINEAGES, 2),
                      UMAP_1 = rep(seq_along(LINEAGES), 2), UMAP_2 = rep(rev(seq_along(LINEAGES)), 2), stringsAsFactors = FALSE)
sc_meta <- read_input(write_tmp(sc_meta, row.names = FALSE))
sc <- calculate_rnhi(prepare_modality_table(sc_meta, 'sc'), health)
sc$plot_data <- sc_umap_data(sc, sc_meta)
check(sc$n == 2L && nrow(sc$plot_data) == nrow(sc_meta) && all(is.finite(sc$plot_data$RNHI)), 'scRNA-seq metadata counted per sample, UMAP retained')
healthy_meta <- data.frame(sample_id = rep(c('H1', 'H2'), each = length(LINEAGES)), celltype = rep(LINEAGES, 2))
check(calculate_rnhi(prepare_modality_table(sc_meta, 'sc'), health, prepare_modality_table(healthy_meta, 'sc'))$nhealth == 2L,
      'healthy reference supplied as cell metadata')
sp_meta <- data.frame(spot_id = rep(sprintf('spot_%02d', 1:4), each = length(LINEAGES)), celltype = rep(LINEAGES, 4),
                      imagecol = rep(rep(1:2, each = 2), each = length(LINEAGES)), imagerow = rep(rep(1:2, 2), each = length(LINEAGES)))
sp <- calculate_rnhi(prepare_modality_table(sp_meta, 'spatial'), health)
sp$plot_data <- spatial_data(sp, sp_meta)
check(sp$n == 4L && nrow(sp$plot_data) == 4L && all(sp$plot_data$coord_type == 'image'), 'Visium metadata counted per spot, image coordinates retained')
sp_long <- read_input('data/simulated_visium_demo.csv'); sp_long$sample_id <- 'section_1'
spl <- calculate_rnhi(prepare_modality_table(sp_long, 'spatial'), health)
spl$plot_data <- spatial_data(spl, sp_long)
check(spl$n == 144L && all(is.finite(spl$plot_data$RNHI)), 'Visium table with section and spot IDs is scored per spot')
sc_counts <- read_input('data/simulated_scrnaseq_demo.csv')
check(is.null(sc_umap_data(calculate_rnhi(sc_counts, health), sc_counts)), 'no UMAP is drawn without supplied coordinates')
conflict <- transform(sc_meta, healthy_reference = c(1L, rep(0L, nrow(sc_meta) - 1L)))
check(fails_with(prepare_modality_table(conflict, 'sc'), 'Conflicting'), 'conflicting healthy_reference flags rejected')
check(fails_with(prepare_modality_table(transform(sc_meta, celltype = 'T'), 'sc'), 'remove other cell types'), 'other cell types are reported')

cat('All tests passed.\n')
