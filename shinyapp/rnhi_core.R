# RNHI core functions
#
# Residual normal hematopoiesis integrity (RNHI) of a sample is the negative
# Aitchison distance between the centered log-ratio (CLR) transform of its
# nine-lineage residual normal composition x and the healthy center c:
#
#   RNHI(x) = -|| clr(x) - c ||_2
#
# Zeros are handled with the multiplicative replacement rule of
# zCompositions::multRepl: each zero is replaced with 0.65 times the detection
# limit of that lineage (its minimum positive proportion) and the non-zero
# parts are rescaled so that the composition still sums to 1. The healthy
# center is the mean CLR vector of the healthy reference donors, with zero
# replacement carried out within the donors.

LINEAGES <- c('MyeProg', 'Mono', 'Gran', 'Mega', 'EarlyEry', 'LateEry', 'B', 'NK', 'DC')
ZERO_REPLACEMENT_FACTOR <- 0.65

# Normal states of the joint healthy/AML reference that form the nine RNHI lineages.
REFERENCE_STATE_MAP <- c(
  Normal_MyeloidProgenitor = 'MyeProg', Normal_Monocyte = 'Mono',
  Normal_Granulocytic = 'Gran', Normal_Megakaryocytic = 'Mega',
  Normal_EarlyErythroid = 'EarlyEry', Normal_LateErythroid = 'LateEry',
  Normal_B = 'B', Normal_NK = 'NK', Normal_Dendritic = 'DC')

# ---------------------------------------------------------------------------
# Compositional operations
# ---------------------------------------------------------------------------

check_factor <- function(factor) {
  if (!is.numeric(factor) || length(factor) != 1L || !is.finite(factor) || factor <= 0 || factor >= 1)
    stop('The zero-replacement factor must be a single number between 0 and 1.')
  factor
}

lineage_matrix <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0))
    stop('Lineage abundances must be finite, non-negative numbers; blank or missing values are not accepted.')
  if (!nrow(x) || ncol(x) != length(LINEAGES))
    stop('A non-empty table with all nine lineages is required.')
  x
}

closure <- function(x) {
  x <- lineage_matrix(x)
  empty <- rowSums(x) <= 0
  if (any(empty)) {
    ids <- if (is.null(rownames(x))) as.character(which(empty)) else rownames(x)[empty]
    stop('No lineage is detected (all nine abundances are zero) in ', paste(utils::head(ids, 6L), collapse = ', '),
         if (length(ids) > 6L) sprintf(' and %d more', length(ids) - 6L) else '', '; RNHI is undefined for such samples.')
  }
  x / rowSums(x)
}

# Detection limit of each lineage: its minimum positive proportion in a set of closed compositions.
detection_limits <- function(comp) {
  comp <- lineage_matrix(comp)
  dl <- vapply(seq_len(ncol(comp)), function(j) {
    positive <- comp[comp[, j] > 0, j]
    if (length(positive)) min(positive) else NA_real_
  }, numeric(1))
  names(dl) <- LINEAGES
  dl
}

# Multiplicative zero replacement of closed compositions. If `dl` is NULL, detection
# limits are estimated from `comp` itself; otherwise the supplied limits are used.
multRepl <- function(comp, factor = ZERO_REPLACEMENT_FACTOR, dl = NULL) {
  factor <- check_factor(factor)
  comp <- lineage_matrix(comp)
  if (is.null(dl)) dl <- detection_limits(comp)
  dl <- as.numeric(dl)
  if (length(dl) != ncol(comp)) stop('Nine detection limits are required.')
  needs_dl <- colSums(comp == 0) > 0
  missing_dl <- needs_dl & (is.na(dl) | dl <= 0)
  if (any(missing_dl))
    stop('no detection limit can be estimated for ', paste(LINEAGES[missing_dl], collapse = ', '),
         ' (zero in every sample of the group).')
  out <- comp
  for (i in seq_len(nrow(comp))) {
    zero <- which(comp[i, ] == 0)
    if (!length(zero)) next
    nonzero <- which(comp[i, ] > 0)
    replacement <- factor * dl[zero]
    scale <- (sum(comp[i, ]) - sum(replacement)) / sum(comp[i, nonzero])
    if (!is.finite(scale) || scale <= 0) stop('Zero replacement leaves no positive mass for sample ', i, '.')
    out[i, zero] <- replacement
    out[i, nonzero] <- comp[i, nonzero] * scale
  }
  out
}

clr <- function(x) {
  lx <- log(x)
  sweep(lx, 1, rowMeans(lx), '-')
}

# Healthy center: mean CLR vector of the healthy donors (zero replacement within the donors).
healthy_center <- function(health, factor = ZERO_REPLACEMENT_FACTOR) {
  colMeans(clr(multRepl(closure(health), factor)))
}

with_context <- function(expr, prefix, suffix = '') {
  tryCatch(expr, error = function(e) stop(prefix, conditionMessage(e), suffix, call. = FALSE))
}

# ---------------------------------------------------------------------------
# Input tables
# ---------------------------------------------------------------------------

ID_COLUMNS <- c('sample_id', 'spot_id', 'barcode', 'orig.ident', 'sample', 'patient_id')
TEXT_COLUMNS <- c(ID_COLUMNS, 'cell_type', 'celltype', 'cell.type', 'cohort', 'gene')
KNOWN_COLUMNS <- c(LINEAGES, TEXT_COLUMNS, 'abundance', 'healthy_reference', 'umap_1', 'umap_2', 'x', 'y',
                   'imagecol', 'imagerow', 'pxl_col_in_fullres', 'pxl_row_in_fullres')

# Column names are matched to the expected names without regard to case.
standard_names <- function(nm) {
  hit <- match(tolower(nm), tolower(KNOWN_COLUMNS))
  out <- ifelse(is.na(hit), nm, KNOWN_COLUMNS[hit])
  if (anyDuplicated(out)) nm else out
}

# Lineage labels are matched to the nine lineage names without regard to case.
standard_lineages <- function(ct) {
  hit <- match(tolower(ct), tolower(LINEAGES))
  ifelse(is.na(hit), ct, LINEAGES[hit])
}

read_input <- function(path, name = basename(path)) {
  sep <- if (grepl('\\.(tsv|txt)$', name, ignore.case = TRUE)) '\t' else ','
  dat <- read.table(path, header = TRUE, sep = sep, check.names = FALSE, colClasses = 'character',
                    stringsAsFactors = FALSE, comment.char = '', quote = '"', fileEncoding = 'UTF-8-BOM')
  if (!nrow(dat)) stop('The uploaded table has no data rows.')
  nm <- trimws(names(dat))
  blank <- is.na(nm) | !nzchar(nm)
  nm[blank] <- paste0('column_', which(blank))
  names(dat) <- standard_names(nm)
  text <- names(dat) %in% TEXT_COLUMNS
  text[1] <- text[1] || !any(LINEAGES %in% names(dat))
  for (j in which(!text)) dat[[j]] <- utils::type.convert(dat[[j]], as.is = TRUE)
  dat
}

first_column <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit)) hit[[1]] else NULL
}

healthy_flags <- function(values, ids, column) {
  if (is.logical(values)) values <- as.integer(values)
  if (is.character(values)) values <- suppressWarnings(as.numeric(trimws(values)))
  if (!is.numeric(values) || anyNA(values) || any(!values %in% c(0, 1)))
    stop(column, ' must contain only 0 (query sample) or 1 (healthy reference).')
  groups <- split(values, ids)
  if (any(vapply(groups, function(v) length(unique(v)) != 1L, logical(1))))
    stop('Conflicting ', column, ' values within one sample.')
  groups
}

# Parse a wide (nine lineage columns) or long (sample_id, cell_type, abundance) composition table.
parse_table <- function(dat) {
  if (anyDuplicated(names(dat))) stop('Duplicate column names are not accepted.')
  if (!'sample_id' %in% names(dat)) {
    id <- first_column(dat, c('spot_id', 'barcode', 'sample'))
    if (is.null(id)) stop('Add a sample_id column (spot_id or barcode for Visium spots).')
    dat$sample_id <- dat[[id]]
  }
  dat$sample_id <- as.character(dat$sample_id)
  if (anyNA(dat$sample_id) || any(!nzchar(trimws(dat$sample_id)))) stop('Sample IDs must not be blank.')
  ids <- unique(dat$sample_id)
  flags <- rep(FALSE, length(ids))
  if ('healthy_reference' %in% names(dat)) {
    groups <- healthy_flags(dat$healthy_reference, dat$sample_id, 'healthy_reference')
    flags <- vapply(ids, function(id) groups[[id]][1] == 1, logical(1), USE.NAMES = FALSE)
  }
  if (all(LINEAGES %in% names(dat))) {
    if (anyDuplicated(dat$sample_id)) stop('A wide table requires one row per sample_id.')
    x <- as.matrix(dat[, LINEAGES, drop = FALSE])
    rownames(x) <- dat$sample_id
  } else if (all(c('cell_type', 'abundance') %in% names(dat))) {
    dat$cell_type <- standard_lineages(trimws(as.character(dat$cell_type)))
    if (anyNA(dat$cell_type) || any(!dat$cell_type %in% LINEAGES))
      stop(unknown_labels(dat$cell_type))
    if (!is.numeric(dat$abundance) || any(!is.finite(dat$abundance)) || any(dat$abundance < 0))
      stop('abundance must contain finite, non-negative numbers.')
    x <- matrix(0, length(ids), length(LINEAGES), dimnames = list(ids, LINEAGES))
    z <- aggregate(dat$abundance, list(sample_id = dat$sample_id, cell_type = dat$cell_type), sum)
    x[cbind(match(z$sample_id, ids), match(z$cell_type, LINEAGES))] <- z$x
  } else stop('Provide the nine lineage columns (', paste(LINEAGES, collapse = ', '),
              '), or the long columns sample_id, cell_type and abundance.')
  list(ids = ids, x = lineage_matrix(x), healthy = flags)
}

SC_ID_COLUMNS <- c('sample_id', 'orig.ident', 'sample', 'patient_id')
SPOT_ID_COLUMNS <- c('spot_id', 'barcode', 'sample_id')

unknown_labels <- function(ct) {
  bad <- unique(ct[is.na(ct) | !ct %in% LINEAGES])
  paste0('Unknown cell-type labels: ', paste(utils::head(bad, 8L), collapse = ', '),
         '. RNHI uses only the nine lineages (', paste(LINEAGES, collapse = ', '),
         '); rename synonymous labels and remove other cell types before upload.')
}

# Annotated scRNA-seq or Visium metadata (one row per annotated cell or spot record).
metadata_to_composition <- function(dat, mode = c('sc', 'spatial')) {
  mode <- match.arg(mode)
  cell_col <- first_column(dat, c('celltype', 'cell_type', 'cell.type'))
  id_candidates <- if (mode == 'sc') SC_ID_COLUMNS else SPOT_ID_COLUMNS
  id_col <- first_column(dat, id_candidates)
  if (is.null(cell_col) || is.null(id_col)) {
    needed <- if (mode == 'sc') 'celltype and sample_id (or orig.ident)' else 'celltype and spot_id (or barcode)'
    stop('This input requires the columns ', needed, ', with one row per annotated cell.')
  }
  sid <- trimws(as.character(dat[[id_col]])); ct <- standard_lineages(trimws(as.character(dat[[cell_col]])))
  if (anyNA(sid) || any(!nzchar(sid))) stop(id_col, ' contains blank IDs.')
  if (anyNA(ct) || any(!nzchar(ct))) stop(cell_col, ' contains blank cell-type labels.')
  if (any(!ct %in% LINEAGES)) stop(unknown_labels(ct))
  d <- data.frame(sample_id = sid, cell_type = ct, abundance = 1, stringsAsFactors = FALSE)
  out <- aggregate(abundance ~ sample_id + cell_type, d, sum)
  hcol <- first_column(dat, c('healthy_reference', 'healthy', 'is_healthy'))
  if (!is.null(hcol)) {
    groups <- healthy_flags(dat[[hcol]], sid, hcol)
    hf <- data.frame(sample_id = names(groups), healthy_reference = vapply(groups, function(v) as.integer(v[1]), integer(1)),
                     stringsAsFactors = FALSE)
    out <- merge(out, hf, by = 'sample_id', all.x = TRUE, sort = FALSE)
  }
  out
}

prepare_modality_table <- function(dat, mode = c('sc', 'spatial')) {
  mode <- match.arg(mode)
  if (mode == 'spatial') {
    id <- first_column(dat, c('spot_id', 'barcode'))
    if (!is.null(id)) dat$sample_id <- dat[[id]]
  }
  if (all(LINEAGES %in% names(dat)) || all(c('cell_type', 'abundance') %in% names(dat))) return(dat)
  metadata_to_composition(dat, mode)
}

# ---------------------------------------------------------------------------
# Bulk expression input: MuSiC estimator against the joint healthy/AML reference
# ---------------------------------------------------------------------------

expression_gene_column <- function(dat) {
  hits <- intersect(c('gene', 'genes', 'gene_id', 'geneid', 'gene_symbol', 'symbol', 'feature', 'features'), names(dat))
  if (length(hits)) hits[[1]] else names(dat)[[1]]
}

is_expression_matrix <- function(dat) {
  if (any(LINEAGES %in% names(dat)) || all(c('cell_type', 'abundance') %in% names(dat))) return(FALSE)
  if (ncol(dat) < 2 || nrow(dat) < 20) return(FALSE)
  cols <- which(names(dat) != expression_gene_column(dat))
  ok <- vapply(cols, function(j) { z <- suppressWarnings(as.numeric(dat[[j]])); any(is.finite(z)) }, logical(1))
  length(cols) >= 1 && all(ok)
}

parse_expression_matrix <- function(dat) {
  gcol <- expression_gene_column(dat)
  cols <- which(names(dat) != gcol); sample_cols <- names(dat)[cols]
  genes <- trimws(as.character(dat[[gcol]]))
  if (anyNA(genes) || any(!nzchar(genes))) stop('The gene column of the expression matrix contains blank IDs.')
  if (anyDuplicated(sample_cols)) stop('Sample columns of the expression matrix must have unique names.')
  mat <- vapply(cols, function(j) suppressWarnings(as.numeric(dat[[j]])), numeric(length(genes)))
  mat <- matrix(mat, nrow = length(genes), dimnames = list(NULL, sample_cols))
  if (any(!is.finite(mat)) || any(mat < 0))
    stop('Expression values must be finite and non-negative (linear-scale values such as counts or TPM).')
  if (max(mat) > 1 && max(mat) < 50)
    stop('Expression values appear to be log-transformed (maximum ', signif(max(mat), 3),
         '); upload linear-scale values such as counts or TPM.')
  if (anyDuplicated(genes)) { mat <- rowsum(mat, group = genes, reorder = FALSE); genes <- rownames(mat) }
  else rownames(mat) <- genes
  list(matrix = mat, sample_ids = colnames(mat), gene_column = gcol)
}

# MuSiC estimator (Wang et al., Nat Commun 2019; default settings): iteratively
# re-weighted non-negative least squares with gene weights
# 1 / (nu + r_g^2 + sum_k (q_k S_k)^2 Sigma_gk), where r are the residuals of the
# previous fit, q the previous coefficients, S the average library sizes and Sigma
# the cross-subject variances of the reference. y is the relative abundance of the
# bulk sample over the reference genes (scaled by 100), restricted to expressed genes.
music_wnnls <- function(y, D, S, Sigma, nu = 1e-4, eps = 0.01, max_iter = 1000L) {
  y <- 100 * y
  gene_weight <- function(q, r) 1 / (nu + as.vector(r)^2 + as.vector(Sigma %*% (q * S)^2))
  weighted_fit <- function(w) nnls::nnls(D * sqrt(w), y * sqrt(w))
  fit <- nnls::nnls(D, y)
  fit <- weighted_fit(gene_weight(fit$x, fit$residuals))
  p <- fit$x / sum(fit$x)
  converged <- NA_integer_
  for (iter in seq_len(max_iter)) {
    new_fit <- weighted_fit(gene_weight(fit$x, fit$residuals))
    p_new <- new_fit$x / sum(new_fit$x)
    fit <- new_fit
    if (sum(abs(p_new - p)) < eps) { p <- p_new; converged <- iter; break }
    p <- p_new
  }
  fitted <- as.vector(D %*% fit$x)
  list(p = setNames(p, colnames(D)),
       r_squared = 1 - stats::var(y - fitted) / stats::var(y),
       iterations = if (is.na(converged)) max_iter else converged,
       converged = !is.na(converged))
}

deconvolve_expression <- function(dat, basis) {
  expr <- parse_expression_matrix(dat)
  common <- intersect(rownames(basis$D), rownames(expr$matrix))
  min_common <- ceiling(0.10 * nrow(basis$D))
  if (length(common) < min_common)
    stop('Too few genes match the reference (', format(length(common), big.mark = ','), ' of ',
         format(nrow(basis$D), big.mark = ','), '; at least ', format(min_common, big.mark = ','),
         ' are required). Use gene symbols in the first column.')
  D <- basis$D[common, , drop = FALSE]; Sigma <- basis$Sigma[common, , drop = FALSE]; S <- basis$S[colnames(D)]
  n <- ncol(expr$matrix)
  comp <- matrix(NA_real_, n, length(LINEAGES), dimnames = list(expr$sample_ids, LINEAGES))
  qc <- data.frame(sample_id = expr$sample_ids, reference_genes = length(common), expressed_genes = NA_integer_,
                   r_squared = NA_real_, iterations = NA_integer_, status = 'failed', stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    y <- expr$matrix[common, i]
    pos <- y > 0
    qc$expressed_genes[i] <- sum(pos)
    if (sum(pos) < ceiling(0.20 * length(common))) { qc$status[i] <- 'failed: too few expressed reference genes'; next }
    fit <- tryCatch(music_wnnls(y[pos] / sum(y), D[pos, , drop = FALSE], S, Sigma[pos, , drop = FALSE]),
                    error = function(e) e)
    if (inherits(fit, 'error') || any(!is.finite(fit$p))) { qc$status[i] <- 'failed: estimation error'; next }
    comp[i, ] <- vapply(LINEAGES, function(l) sum(fit$p[names(REFERENCE_STATE_MAP)[REFERENCE_STATE_MAP == l]]), numeric(1))
    qc$r_squared[i] <- fit$r_squared; qc$iterations[i] <- fit$iterations
    qc$status[i] <- paste0(if (fit$converged) 'converged' else 'maximum iterations reached',
                           if (fit$r_squared < 0) '; poor fit (R² < 0)' else '')
  }
  failed <- grepl('^failed', qc$status)
  if (any(failed)) {
    detail <- paste0(qc$sample_id[failed], ' (', sub('^failed: ', '', qc$status[failed]), ')')
    stop('Deconvolution failed for ', sum(failed), ' sample(s): ', paste(utils::head(detail, 6L), collapse = '; '),
         if (sum(failed) > 6L) paste0('; and ', sum(failed) - 6L, ' more') else '', '. No result was returned.')
  }
  data <- data.frame(sample_id = expr$sample_ids, comp, healthy_reference = 0L,
                     stringsAsFactors = FALSE, check.names = FALSE)
  rownames(data) <- NULL
  list(data = data, qc = qc, reference_genes = length(common), total_reference_genes = nrow(basis$D))
}

# ---------------------------------------------------------------------------
# RNHI calculation
# ---------------------------------------------------------------------------

# query_dl: detection limits applied to query samples (e.g. the fixed limits of the bulk AML
# compendium); NULL estimates them from the uploaded query samples. Healthy reference
# samples always use limits estimated within the healthy reference set.
calculate_rnhi <- function(dat, default_health, custom_health = NULL, factor = ZERO_REPLACEMENT_FACTOR, query_dl = NULL) {
  factor <- check_factor(factor)
  q <- parse_table(dat)
  composition <- closure(q$x)
  if (!is.null(custom_health) && any(q$healthy))
    stop('Use either rows flagged healthy_reference = 1 or a separate healthy reference table, not both.')
  health <- if (any(q$healthy)) q$x[q$healthy, , drop = FALSE]
            else if (!is.null(custom_health)) parse_table(custom_health)$x
            else default_health
  health <- lineage_matrix(health)
  healthy_msg <- ' Every lineage must be detected in at least one healthy reference sample.'
  center <- with_context(healthy_center(health, factor), 'Healthy reference: ', healthy_msg)
  replaced <- composition
  if (any(q$healthy))
    replaced[q$healthy, ] <- with_context(multRepl(composition[q$healthy, , drop = FALSE], factor), 'Healthy reference: ', healthy_msg)
  if (any(!q$healthy))
    replaced[!q$healthy, ] <- with_context(multRepl(composition[!q$healthy, , drop = FALSE], factor, dl = query_dl), 'Query samples: ',
      ' Use the fixed detection limits, or include samples in which this lineage is detected.')
  deviation <- sweep(clr(replaced), 2, center, '-')
  distance <- sqrt(rowSums(deviation^2))
  detected <- as.integer(rowSums(q$x > 0))
  out <- data.frame(sample_id = q$ids, healthy_reference = as.integer(q$healthy), RNHI = -distance,
                    aitchison_distance = distance, detected_lineages = detected,
                    lineage_total = rowSums(q$x),
                    qc_flag = ifelse(detected == 1, 'single lineage detected', 'ok'),
                    stringsAsFactors = FALSE)
  prop <- composition; colnames(prop) <- paste0('prop_', LINEAGES)
  dev <- deviation; colnames(dev) <- paste0('clr_diff_', LINEAGES)
  out <- cbind(out, as.data.frame(prop, row.names = NULL), as.data.frame(dev, row.names = NULL))
  rownames(out) <- NULL
  list(data = out, composition = composition, deviation = deviation, center = center,
       n = nrow(out), nhealth = nrow(health), factor = factor,
       dl_source = if (is.null(query_dl)) 'estimated' else 'fixed')
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

distribution_plot <- function(r) {
  d <- r$data
  d$population <- factor(ifelse(d$healthy_reference == 1, 'Healthy reference', 'Query samples'),
                         levels = c('Healthy reference', 'Query samples'))
  cols <- c('Healthy reference' = '#008574', 'Query samples' = '#2454a6')
  p <- ggplot2::ggplot(d, ggplot2::aes(RNHI, fill = population, color = population))
  dens_ok <- tapply(d$RNHI, d$population, function(v) length(v) > 1 && length(unique(v)) > 1)
  dens_ok[is.na(dens_ok)] <- FALSE
  if (any(dens_ok))
    p <- p + ggplot2::geom_density(data = d[d$population %in% names(dens_ok)[dens_ok], , drop = FALSE],
                                   alpha = .18, linewidth = .9)
  p + ggplot2::geom_rug(alpha = .35) +
    ggplot2::scale_color_manual(values = cols, drop = TRUE) + ggplot2::scale_fill_manual(values = cols, drop = TRUE) +
    ggplot2::theme_classic(base_size = 13) + ggplot2::theme(legend.position = 'top') +
    ggplot2::labs(x = 'RNHI (negative Aitchison distance to the healthy center)', y = 'Density', fill = NULL, color = NULL)
}

sample_heatmap <- function(r, values, fill_scale, fill_label) {
  ord <- order(r$data$RNHI, decreasing = TRUE)
  x <- values[ord, , drop = FALSE]
  z <- data.frame(rank = rep(seq_len(nrow(x)), ncol(x)),
                  lineage = factor(rep(LINEAGES, each = nrow(x)), levels = rev(LINEAGES)),
                  value = as.vector(x))
  ggplot2::ggplot(z, ggplot2::aes(rank, lineage, fill = value)) + ggplot2::geom_raster() + fill_scale +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(), legend.position = 'right') +
    ggplot2::labs(x = 'Samples ordered by RNHI (high to low)', y = NULL, fill = fill_label)
}

composition_plot <- function(r) {
  sample_heatmap(r, r$composition,
    ggplot2::scale_fill_gradientn(colours = c('#F7FBFF', '#DEEBF7', '#9ECAE1', '#4292C6', '#08519C'),
                                  limits = c(0, 1), oob = scales::squish, labels = scales::percent),
    'Proportion\nwithin sample')
}

profile_plot <- function(r) {
  sample_heatmap(r, r$deviation,
    ggplot2::scale_fill_gradient2(low = '#2166AC', mid = '#F7F7F7', high = '#B2182B', midpoint = 0,
                                  limits = c(-4, 4), oob = scales::squish),
    'CLR difference\nfrom healthy\ncenter')
}

# Cell-level view for scRNA-seq metadata that include UMAP coordinates (umap_1, umap_2).
sc_umap_data <- function(r, dat) {
  if (is.null(dat)) return(NULL)
  id_col <- first_column(dat, SC_ID_COLUMNS)
  u1 <- first_column(dat, 'umap_1'); u2 <- first_column(dat, 'umap_2')
  if (is.null(id_col) || is.null(u1) || is.null(u2)) return(NULL)
  z <- data.frame(sample_id = trimws(as.character(dat[[id_col]])),
                  umap_1 = suppressWarnings(as.numeric(dat[[u1]])),
                  umap_2 = suppressWarnings(as.numeric(dat[[u2]])), stringsAsFactors = FALSE)
  if (any(!is.finite(z$umap_1)) || any(!is.finite(z$umap_2)))
    stop('UMAP coordinates (umap_1, umap_2) must be finite numbers.')
  merge(z, r$data[, c('sample_id', 'RNHI')], by = 'sample_id', all.x = TRUE, sort = FALSE)
}

RNHI_MAP_COLORS <- c('#A8C7DD', '#C8DCEA', '#F0E0C2', '#F2B18F', '#D36B67')

sc_umap_plot <- function(r) {
  ggplot2::ggplot(r$plot_data, ggplot2::aes(umap_1, umap_2, color = RNHI)) +
    ggplot2::geom_point(size = .9, alpha = .7) +
    ggplot2::scale_color_gradientn(colours = RNHI_MAP_COLORS, name = 'Sample RNHI') +
    ggplot2::coord_equal() + ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Uploaded cells colored by the RNHI of their sample')
}

# Spot-level map for Visium data with spot coordinates (x/y, imagecol/imagerow or
# pxl_col_in_fullres/pxl_row_in_fullres). Image coordinates are drawn with the y axis reversed.
spatial_data <- function(r, dat) {
  if (is.null(dat)) return(NULL)
  id_col <- first_column(dat, SPOT_ID_COLUMNS)
  pairs <- list(c('x', 'y', 'cartesian'), c('imagecol', 'imagerow', 'image'),
                c('pxl_col_in_fullres', 'pxl_row_in_fullres', 'image'))
  hit <- Filter(function(p) all(p[1:2] %in% names(dat)), pairs)
  if (is.null(id_col) || !length(hit)) return(NULL)
  p <- hit[[1]]
  z <- data.frame(sample_id = trimws(as.character(dat[[id_col]])),
                  x = suppressWarnings(as.numeric(dat[[p[1]]])),
                  y = suppressWarnings(as.numeric(dat[[p[2]]])), stringsAsFactors = FALSE)
  z <- z[!duplicated(z$sample_id), , drop = FALSE]
  if (any(!is.finite(z$x)) || any(!is.finite(z$y)))
    stop('Spot coordinates must be finite numbers.')
  z$coord_type <- p[3]
  merge(z, r$data[, c('sample_id', 'RNHI')], by = 'sample_id', all.x = TRUE, sort = FALSE)
}

spatial_plot <- function(r) {
  z <- r$plot_data
  rng <- range(z$RNHI, finite = TRUE)
  if (!all(is.finite(rng))) stop('The spatial map has no finite RNHI values.')
  rng <- c(rng[1], max(0, rng[2]))
  if (diff(rng) == 0) rng <- rng + c(-.5, .5)
  span <- max(diff(range(z$x)), diff(range(z$y)))
  spot_size <- if (is.finite(span) && span <= 20) 8.2 else 5.4
  y_scale <- if (any(z$coord_type == 'image')) ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = .13))
             else ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = .13))
  ggplot2::ggplot(z, ggplot2::aes(x, y, fill = RNHI)) +
    ggplot2::geom_point(shape = 21, size = spot_size, alpha = .96, colour = '#fffdf8', stroke = .32) +
    ggplot2::scale_fill_gradientn(colours = RNHI_MAP_COLORS, limits = rng, oob = scales::squish, name = 'RNHI\n(0 = healthy-like)',
                                  guide = ggplot2::guide_colorbar(barheight = grid::unit(55, 'mm'), barwidth = grid::unit(5, 'mm'))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .13)) + y_scale +
    ggplot2::coord_equal() + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = 'Spot-level RNHI', subtitle = 'Higher (closer to 0) indicates a composition closer to healthy bone marrow',
                  x = NULL, y = NULL) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = .5, face = 'bold'),
                   plot.subtitle = ggplot2::element_text(hjust = .5, size = 10, colour = '#5b6862'))
}
