script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_args)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_args[1]), mustWork = TRUE))
  normalizePath(file.path(script_dir, "../../.."), mustWork = TRUE)
} else normalizePath(Sys.getenv("PROJECT_ROOT", "."), mustWork = TRUE)
setwd(project_root)
anonymous_ids <- function(x) {
  x <- as.character(x)
  if (anyNA(x) || any(!grepl("^[1-9][0-9]*$", x))) stop("Anonymous positive-integer identifiers are required.")
  x
}

suppressMessages({ library(ComplexHeatmap); library(circlize) })
args <- commandArgs(trailingOnly = FALSE)
src <- "data/figure_1/marker_heatmap"
a <- commandArgs(trailingOnly = TRUE)
if (length(a) >= 1 && nzchar(a[1])) src <- a[1]
outd <- if (length(a) >= 2 && nzchar(a[2])) a[2] else "output/figure_1/marker_heatmap"
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

inp <- src
matf <- file.path(inp, "heatmap_matrix_v4.csv"); metaf <- file.path(inp, "gene_meta_v4.csv"); annf <- file.path(inp, "col_annotation.csv")
stopifnot(file.exists(matf), file.exists(metaf), file.exists(annf))
M <- read.csv(matf, check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
rownames(M) <- M$feature; M$feature <- NULL
meta <- read.csv(metaf, check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
ann <- read.csv(annf, check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
ann <- ann[match(colnames(M), ann$sample), , drop = FALSE]
stopifnot(all(ann$sample == colnames(M)))

STATE_ORDER <- colnames(M)
meta <- meta[match(rownames(M), meta$feature), , drop = FALSE]
stopifnot(all(meta$feature == rownames(M)))
Z_full <- t(scale(t(as.matrix(M)))); Z_full[!is.finite(Z_full)] <- 0

drop_cols <- trimws(strsplit(Sys.getenv("DROP_COLS", ""), ",")[[1]])
drop_cols <- drop_cols[nzchar(drop_cols)]
keep <- setdiff(STATE_ORDER, drop_cols)
stopifnot(length(keep) >= 2, all(drop_cols %in% STATE_ORDER))
if (length(drop_cols)) cat("dropping columns:", paste(drop_cols, collapse=","), "\n")
M <- M[, keep, drop = FALSE]; ann <- ann[ann$sample %in% keep, , drop = FALSE]
ann <- ann[match(colnames(M), ann$sample), , drop = FALSE]
STATE_ORDER <- colnames(M)
Z <- Z_full[, keep, drop = FALSE]

block_levels <- unique(meta$block)
block_levels <- c(intersect(STATE_ORDER, block_levels), setdiff(block_levels, STATE_ORDER))
ord <- integer(0)
for (blk in block_levels) {
  i <- which(meta$block == blk)
  if (!length(i)) next
  j <- match(blk, STATE_ORDER)
  i <- i[order(-Z[i, j])]
  ord <- c(ord, i)
}
Z <- Z[ord, , drop = FALSE]; meta <- meta[ord, , drop = FALSE]
row_split <- factor(meta$block, levels = block_levels)
lab_idx <- which(meta$label == 1)
lab_names <- meta$feature[lab_idx]
cat(sprintf("rows=%d blocks=%d labelled=%d\n", nrow(Z), length(block_levels), length(lab_idx)))

FS_GENE  <- as.numeric(Sys.getenv("FS_GENE",  "13"))

X_ROT    <- as.numeric(Sys.getenv("X_ROT",    "45"))
PAD_MM   <- as.numeric(Sys.getenv("PAD_MM",   "1.2"))
LINK_W   <- as.numeric(Sys.getenv("LINK_W",   "4"))
LINK_H   <- as.numeric(Sys.getenv("LINK_H",   "2.5"))
ROW_GAP  <- as.numeric(Sys.getenv("ROW_GAP",  "1.2"))
WIDTH_IN <- as.numeric(Sys.getenv("WIDTH_IN", "11.5"))
H_IN     <- suppressWarnings(as.numeric(Sys.getenv("H_IN", "")))

cols2 <- c(Normal = "#0072B2", Leukemia = "#D55E00")
col_fun <- colorRamp2(c(min(Z), 0, max(Z)), c("#2166AC", "white", "#B2182B"))
ha_top <- HeatmapAnnotation(group = factor(ann$group, levels = c("Normal", "Leukemia")),
                            col = list(group = cols2), annotation_name_gp = gpar(fontsize = 8),
                            simple_anno_size = unit(0.28, "cm"))
ha_grp <- rowAnnotation(`Gene group` = factor(ifelse(grepl("^AML_", meta$block), "Leukemia", "Normal"),
                                              levels = c("Normal", "Leukemia")),
                        col = list(`Gene group` = cols2), annotation_name_gp = gpar(fontsize = 8),
                        annotation_name_side = "top", simple_anno_size = unit(0.25, "cm"))
ha_mark <- rowAnnotation(` ` = anno_mark(at = lab_idx, labels = lab_names, which = "row", side = "left",
                                         labels_gp = gpar(fontsize = FS_GENE),
                                         link_width = unit(LINK_W, "mm"), link_height = unit(LINK_H, "mm"),
                                         link_gp = gpar(lwd = 0.6, col = "grey20"), padding = unit(PAD_MM, "mm")))
nrow_ <- nrow(Z)
ht <- Heatmap(Z, name = "row z(MTheta)", col = col_fun,
              cluster_rows = FALSE, row_order = seq_len(nrow_), row_split = row_split,
              cluster_row_slices = FALSE, row_title = NULL,
              top_annotation = ha_top, left_annotation = c(ha_mark, ha_grp),
              cluster_columns = FALSE, show_row_names = FALSE,
              show_column_names = TRUE, column_names_gp = gpar(fontsize = 8),
              column_names_rot = X_ROT, column_names_centered = FALSE,
              column_names_max_height = unit(4.2, "cm"),
              row_gap = unit(ROW_GAP, "mm"),
              heatmap_legend_param = list(title = "row z-score", title_gp = gpar(fontsize = 7),
                                          labels_gp = gpar(fontsize = 6)))

need_lab <- (length(lab_idx) * 55 + 700) / 300
height_in <- if (!is.na(H_IN) && H_IN > 0) H_IN else max(9, min(28, max(nrow_ * 0.07, need_lab)))
for (dev_type in c("png", "pdf")) {
  tag <- if (length(drop_cols)) paste0("_v4b_drop_", paste(drop_cols, collapse="-"), "columns") else "_v4_canonical_markers"
  f <- file.path(outd, paste0("marker_heatmap", tag, ".", dev_type))
  if (dev_type == "png") png(f, width = WIDTH_IN, height = height_in, units = "in", res = 300, bg = "white") else cairo_pdf(f, width = WIDTH_IN, height = height_in)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = "Joint-reference canonical-marker heatmap (Normal vs Leukemia)",
       column_title_gp = gpar(fontsize = 10))
  dev.off()
}
write.csv(data.frame(feature = meta$feature, block = meta$block, is_marker = meta$is_marker,
                     forced = meta$forced, label = meta$label, row_order = seq_len(nrow_)),
          file.path(outd, "marker_row_order_v4.csv"), row.names = FALSE)
cat("Output completed.\n")
