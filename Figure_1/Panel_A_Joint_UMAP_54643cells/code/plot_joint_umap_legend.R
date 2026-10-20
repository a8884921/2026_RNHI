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

suppressMessages(library(grDevices))
args <- commandArgs(trailingOnly = TRUE)
OUT <- "output/figure_1/joint_umap"
CSV <- "data/figure_1/joint_umap/harmony_umap_table.csv"
if (length(args) >= 1 && file.exists(args[1])) CSV <- args[1]
if (length(args) >= 2) OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(CSV, check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
stopifnot(all(c("UMAP_1", "UMAP_2", "reference_state", "batch") %in% names(d)))
order_lv <- c("Normal_HSC_MPP","Normal_LMPP","Normal_MyeloidProgenitor","Normal_Monocyte","Normal_Granulocytic","Normal_Megakaryocytic",
              "Normal_EarlyErythroid","Normal_LateErythroid","Normal_B","Normal_T","Normal_NK","Normal_Dendritic",
              "AML_LSPC_Quiescent","AML_LSPC_Primed","AML_LSPC_Cycling","AML_GMP_like","AML_ProMono_like","AML_Mono_like","AML_cDC_like")
lv <- intersect(order_lv, unique(d$reference_state))
lv <- c(lv, setdiff(unique(d$reference_state), lv))
d$state <- factor(d$reference_state, levels = lv)

normal_col <- c("#2166AC","#4393C3","#92C5DE","#4A7BBB","#6BAED6","#9ECAE1","#2C7FB8","#74ADD1","#6A51A3","#807DBA","#9E9AC8","#BCBDDC")
leuk_col   <- c("#D55E00","#E08214","#FDB863","#CB181D","#FB6A4A","#FC9272","#FF7F00")
cols <- c(setNames(normal_col, order_lv[grepl("^Normal_", order_lv)]),
          setNames(leuk_col, order_lv[grepl("^AML_", order_lv)]))
cols <- cols[lv]

grp_cells <- table(d$state)
grp_samp  <- vapply(lv, function(z) length(unique(d$batch[d$state == z])), integer(1))
lab <- sprintf("%s  (n=%s cells, %s samples)", lv, formatC(grp_cells[lv], big.mark=","), grp_samp[lv])
total_n <- nrow(d); total_s <- length(unique(d$batch))
colv <- cols[as.character(d$state)]
colpt <- grDevices::adjustcolor(colv, alpha.f = 0.55)

png(file.path(OUT, "Joint_UMAP_Harmony_v3_bottom_legend.png"), width = 5, height = 5, units="in", res=300, bg="white")
nf <- layout(matrix(c(1,2), 2, 1), heights = c(1.25, 1))
par(mar = c(2.5, 3, 2.2, 1), cex.axis = 0.75, cex.lab = 0.8)
plot(d$UMAP_1, d$UMAP_2, col = colpt, pch = 16, cex = 0.35, axes = TRUE,
     xlab = "Harmony UMAP 1", ylab = "Harmony UMAP 2", bty = "l")
title(main = sprintf("Joint reference UMAP (Harmony, n_total=%s cells, %s samples)", formatC(total_n, big.mark=","), total_s),
      cex.main = 1.0, line = 0.6, adj = 0)
par(mar = c(0.2, 1, 0.2, 1))
plot.new()
legend(x = "center", legend = lab, col = cols, pch = 16, bty = "n", ncol = 4, cex = 0.62, x.intersp = 0.8, y.intersp = 0.9)
dev.off()
cairo_pdf(file.path(OUT, "Joint_UMAP_Harmony_v3_bottom_legend.pdf"), width = 11, height = 13)
nf <- layout(matrix(c(1,2), 2, 1), heights = c(1.25, 1))
par(mar = c(2.5, 3, 2.2, 1), cex.axis = 0.75, cex.lab = 0.8)
plot(d$UMAP_1, d$UMAP_2, col = colpt, pch = 16, cex = 0.35, axes = TRUE,
     xlab = "Harmony UMAP 1", ylab = "Harmony UMAP 2", bty = "l")
title(main = sprintf("Joint reference UMAP (Harmony, n_total=%s cells, %s samples)", formatC(total_n, big.mark=","), total_s),
      cex.main = 1.0, line = 0.6, adj = 0)
par(mar = c(0.2, 1, 0.2, 1))
plot.new()
legend(x = "center", legend = lab, col = cols, pch = 16, bty = "n", ncol = 4, cex = 0.62, x.intersp = 0.8, y.intersp = 0.9)
dev.off()
cat("Output completed.\n")
