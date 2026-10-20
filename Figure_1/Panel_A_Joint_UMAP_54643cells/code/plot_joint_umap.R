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

suppressMessages({library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- "output/figure_1/joint_umap"
CSV <- "data/figure_1/joint_umap/harmony_umap_table.csv"
if (length(args) >= 1 && file.exists(args[1])) CSV <- args[1]
if (length(args) >= 2) OUT <- args[2]
d <- read.csv(CSV, check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
order_lv <- c("Normal_HSC_MPP","Normal_LMPP","Normal_MyeloidProgenitor","Normal_Monocyte","Normal_Granulocytic","Normal_Megakaryocytic",
              "Normal_EarlyErythroid","Normal_LateErythroid","Normal_B","Normal_T","Normal_NK","Normal_Dendritic",
              "AML_LSPC_Quiescent","AML_LSPC_Primed","AML_LSPC_Cycling","AML_GMP_like","AML_ProMono_like","AML_Mono_like","AML_cDC_like")
lv <- intersect(order_lv, unique(d$reference_state)); lv <- c(lv, setdiff(unique(d$reference_state), lv))
d$state <- factor(d$reference_state, levels = lv)
cols <- c(setNames(c("#2166AC","#4393C3","#92C5DE","#4A7BBB","#6BAED6","#9ECAE1","#2C7FB8","#74ADD1","#6A51A3","#807DBA","#9E9AC8","#BCBDDC"), grep("^Normal_", order_lv, value=TRUE)),
          setNames(c("#D55E00","#E08214","#FDB863","#CB181D","#FB6A4A","#FC9272","#FF7F00"), grep("^AML_", order_lv, value=TRUE)))
cols <- cols[lv]
grp_cells <- table(d$state); grp_samp <- vapply(lv, function(z) length(unique(d$batch[d$state == z])), integer(1))
short <- sub("^Normal_","",sub("^AML_","",lv))
lab <- sprintf("%s (n=%s,%s)", short, formatC(grp_cells[lv], big.mark=","), grp_samp[lv])
total_n <- nrow(d); total_s <- length(unique(d$batch))
p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = state)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = cols, labels = lab, breaks = lv, name = NULL) +
  labs(x = "Harmony UMAP 1", y = "Harmony UMAP 2",
       title = sprintf("Joint UMAP (Harmony)\nn=%s cells, %s samples", formatC(total_n, big.mark=","), total_s)) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 2)),
        axis.title = element_text(size = 13), axis.text = element_text(size = 10.5),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8.2),
        legend.key.height = unit(0.32, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
        plot.margin = margin(4, 4, 2, 4)) +
  guides(colour = guide_legend(ncol = 3, byrow = FALSE, override.aes = list(size = 2, alpha = 1)))
ggsave(file.path(OUT, "Joint_UMAP_Harmony_v5_5x5.png"), p, width = 5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Joint_UMAP_Harmony_v5_5x5.pdf"), p, width = 5, height = 5, device = cairo_pdf)
cat("done v5 5x5\n")
