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

fread <- function(...) {
  x <- data.table::fread(...)
  id_columns <- intersect(c("sample_id", "sample", "donor", "code", "cell_id"), names(x))
  for (key in id_columns) x[[key]] <- anonymous_ids(x[[key]])
  x
}

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
options(warn = -1)
W <- "output/figure_5/neighbourhood"
dir.create(W, recursive = TRUE, showWarnings = FALSE)
O <- "output/figure_5/neighbourhood"
dir.create(O, recursive = TRUE, showWarnings = FALSE)
FONT <- Sys.getenv("FIG_FONT","Arial")

E1 <- fread(file.path(W,"preserved_vs_impaired_neighborhood_by_patient.csv"))
E2 <- fread(file.path(W,"enrichment_impaired_vs_preserved_richneighbor.csv"))
E3 <- fread(file.path(W,"neighborhood_niche_vs_srnhi_by_patient.csv"))
labc <- setNames(as.character(1:4), as.character(1:4))
nlab <- c(MSC="MSC (CAR-like)",Endo="Endothelial",Adipo="Adipocyte")

D <- rbindlist(list(
  E1[, .(code, niche, stat="a", effect=diff_imp_minus_pres, p=wilcox_p)],
  E2[, .(code, niche, stat="b", effect=diff_imp_minus_pres, p=perm_p)],
  E3[, .(code, niche, stat="c", effect=spearman_rho,          p=spearman_p)]))
D[, fdr := p.adjust(p,"BH"), by=stat]
D[, sig := fdr < 0.05]
nref <- unique(E1[, .(code, n_pres, n_imp)])
D <- merge(D, nref[, .(code, np=n_pres, ni=n_imp)], by="code", all.x=TRUE)
D[, ylab := factor(sprintf("%s\nn = %d / %d", labc[code], np, ni),
                   levels = rev(sprintf("%s\nn = %d / %d", labc[nref$code], nref$n_pres, nref$n_imp)))]
D[, niche_lab := factor(nlab[niche], levels=unname(nlab[c("MSC","Endo","Adipo")]))]
STAT <- c(a="Mean neighbourhood abundance\n(impaired \u2212 preserved)",
          b="Niche-rich neighbours\n(fraction; impaired \u2212 preserved)",
          c="Spearman \u03c1\n(neighbourhood vs sRNHI)")
fwrite(D[, .(code, niche, stat, effect=round(effect,5), p=signif(p,3), fdr=signif(fdr,3), sig)],
       file.path(O,"FigE_combined_statistics.csv"))
pal <- c("MSC (CAR-like)"="#C1272D","Endothelial"="#2166AC","Adipocyte"="#E69F00")

build <- function(cols, cap, nrow_facet=TRUE) {
  d <- D[stat %in% cols]; d[, stat_lab := factor(STAT[stat], levels=unname(STAT[cols]))]
  sl <- d[sig==TRUE]
  ggplot(d, aes(effect, ylab)) +
    geom_vline(xintercept=0, linetype="22", colour="grey55", linewidth=.3) +
    geom_segment(aes(x=0, xend=effect, y=ylab, yend=ylab, colour=niche_lab), linewidth=.75) +
    geom_point(shape=21, fill="white", colour="grey45", size=2.0, stroke=.45) +
    geom_point(data=sl, aes(fill=niche_lab), shape=21, colour="grey20", size=2.05, stroke=.45) +
    geom_text(data=sl, aes(x=effect, y=ylab, label="*"), inherit.aes=FALSE,
              vjust=-0.75, size=2.9, colour="grey15", family=FONT) +
    facet_grid(niche_lab ~ stat_lab, scales="free_x", switch="y") +
    scale_x_continuous(expand=expansion(mult=0.16)) +
    scale_fill_manual(values=pal, guide="none") + scale_colour_manual(values=pal, guide="none") +
    labs(x=NULL, y=NULL, caption=paste(strwrap(cap, width=150), collapse="\n")) +
    theme_bw(base_size=8, base_family=FONT) +
    theme(strip.background=element_rect(fill="grey94", colour="grey70", linewidth=.3),
          strip.text=element_text(size=6.6, face="bold", colour="grey15"),
          strip.text.y.left=element_text(angle=0, hjust=1, size=6.8),
          strip.text.x=element_text(size=6.3, lineheight=1.05),
          axis.text.y=element_text(size=5.8, colour="grey15", lineheight=.95),
          axis.text.x=element_text(size=5.8, colour="grey35"),
          panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.18, colour="grey92"),
          panel.spacing.x=unit(3.2,"mm"), panel.spacing.y=unit(2.0,"mm"),
          plot.caption=element_text(size=5.4, colour="grey35", hjust=0, lineheight=1.15),
          plot.caption.position="plot", plot.margin=margin(4,5,3,4))
}
CAP_MAIN <- "Filled symbols: BH-FDR < 0.05."
CAP_SUP <- "Filled symbols: BH-FDR < 0.05."

g1 <- build(c("a","c"), CAP_MAIN)
grDevices::cairo_pdf(file.path(O,"FigE1_neighborhood_summary.pdf"), width=6.0, height=5.4, family=FONT)
print(g1); dev.off()
ragg::agg_png(file.path(O,"FigE1_neighborhood_summary.png"), width=6.0, height=5.4, units="in", res=400, background="white")
print(g1); dev.off()

if (Sys.getenv("FIGE_S1","0") == "1") {
  g2 <- build("b", CAP_SUP)
  grDevices::cairo_pdf(file.path(O,"FigE_S1_richneighbour_fraction.pdf"), width=3.4, height=5.0, family=FONT)
  print(g2); dev.off()
  ragg::agg_png(file.path(O,"FigE_S1_richneighbour_fraction.png"), width=3.4, height=5.0, units="in", res=400, background="white")
  print(g2); dev.off()
}

cat("[done] FigE1_neighborhood_summary (a+c, 152 x 137 mm); enable FIGE_S1=1 for panel b.\n")
cat(sprintf("[significance] BH-FDR<0.05 cells: a=%d b=%d c=%d / 12\n",
            D[stat=="a"&sig==TRUE,.N], D[stat=="b"&sig==TRUE,.N], D[stat=="c"&sig==TRUE,.N]))
