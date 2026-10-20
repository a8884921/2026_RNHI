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

suppressPackageStartupMessages({library(ggplot2); library(data.table); library(patchwork); library(ggsignif)})

res <- "output/figure_4/continuity"
cells <- fread("data/figure_4/6c_cells.tsv", sep="\t")
ct    <- fread(file.path(res, "Figure6C_patient_continuity.csv"))
tt    <- fread(file.path(res, "Figure6C_continuity_tests.csv"))
fmt2  <- function(x) formatC(x, format="e", digits=1)

grpDisp <- c("Healthy"="Healthy","AML_RNHI-high"="RNHI-high AML","AML_RNHI-low"="RNHI-low AML")
cells[, gd := factor(grpDisp[group], levels=unname(grpDisp))]
cells[, pt := suppressWarnings(as.numeric(pseudotime))]
set.seed(6)
plotc <- rbindlist(list(cells[group=="Healthy"][sample(.N, min(.N, 80000L))],
                        cells[group!="Healthy"]))
n_disp <- plotc[, .(n=.N), by=.(group, gd)]
subA <- paste(sprintf("%s n=%d", grpDisp[n_disp$group], n_disp$n), collapse="   ")
pseudo_panels <- function(glev, colname) {
  d <- plotc[gd==glev]
  d <- d[is.finite(pt)]
  ggplot(d, aes(u1, u2, colour=pt)) +
    geom_point(size=0.10, stroke=0, alpha=0.55) +
    scale_colour_viridis_c(limits=c(0,22.5), breaks=seq(0,20,5), name="pseudotime",
                           guide=guide_colourbar(barwidth=unit(3.2,"cm"), barheight=unit(0.5,"cm"), title.position="top", title.hjust=0.5)) +
    labs(x="BoneMarrowMap UMAP 1", y="BoneMarrowMap UMAP 2", title=glev,
         subtitle=sprintf("normal-like cells, n=%s", formatC(nrow(d), big.mark=","))) +
    theme_minimal(base_size=11) +
    theme(panel.grid=element_blank(), axis.line=element_line(colour="grey40", linewidth=0.4),
          axis.text=element_blank(), axis.ticks=element_blank(),
          axis.title=element_text(size=9),
          plot.title=element_text(size=12, face="bold", hjust=0.5),
          plot.subtitle=element_text(size=8, hjust=0.5, colour="grey30"),
          legend.position="bottom", legend.text=element_text(size=7.6),
          legend.title=element_text(size=8.6), legend.key.size=unit(0.3,"cm"),
          legend.margin=margin(-2,0,0,0), plot.margin=margin(4,6,2,4))
}
pH <- pseudo_panels("Healthy","g")
pHi<- pseudo_panels("RNHI-high AML","g")
pLo<- pseudo_panels("RNHI-low AML","g")

d <- ct[path=="P2" & !is.na(score)]
grp_cols <- c("Healthy"="#009E73","RNHI-high AML"="#D55E00","RNHI-low AML"="#CC79A7")
d[, grp := factor(grpDisp[group], levels=unname(grpDisp))]
d[, unit := factor(unit, levels=unit)]
kw_t <- tt[metric=="score"]
p_txt <- sprintf("KW p = %s", fmt2(kw_t$kw_p))
sig   <- tt[metric=="score"]
sigC  <- character(0)
if (!is.na(sig$p_HvHi_bh) && sig$p_HvHi_bh<0.05) sigC["Healthy\u2212RNHI-high"] <- sprintf("BH p = %s", fmt2(sig$p_HvHi_bh))
if (!is.na(sig$p_HvLo_bh) && sig$p_HvLo_bh<0.05) sigC["Healthy\u2212RNHI-low"] <- sprintf("BH p = %s", fmt2(sig$p_HvLo_bh))
if (!is.na(sig$p_HiLo_bh) && sig$p_HiLo_bh<0.05) sigC["RNHI-high\u2212RNHI-low"] <- sprintf("BH p = %s", fmt2(sig$p_HiLo_bh))
pB <- ggplot(d, aes(grp, score, fill=grp)) +
  geom_violin(scale="width", width=0.85, alpha=0.30, linewidth=0.4) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.85, linewidth=0.4) +
  geom_jitter(aes(colour=grp), width=0.10, size=1.5, alpha=0.9) +
  scale_fill_manual(values=grp_cols, name=NULL, labels=c("Healthy"="Healthy (n=54)","RNHI-high AML"="RNHI-high AML (n=5)","RNHI-low AML"="RNHI-low AML (n=5)")) +
  scale_colour_manual(values=grp_cols, guide="none") +
  scale_x_discrete(labels=c("Healthy"="Healthy","RNHI-high AML"="RNHI-high","RNHI-low AML"="RNHI-low")) +
  scale_y_continuous(breaks=seq(0,1,0.25)) +
  coord_cartesian(ylim=c(0,1.12), clip="off") +
  labs(x=NULL, y="Continuity score (0\u20131)",
       title=sprintf("Continuity score by group (%s)", p_txt),
       subtitle="P2 myelo-erythroid axis (HSC/MPP -> MyeProg -> Mono/Gran/Ery/Mega -> LateEry)\nper-patient pseudotime coverage & gap; n=300 downsampling x50; points = patients.",
       caption="Pairwise Wilcoxon, global BH: Healthy\u2212low BH p = 3.5e-03 (bracket)\nHealthy\u2212high ns & high\u2212low ns (raw p 2.1e-01 / 2.1e-01).") +
  theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=9.5, colour="black"),
        axis.text.y=element_text(size=8.6, colour="black"),
        axis.title=element_text(size=9.6),
        plot.title=element_text(size=11.5, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=7.8, hjust=0.5, colour="grey25", lineheight=1.05),
        plot.caption=element_text(size=7, hjust=0.5, colour="grey35"),
        legend.position="bottom", legend.text=element_text(size=8.2),
        legend.key.size=unit(0.4,"cm"), legend.margin=margin(1,0,0,0),
        plot.margin=margin(6,10,3,8))

if ("Healthy\u2212RNHI-low" %in% names(sigC)) {
  pB <- pB + ggsignif::geom_signif(comparisons=list(c("Healthy","RNHI-low AML")),
                                   annotations=sigC[["Healthy\u2212RNHI-low"]], y_position=1.04,
                                   tip_length=0.02, textsize=3.2, vjust=0, colour="#A50F15",
                                   fontface="bold")
}

emit <- function(p, base, w, h) {
  pdf_f<-file.path(res,paste0(base,".pdf")); png_f<-file.path(res,paste0(base,".png"))
  grDevices::cairo_pdf(pdf_f,width=w,height=h,family="Helvetica"); print(p); grDevices::dev.off()
  grDevices::png(png_f,width=w,height=h,units="in",res=300,type="cairo"); print(p); grDevices::dev.off()
  cat("Output completed.\n")
}
emit(pH,  "Figure6C_panelA_UMAP_Healthy", 4.5, 4.5)
emit(pHi, "Figure6C_panelA_UMAP_RNHI-high", 4.5, 4.5)
emit(pLo, "Figure6C_panelA_UMAP_RNHI-low", 4.5, 4.5)
rowA <- pH + pHi + pLo + plot_layout(ncol=3) +
  plot_annotation(title="Residual normal haematopoiesis coloured by BoneMarrowMap pseudotime",
                  theme=theme(plot.title=element_text(size=13, face="bold", hjust=0.5)))
emit(rowA, "Figure6C_panelA_pseudotime_UMAP", 12.6, 4.9)
emit(pB, "Figure6C_panelB_continuity_score", 4.5, 4.5)
comp <- (rowA / pB) + plot_layout(heights=c(1.15,1)) +
  plot_annotation(title="Figure 6C  Residual normal differentiation continuity: Healthy vs RNHI-high/low AML",
                  theme=theme(plot.title=element_text(size=13.5, face="bold", hjust=0.5)))
emit(comp, "Figure6C_composite", 12.6, 9.6)
invisible(list(A=rowA,B=pB))
