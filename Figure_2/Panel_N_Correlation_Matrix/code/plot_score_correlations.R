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

options(warn=-1)
suppressMessages({library(ggplot2); library(ComplexHeatmap); library(circlize); library(grid); library(cowplot)})
R <- "data/figure_2/correlation"
D43 <- R
D31 <- R
D32 <- "output/figure_2/correlation"
OUT <- D32
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
st <- readRDS(file.path(D43, "Figure2G_stats.rds")); poolS <- st$poolS
mat <- readRDS(file.path(D31, "Figure2G_heatmap_input.rds"))$mat

cohort_cols <- c("cohort_1"="#E15759","cohort_2"="#4E79A7","cohort_3"="#59A14F","cohort_4"="#F28E2B","cohort_5"="#B07AA1")
poolS$cohort <- factor(poolS$cohort, levels=names(cohort_cols))
fmt_p <- function(p) ifelse(p < 1e-99, "<1.0e-99", formatC(signif(p,2), format="e", digits=1))

mk_scatter <- function(xv, main){
  ok <- is.finite(poolS[[xv]]) & is.finite(poolS$z_RNHI); dd <- poolS[ok,]
  ct <- suppressWarnings(cor.test(dd[[xv]], dd$z_RNHI, method="spearman", exact=FALSE))
  rho <- unname(ct$estimate); p <- ct$p.value; n <- nrow(dd)
  rho_show <- if (abs(rho) < 0.005) sprintf("%.3f", rho) else sprintf("%.2f", rho)
  p_show <- if (p < 5e-324) "<1.0e-99" else fmt_p(p)
  xlab <- if (grepl("iScore", xv)) "iScore (within-cohort z)" else "LSC17 (within-cohort z)"
  ggplot(dd, aes(.data[[xv]], z_RNHI)) +
    geom_point(aes(colour=cohort), size=1.7, alpha=0.55) +
    geom_smooth(method="lm", se=FALSE, colour="grey25", linewidth=0.7, linetype="dashed") +
    scale_colour_manual(values=cohort_cols, labels=paste("Cohort", seq_along(cohort_cols)), name="Cohort") +
    labs(x=xlab, y="RNHI (within-cohort z)", title=main,
         subtitle=bquote(italic(r)[s] * " = " * .(rho_show) * ",  " * italic(p) * " = " * .(p_show) * ";  n = " * .(n))) +
    theme_minimal(base_size=15) +
    theme(panel.grid=element_blank(), panel.border=element_rect(colour="grey40", fill=NA, linewidth=0.6),
          axis.title=element_text(size=14, face="bold"), axis.text=element_text(size=12),
          plot.title=element_text(size=15.5, face="bold", hjust=0.5),
          plot.subtitle=element_text(size=13, hjust=0.5, colour="#333333"),
          legend.position="none",
          plot.margin=margin(6,8,4,6))
}
pB <- mk_scatter("z_iScore_raw", "RNHI vs iScore")
pC <- mk_scatter("z_LSC17_raw", "RNHI vs LSC17")
pCleg <- pC + theme(legend.position="bottom", legend.title=element_text(size=12,face="bold"),
                    legend.text=element_text(size=11)) +
  guides(colour=guide_legend(nrow=1, override.aes=list(size=3.5, alpha=1)))
leg <- cowplot::get_legend(pCleg)

col_fun <- colorRamp2(c(-1,0,1), c("#2166AC","#F7F7F7","#B2182B"))
ht <- Heatmap(mat, name="rho", col=col_fun,
  rect_gp=gpar(col="white", lwd=2),
  cell_fun=function(j,i,x,y,w,h,fill){ grid.rect(x,y,w,h, gp=gpar(fill=fill, col="white", lwd=2)); grid.text(sprintf("%.2f", mat[i,j]), x, y, gp=gpar(fontsize=14, col="#1A1A1A", fontface="bold")) },
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_names_gp=gpar(fontsize=13, fontface="bold"), column_names_gp=gpar(fontsize=13, fontface="bold"),
  column_names_rot=30,
  heatmap_legend_param=list(title="Spearman rho", title_gp=gpar(fontsize=12,fontface="bold"),
                            labels_gp=gpar(fontsize=11), legend_height=unit(2.2,"in")),
  width=unit(3.4,"in"), height=unit(2.9,"in"))

render_heatmap_grob <- function(){
  grid.grabExpr(draw(ht, newpage=TRUE, heatmap_legend_side="bottom",
                     padding=unit(c(8,4,8,4),"mm")))
}
hm_grob <- render_heatmap_grob()

ggB <- ggplotGrob(pB); ggC <- ggplotGrob(pC)
row1 <- cowplot::plot_grid(hm_grob, ggB, ggC, ncol=3, labels=c("A","B","C"),
                           label_size=19, label_fontface="bold",
                           rel_widths=c(1.06,1,1), align="hv", axis="tb")
final <- cowplot::plot_grid(row1, leg, ncol=1, rel_heights=c(10,0.55))
ggsave(file.path(OUT,"Figure2G_main_A_heatmap_BC_scatter.pdf"), final, width=15.6, height=6.0, device=cairo_pdf, dpi=300)
ggsave(file.path(OUT,"Figure2G_main_A_heatmap_BC_scatter.png"), final, width=15.6, height=6.0, dpi=300)
cat("saved composite v2 pdf/png\n")
