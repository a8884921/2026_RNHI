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

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table); library(ggbeeswarm)})
options(warn=-1)
W <- "output/figure_5/fragmentation"
O <- "output/figure_5/fragmentation"
dir.create(O, recursive = TRUE, showWarnings = FALSE)
fr <- fread(file.path(W,"per_sample_fragmentation_final.csv"))
pt <- fread(file.path(O,"figD2_group_tests_all.csv"))
metrics <- c("preserved_area_fraction","largest_cc_fraction","mean_component_size","boundary_density","morans_I_sRNHI")
titles  <- c("Preserved area fraction","Largest component / preserved","Mean island size","Boundary density","Moran's I of sRNHI")
ylabs   <- c("Fraction preserved","Fraction in largest CC","Spots per island","Boundary contacts per preserved spot","Moran's I")
pval <- setNames(pt$Welch_p, pt$metric)
fancy <- function(p) ifelse(p<0.001,"p = <1e-03",paste0("p = ",formatC(p,format="e",digits=2)))

grp_lab <- c(Healthy="Healthy (n = 3)", AML="AML (n = 4)")
pl <- list()
for(i in seq_along(metrics)){
  m <- metrics[i]
  d <- data.frame(group=fr$group, value=fr[[m]])
  d$grp <- factor(grp_lab[d$group], levels=unname(grp_lab[c("Healthy","AML")]))
  ymax <- max(d$value); pv <- pval[[m]]
  lbl <- if(pv<0.05) paste0(fancy(pv),"   *") else fancy(pv)
  p <- ggplot(d, aes(grp,value)) +
    geom_boxplot(aes(fill=group), width=0.42, outlier.shape=NA, alpha=0.18, linewidth=0.8, colour="grey25") +
    geom_beeswarm(aes(colour=group), cex=3.4, size=4, priority="density", alpha=0.95) +
    scale_fill_manual(values=c(Healthy="#00857C",AML="#C1272D"), guide="none") +
    scale_colour_manual(values=c(Healthy="#00857C",AML="#C1272D"), guide="none") +
    coord_cartesian(ylim=c(0, ymax*1.30), clip="off") +
    theme_classic(base_size=18) +
    theme(plot.title=element_text(size=15,face="bold",hjust=0.5,margin=margin(b=8)),
          axis.text.x=element_text(size=13,face="bold",colour="black"),
          axis.text.y=element_text(size=13), axis.title.x=element_blank(),
          axis.title.y=element_text(size=13.5), plot.margin=margin(10,10,6,8)) +
    labs(y=ylabs[i], title=titles[i]) +
    annotate("text", x=1.5, y=ymax*1.24, size=5.8, fontface="bold", label=lbl)
  pl[[i]] <- p
}
g <- wrap_plots(pl, ncol=2) + plot_layout(guides="collect")
ggsave(file.path(O,"FigD2_fragmentation_nature.pdf"), g, width=9.2, height=12.5, device=cairo_pdf)
ggsave(file.path(O,"FigD2_fragmentation_nature.png"), g, width=9.2, height=12.5, dpi=300, device=png, type="cairo")
cat("FIGD2 NATURE v4 DONE\n")
