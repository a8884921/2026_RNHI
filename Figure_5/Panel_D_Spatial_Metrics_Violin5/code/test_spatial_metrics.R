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

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ggpubr); library(data.table)})
options(warn=-1)
W <- "output/figure_5/fragmentation"
O <- "output/figure_5/fragmentation"
dir.create(O, recursive = TRUE, showWarnings = FALSE)
fr <- fread(file.path(W,"per_sample_fragmentation_final.csv"))
fr$lab <- anonymous_ids(fr$sample_id)
metrics <- list(
  preserved_area_fraction="Preserved area fraction",
  largest_cc_fraction="Largest component / preserved",
  mean_component_size="Mean island size (official CC)",
  boundary_density="Boundary density",
  morans_I_sRNHI="Moran's I (sRNHI, Squidpy)")
fmtp <- function(p) formatC(p, format="e", digits=2)

plots <- list(); ptab <- list()
for(m in names(metrics)){
  d <- data.frame(group=factor(fr$group, levels=c("Healthy","AML")), value=fr[[m]], lab=fr$lab)

  wt <- wilcox.test(value~group, data=d, exact=FALSE)
  tt <- tryCatch(t.test(value~group, data=d), error=function(e) NULL)
  ptab[[m]] <- data.frame(metric=metrics[[m]], wilcox_p=wt$p.value,
                          welch_p=if(!is.null(tt)) tt$p.value else NA_real_,
                          healthy_n=sum(d$group=="Healthy"), aml_n=sum(d$group=="AML"))

  shown <- if(is.finite(wt$p.value)) wt$p.value else NA
  p <- ggplot(d, aes(group, value, colour=group)) +
    geom_jitter(width=0.08, size=3, alpha=0.85) +
    stat_summary(aes(group=group), fun=median, geom="crossbar", width=0.35, linewidth=0.5, colour="grey30") +
    geom_text(aes(label=lab), nudge_x=0.22, size=3.2, colour="grey25") +
    scale_colour_manual(values=c(Healthy="#00A087",AML="#E64B35")) +
    theme_minimal(base_size=13) +
    theme(legend.position="none", panel.grid.minor=element_blank(),
          axis.title.y=element_text(size=11), axis.text=element_text(size=11)) +
    labs(x="", y=metrics[[m]]) +
    annotate("text", x=1.5, y=max(d$value)*1.12,
             label=paste0("Wilcoxon p = ", fmtp(shown)), size=4.5, fontface="bold")
  plots[[m]] <- p
}
pt <- rbindlist(ptab)
pt$wilcox_p_e <- sapply(pt$wilcox_p, fmtp)
pt$welch_p_e  <- sapply(pt$welch_p, fmtp)
fwrite(pt, file.path(O,"figD2_group_tests.csv"))
print(pt[, .(metric, wilcox_p_e, welch_p_e, healthy_n, aml_n)])
g <- wrap_plots(plots, ncol=2)
ggsave(file.path(O,"FigD2_fragmentation_stats_pval.pdf"), g, width=9, height=13.5, device=cairo_pdf)
ggsave(file.path(O,"FigD2_fragmentation_stats_pval.png"), g, width=9, height=13.5, dpi=300, device=png, type="cairo")
cat("FIGD2 STATS DONE\n")
