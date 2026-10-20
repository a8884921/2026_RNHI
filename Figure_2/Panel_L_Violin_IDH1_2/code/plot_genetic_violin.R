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
suppressMessages({library(ggplot2)})
IN <- "data/figure_2/Figure2E_input_patients.csv"
OUT <- "output/figure_2/genetic_groups"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

PANEL_IN <- as.numeric(Sys.getenv("PANEL_IN","3.2")); SC <- PANEL_IN/4.5
d <- utils::read.csv(IN, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
d <- d[, c("RNHI", "TP53", "NPM1", "FLT3_ITD", "IDH1_2", "ELN2017"), drop = FALSE]
d$RNHI <- as.numeric(d$RNHI)
na <- function(x) is.na(x) | trimws(as.character(x)) %in% c("","NA","nan")
fmt_p <- function(p) formatC(signif(p,2), format="e", digits=1)

spec <- list(
  TP53     = list(col="TP53",    levels=c("Mut","WT"), labels=c("TP53-mutated","TP53-WT"), cols=c("#E8836B","#7FB6D9")),
  NPM1     = list(col="NPM1",    levels=c("Mut","WT"), labels=c("NPM1-mutated","NPM1-WT"), cols=c("#E8836B","#7FB6D9")),
  FLT3_ITD = list(col="FLT3_ITD",levels=c("Positive","Negative"), labels=c("FLT3-ITD+","FLT3-ITD-"), cols=c("#E8836B","#7FB6D9")),
  IDH1_2   = list(col="IDH1_2",  levels=c("Mut","WT"), labels=c("IDH1/2-mutated","IDH1/2-WT"), cols=c("#E8836B","#7FB6D9")),
  ELN2017  = list(col="ELN2017", levels=c("Favorable","Intermediate","Adverse"),
                  labels=c("ELN favorable","ELN intermediate","ELN adverse"),
                  cols=c("#7FB6D9","#F0C36D","#E8836B"), xrot=45)
)

spec <- lapply(spec, function(z) { if (is.null(z$xrot)) z$xrot <- 0; z })

vals_all <- numeric(0)
for (nm in names(spec)) { cc <- !na(d[[spec[[nm]]$col]]) & is.finite(d$RNHI); vals_all <- c(vals_all, d$RNHI[cc]) }
rng <- range(vals_all); pad <- diff(rng)*0.03
YLIM <- c(rng[1]-pad, rng[2]+pad)
cat("shared YLIM:", round(YLIM,3), "\n")

pl <- list(); meta <- list()
for (nm in names(spec)) {
  s <- spec[[nm]]
  cc <- !na(d[[s$col]]) & is.finite(d$RNHI)
  sub <- data.frame(value=d$RNHI[cc], grp=d[[s$col]][cc])
  sub$grp <- factor(sub$grp, levels=s$levels)
  n <- as.integer(table(sub$grp)); names(n) <- s$levels
  med <- vapply(s$levels, function(g) median(sub$value[sub$grp==g], na.rm=TRUE), numeric(1))
  q1  <- vapply(s$levels, function(g) quantile(sub$value[sub$grp==g], .25, na.rm=TRUE, names=FALSE), numeric(1))
  q3  <- vapply(s$levels, function(g) quantile(sub$value[sub$grp==g], .75, na.rm=TRUE, names=FALSE), numeric(1))
  if (length(s$levels)==2) {
    w <- suppressWarnings(wilcox.test(value~grp, data=sub, exact=FALSE)); ptxt <- fmt_p(w$p.value); tname <- "Wilcoxon p"
  } else {
    k <- kruskal.test(value~grp, data=sub); ptxt <- fmt_p(k$p.value); tname <- "Kruskal-Wallis p"
  }
  sub$g <- factor(as.character(sub$grp), levels=s$levels, labels=s$labels)
  xlab_txt <- sprintf("%s (n=%d)", s$labels, n)

  px <- if (length(s$levels)==2) 1.5 else 2
  py <- YLIM[2] - diff(YLIM)*0.045
  gg <- ggplot(sub, aes(g, value, fill=g)) +
    geom_violin(trim=FALSE, width=0.92, colour="#8C8C8C", linewidth=0.4, alpha=0.5) +
    geom_jitter(aes(colour=g), width=0.16, size=0.9*SC, alpha=0.5) +
    stat_summary(fun=median, geom="point", shape=18, size=3.2*SC, colour="#404040") +
    stat_summary(fun.min=function(z) quantile(z,.25), fun.max=function(z) quantile(z,.75),
                 geom="errorbar", width=0.10, linewidth=0.8*SC, colour="#404040") +
    scale_fill_manual(values=setNames(s$cols, s$labels), drop=FALSE) +
    scale_colour_manual(values=setNames(s$cols, s$labels), drop=FALSE) +
    scale_x_discrete(labels=xlab_txt) +
    scale_y_continuous(limits=YLIM, breaks=seq(ceiling(YLIM[1]), floor(YLIM[2]), by=1)) +
    annotate("text", x=px, y=py, label=sprintf("%s = %s", tname, ptxt),
             hjust=0.5, vjust=1, size=4.6*SC, colour="#1F1F1F") +
    labs(x=NULL, y="RNHI") +
    theme_bw(base_size=17*SC) +
    theme(text=element_text(colour="#1F1F1F"),
          axis.title.y=element_text(size=18*SC),
          axis.text.x=element_text(size=12.5*SC, colour="#1F1F1F",
                                   angle=s$xrot, hjust=if (s$xrot>0) 1 else 0.5,
                                   vjust=if (s$xrot>0) 1 else 0.5),
          axis.text.y=element_text(size=12.5*SC),
          axis.ticks=element_line(colour="#1F1F1F", linewidth=0.5),
          panel.border=element_rect(colour="#1F1F1F", fill=NA, linewidth=0.7),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          legend.position="none", plot.margin=margin(4,8,4,8))
  ggsave(file.path(OUT, sprintf("Figure2E_violin_%s.pdf", nm)), gg, width=PANEL_IN, height=PANEL_IN, device=cairo_pdf, dpi=300)
  ggsave(file.path(OUT, sprintf("Figure2E_violin_%s.png", nm)), gg, width=PANEL_IN, height=PANEL_IN, dpi=300)
  saveRDS(gg, file.path(OUT, sprintf("Figure2E_violin_%s.rds", nm)))
  pl[[nm]] <- gg
  meta[[nm]] <- data.frame(variable=nm, group=s$labels, n=n, median=round(med,3), q1=round(q1,3), q3=round(q3,3),
                           test=tname, p_value=ptxt, stringsAsFactors=FALSE)
  cat("saved", nm, "|", tname, "=", ptxt, "\n")
}
m <- do.call(rbind, meta); write.csv(m, file.path(OUT,"Figure2E_violin_annot.csv"), row.names=FALSE)
saveRDS(list(pl=pl, meta=m), file.path(OUT,"Figure2E_panels.rds"))
cat("Output completed.\n")
