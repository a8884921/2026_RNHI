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

suppressPackageStartupMessages({library(data.table)})
argv <- commandArgs(trailingOnly=TRUE); PATH <- if (length(argv)) argv[1] else "P1"
res <- "output/figure_4/continuity"
d <- fread(file.path(res, "Figure6C_patient_continuity.csv"))[path==PATH & !is.na(score)]
G <- c("Healthy","AML_RNHI-high","AML_RNHI-low")
d[, grp := factor(group, levels=G)]
grpL <- c("Healthy"="Healthy","AML_RNHI-high"="RNHI-high AML","AML_RNHI-low"="RNHI-low AML")
metrics <- c("score","cov","maxgap")
fmt2 <- function(x) formatC(x, format="e", digits=1)
summ <- d[, .(n= .N, mean=mean(get("score")), med=median(get("score"))), by=group]
summ <- rbindlist(lapply(metrics, function(mt) {
  d[, .(metric=mt, n_units=.N, mean=mean(get(mt)), median=median(get(mt)), sd=sd(get(mt))), by=.(group)]
}))
fwrite(summ, file.path(res,"Figure6C_metric_group_summary.csv"))
pairs <- list(c("Healthy","AML_RNHI-high"),c("Healthy","AML_RNHI-low"),c("AML_RNHI-high","AML_RNHI-low"))
pairCol <- c("p_HvHi","p_HvLo","p_HiLo")
out <- list()
for (mt in metrics) {
  vals <- split(d[[mt]], d$grp)
  kw <- suppressWarnings(kruskal.test(vals[G]))$p.value; if(is.na(kw)) kw<-1
  r <- data.table(metric=mt, kw_p=kw, kw_sig=kw<0.05)
  r[, (pairCol) := lapply(pairs, function(pr){x<-vals[[pr[1]]];y<-vals[[pr[2]]]; if(all(x==y)) 1 else suppressWarnings(wilcox.test(x,y,exact=FALSE)$p.value)})]
  out[[length(out)+1L]] <- r
}
t <- rbindlist(out)
mat <- as.matrix(t[kw_sig==TRUE, ..pairCol]); allp <- as.vector(t(mat)); bh <- p.adjust(allp[!is.na(allp)], "BH"); k <- 0L
for (i in seq_len(nrow(t))) if (t$kw_sig[i]) for (j in seq_along(pairCol)) { k<-k+1L; t[i, paste0(pairCol[j],"_bh") := bh[k]] } else for (j in seq_along(pairCol)) t[i, paste0(pairCol[j],"_bh") := NA_real_]
for (pc in pairCol) {
  bc<-paste0(pc,"_bh"); t[, (paste0(pc,"_txt")) := sapply(seq_len(nrow(t)), function(i)
    if(!t$kw_sig[i] || is.na(t[[bc]][i])) "na" else if(t[[bc]][i]<0.05) fmt2(t[[bc]][i]) else "ns")]
  t[, (paste0(pc,"_raw_txt")) := sapply(t[[pc]], fmt2)]
}
t[, kw_p_txt := sapply(kw_p, fmt2)]
fwrite(t, file.path(res,"Figure6C_continuity_tests.csv"))
print(t[, .(metric, kw_p_txt, p_HvHi_raw_txt, p_HvHi_txt, p_HvLo_raw_txt, p_HvLo_txt, p_HiLo_raw_txt, p_HiLo_txt)])
cat("[out] Figure6C_continuity_tests.csv / Figure6C_metric_group_summary.csv\n")
