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

cells <- fread("data/figure_4/6c_cells.tsv", sep="\t")
cells$donor <- anonymous_ids(cells$donor)
res   <- "output/figure_4/continuity"; dir.create(res, recursive=TRUE, showWarnings=FALSE)

P1 <- c("Normal_HSC_MPP","Normal_MyeloidProgenitor","Normal_EarlyErythroid","Normal_Megakaryocytic","Normal_LateErythroid")
P2 <- c(P1, "Normal_Monocyte","Normal_Granulocytic")
PTMIN <- 50; NBINS <- 50L; PMAX <- 22.5; NREP <- 50L; NSAMP <- 300L
w <- PMAX / NBINS

metric1 <- function(pt) {
  if (!length(pt)) return(list(cov=NA_real_, maxgap=NA_real_, score=NA_real_))
  h <- table(cut(pt, breaks=seq(0, PMAX, length.out=NBINS+1), include.lowest=TRUE))
  occ <- as.integer(h) > 0
  cov <- mean(occ)
  runs <- rle(occ)
  maxgap <- max(runs$lengths[!runs$values], default=0L)
  score <- cov * (1 - maxgap / NBINS)
  list(cov=cov, maxgap=maxgap, score=score)
}
calc_patient <- function(sub, onpath) {
  pt <- sub[state12 %in% onpath & is.finite(pseudotime), pseudotime]
  n  <- length(pt)
  if (n < PTMIN) return(data.table(cells=n, reps=0L, cov=NA_real_, maxgap=NA_real_, score=NA_real_))
  rep <- min(NREP, if (n > NSAMP) NREP else 1L)
  M <- rbindlist(lapply(seq_len(rep), function(i) {
    s <- if (n > NSAMP) sample(pt, NSAMP) else pt
    as.data.table(metric1(s))
  }))
  data.table(cells=n, reps=rep, cov=mean(M$cov, na.rm=TRUE), maxgap=mean(M$maxgap, na.rm=TRUE), score=mean(M$score, na.rm=TRUE))
}

units <- cells[, .(unit=unique(donor)), by=.(group, kind)][, .(group, kind, unit)]

cells[, unit := donor]

resl <- list()
for (pn in c("P1","P2")) {
  onpath <- if (pn=="P1") P1 else P2
  for (i in seq_len(nrow(units))) {
    sub <- cells[group==units$group[i] & kind==units$kind[i] & donor==units$unit[i]]
    m <- calc_patient(sub, onpath)
    resl[[length(resl)+1L]] <- data.table(path=pn, group=units$group[i], kind=units$kind[i], unit=units$unit[i], m)
  }
}
out <- rbindlist(resl)
setorder(out, path, group, unit)
fwrite(out, file.path(res, "Figure6C_patient_continuity.csv"))
fwrite(cells[, .N, by=.(group, state12)], file.path(res, "Figure6C_onpath_counts.csv"))

d <- dcast(out[!is.na(score)], path + group ~ ., fun.aggregate=list(n=function(x) sum(!is.na(x))), value.var="score")
cat("eligible patients by path/group (score non-NA):\n"); print(d, topn=20)
cat("cells kept total:", nrow(cells), "\n")
cat("[out] Figure6C_patient_continuity.csv\n[out] Figure6C_onpath_counts.csv\n")
