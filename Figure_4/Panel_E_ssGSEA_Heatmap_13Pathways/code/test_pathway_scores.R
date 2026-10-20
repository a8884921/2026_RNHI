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

res <- "output/figure_4/ssgsea"
COMP <- c("Normal_HSC_MPP","Normal_MyeloidProgenitor","Normal_EarlyErythroid","Normal_EarlyMyeloidProg")
fmt2 <- function(x) formatC(x, format="e", digits=1)
pc <- c("p_HvHi","p_HvLo","p_HiLo")

pbd  <- "data/figure_4/pseudobulk"
meta <- do.call(rbind, lapply(list.files(pbd, pattern="_meta.csv$", full.names=TRUE),
                              function(f) read.csv(f, stringsAsFactors=FALSE)))
stopifnot(all(c("sample", "group", "cells", "comp") %in% names(meta)))
meta$sample <- anonymous_ids(meta$sample)
rows <- list()
for (cmp in COMP) {
  f <- file.path(res, sprintf("Figure6D_ssgsea_%s.csv", cmp)); if (!file.exists(f)) next
  sc <- read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  paths <- sc[,1]; m <- as.matrix(sc[,-1]); rownames(m) <- paths
  cm <- meta[meta$comp == cmp, c("sample", "group", "cells")]
  stopifnot(!anyDuplicated(cm$sample))
  ids <- anonymous_ids(colnames(m))
  idx <- match(ids, cm$sample)
  stopifnot(!anyNA(idx))
  g <- unname(c("Healthy"="H", "RNHI-high AML"="hi", "RNHI-low AML"="lo")[cm$group[idx]])
  stopifnot(!anyNA(g))
  if (length(unique(g)) < 3) { cat("skip stats", cmp, "\n"); next }
  cl <- cm$cells[idx]
  n_H <- sum(g=="H"); n_hi <- sum(g=="hi"); n_lo <- sum(g=="lo")
  min_hi <- if (n_hi) min(cl[g=="hi"], na.rm=TRUE) else NA_integer_
  min_lo <- if (n_lo) min(cl[g=="lo"], na.rm=TRUE) else NA_integer_
  for (pw in paths) {
    x <- split(m[pw,], g)
    if (length(x$H) < 3 || length(x$hi) < 2 || length(x$lo) < 2) next
    kw <- suppressWarnings(kruskal.test(x[c("H","hi","lo")])$p.value); if (is.na(kw)) kw <- 1
    pr <- c(suppressWarnings(wilcox.test(x$H,x$hi,exact=FALSE)$p.value),
            suppressWarnings(wilcox.test(x$H,x$lo,exact=FALSE)$p.value),
            suppressWarnings(wilcox.test(x$hi,x$lo,exact=FALSE)$p.value))
    pr[is.na(pr)] <- 1
    rows[[length(rows)+1L]] <- data.frame(comp=cmp, pathway=pw, kw_p=kw, kw_sig=kw<0.05,
                                          p_HvHi=pr[1], p_HvLo=pr[2], p_HiLo=pr[3],
                                          n_healthy=n_H, n_high=n_hi, n_low=n_lo,
                                          min_cells_high=min_hi, min_cells_low=min_lo,
                                          stringsAsFactors=FALSE)
  }
}
t <- do.call(rbind, rows)
t$p_HvHi_bh <- NA_real_; t$p_HvLo_bh <- NA_real_; t$p_HiLo_bh <- NA_real_
for (cmp in unique(t$comp)) {
  id <- which(t$comp==cmp & t$kw_sig)
  if (!length(id)) next
  vec <- as.vector(t(as.matrix(t[id, pc]))); bh <- p.adjust(vec, "BH")
  t[id, c("p_HvHi_bh","p_HvLo_bh","p_HiLo_bh")] <- matrix(bh, ncol=3, byrow=TRUE)
}
txt <- function(raw, bh, sig) ifelse(!sig | is.na(bh), "na", ifelse(bh<0.05, fmt2(bh), "ns"))
t$kw_p_txt <- fmt2(t$kw_p)
for (j in seq_along(pc)) { b <- pc[j]; t[[paste0(b,"_raw_txt")]] <- fmt2(t[[b]]); t[[paste0(b,"_txt")]] <- txt(t[[b]], t[[paste0(b,"_bh")]], t$kw_sig) }
write.csv(t, file.path(res,"Figure6D_pathway_tests.csv"), row.names=FALSE)
cat("tests rows:", nrow(t), "\n")
o <- t[t$comp!="Normal_EarlyErythroid", c("comp","pathway","kw_p_txt","p_HvHi_txt","p_HvLo_txt","p_HiLo_txt")]
o$pathway <- sub("HALLMARK_","",o$pathway)
print(o, row.names=FALSE)
cat("[out] Figure6D_pathway_tests.csv\n")
