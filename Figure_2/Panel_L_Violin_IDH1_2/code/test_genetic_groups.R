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
IN <- "data/figure_2/Figure2E_input_patients.csv"
OUT <- "output/figure_2/genetic_groups"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

d <- utils::read.csv(IN, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
d <- d[, c("RNHI", "TP53", "NPM1", "FLT3_ITD", "IDH1_2", "ELN2017"), drop = FALSE]
d$RNHI <- as.numeric(d$RNHI)
fmt_p <- function(p) { p <- as.numeric(p); ifelse(is.na(p), "NA", formatC(signif(p,2), format="e", digits=1)) }
na <- function(x) is.na(x) | trimws(as.character(x)) %in% c("","NA","nan")

spec <- list(
  TP53     = list(col="TP53",   levels=c("Mut","WT"), three=FALSE),
  NPM1     = list(col="NPM1",   levels=c("Mut","WT"), three=FALSE),
  FLT3_ITD = list(col="FLT3_ITD",levels=c("Positive","Negative"), three=FALSE),
  IDH1_2   = list(col="IDH1_2", levels=c("Mut","WT"), three=FALSE),
  ELN2017  = list(col="ELN2017",levels=c("Favorable","Intermediate","Adverse"), three=TRUE)
)
sum_rows <- list(); test_rows <- list(); pair_rows <- list()
for (nm in names(spec)) {
  s <- spec[[nm]]; cc <- !na(d[[s$col]]) & is.finite(d$RNHI)
  sub <- d[cc, c("RNHI", s$col)]; names(sub)[2] <- "grp"
  sub$grp <- factor(sub$grp, levels=s$levels)
  n <- table(sub$grp); med <- tapply(sub$RNHI, sub$grp, median, na.rm=TRUE)
  q1 <- tapply(sub$RNHI, sub$grp, quantile, .25, na.rm=TRUE); q3 <- tapply(sub$RNHI, sub$grp, quantile, .75, na.rm=TRUE)
  for (g in s$levels) sum_rows[[length(sum_rows)+1L]] <- data.frame(variable=nm, group=g, n=as.integer(n[g]),
      median=med[[g]], q1=q1[[g]], q3=q3[[g]], stringsAsFactors=FALSE)
  if (s$three) {
    kv <- kruskal.test(RNHI ~ grp, data=sub)
    test_rows[[length(test_rows)+1L]] <- data.frame(variable=nm, test="kruskal", comparison="overall",
        n=sum(cc), p=kv$p.value, p_format=fmt_p(kv$p.value), stringsAsFactors=FALSE)
    pairs <- combn(s$levels, 2, simplify=FALSE)
    pv <- sapply(pairs, function(pr) {
      w <- suppressWarnings(wilcox.test(RNHI ~ grp, data=sub[sub$grp %in% pr,], exact=FALSE))
      w$p.value })
    padj <- p.adjust(pv, method="BH")
    for (i in seq_along(pairs)) pair_rows[[length(pair_rows)+1L]] <- data.frame(variable=nm,
        group1=pairs[[i]][1], group2=pairs[[i]][2], p_raw=pv[i], p_raw_format=fmt_p(pv[i]),
        p_BH=padj[i], p_BH_format=fmt_p(padj[i]), stringsAsFactors=FALSE)
  } else {
    w <- suppressWarnings(wilcox.test(RNHI ~ grp, data=sub, exact=FALSE))
    test_rows[[length(test_rows)+1L]] <- data.frame(variable=nm, test="wilcoxon", comparison=paste(s$levels, collapse=" vs "),
        n=sum(cc), p=w$p.value, p_format=fmt_p(w$p.value), stringsAsFactors=FALSE)
  }
}
sum_t <- do.call(rbind, sum_rows); test_t <- do.call(rbind, test_rows)
pair_t <- if(length(pair_rows)) do.call(rbind, pair_rows) else data.frame()
write.csv(sum_t, file.path(OUT,"Figure2E_summary_by_group.csv"), row.names=FALSE)
write.csv(test_t, file.path(OUT,"Figure2E_group_tests.csv"), row.names=FALSE)
if(nrow(pair_t)) write.csv(pair_t, file.path(OUT,"Figure2E_ELN_pairwise_BH.csv"), row.names=FALSE)
cat("== group summary ==\n"); print(sum_t, row.names=FALSE)
cat("\n== tests ==\n"); print(test_t, row.names=FALSE)
if(nrow(pair_t)){cat("\n== ELN pairwise ==\n"); print(pair_t, row.names=FALSE)}
cat("Output completed.\n")
