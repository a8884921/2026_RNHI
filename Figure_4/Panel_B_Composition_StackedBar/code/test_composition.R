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

res <- "output/figure_4/composition"
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dt  <- fread(file.path(res, "Figure6B_patient_lineage_prop.csv"))
dt$unit <- anonymous_ids(dt$unit)
states <- c("Normal_HSC_MPP","Normal_LMPP","Normal_MyeloidProgenitor","Normal_Monocyte",
            "Normal_Granulocytic","Normal_Megakaryocytic","Normal_EarlyErythroid","Normal_LateErythroid",
            "Normal_B","Normal_T","Normal_NK","Normal_Dendritic")
lab12 <- c("Normal_HSC_MPP"="HSC/MPP","Normal_LMPP"="LMPP","Normal_MyeloidProgenitor"="Myeloid prog.",
           "Normal_Monocyte"="Monocyte","Normal_Granulocytic"="Granulocytic","Normal_Megakaryocytic"="Megakaryocytic",
           "Normal_EarlyErythroid"="Early erythroid","Normal_LateErythroid"="Late erythroid","Normal_B"="B",
           "Normal_T"="T","Normal_NK"="NK","Normal_Dendritic"="Dendritic")
grpLv <- c("Healthy","RNHI-high AML","RNHI-low AML")
dt[, grpDisp := factor(grpDisp, levels = grpLv)]
pairs <- list(c("Healthy","RNHI-high AML"), c("Healthy","RNHI-low AML"),
              c("RNHI-high AML","RNHI-low AML"))
pairsCol <- c("p_HvHi","p_HvLo","p_HiLo")

summ <- dt[, .(n_units = .N,
               mean_prop = mean(prop), median_prop = median(prop),
               sd_prop = sd(prop), q1 = quantile(prop, .25), q3 = quantile(prop, .75)),
           by = .(state12, grpDisp)]
summ[, mean_pct := 100 * mean_prop][, median_pct := 100 * median_prop]
summ[, iqr_prop := q3 - q1]
summ[, st_idx := match(state12, states)][, grp_idx := match(as.character(grpDisp), grpLv)]
setorder(summ, st_idx, grp_idx)
summ[, `:=`(st_idx = NULL, grp_idx = NULL)]
fwrite(summ, file.path(res, "Figure6B_lineage_group_summary.csv"))

fmt2 <- function(x) formatC(x, format = "e", digits = 1)
out <- list()
for (st in states) {
  sub <- dt[state12 == st]
  vals <- split(sub$prop, sub$grpDisp)
  kw  <- suppressWarnings(kruskal.test(vals[grpLv]))$p.value
  if (is.na(kw)) kw <- 1
  row <- data.table(state12 = st, kw_p = kw)
  row[, (pairsCol) := lapply(pairs, function(pr) {
    x <- vals[[pr[1]]]; y <- vals[[pr[2]]]
    if (all(x == y)) return(1)
    suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
  })]
  row[, kw_sig := kw < 0.05]
  out[[length(out) + 1L]] <- row
}
tab <- rbindlist(out)
tab[, kw_sig := kw_p < 0.05]

mat <- as.matrix(tab[kw_sig == TRUE, ..pairsCol])
allp <- as.vector(t(mat))
allp[is.na(allp)] <- 1
bh <- p.adjust(allp, method = "BH")

idx <- 0L
for (i in seq_len(nrow(tab))) {
  if (tab$kw_sig[i]) {
    for (j in seq_along(pairsCol)) {
      idx <- idx + 1L
      tab[i, paste0(pairsCol[j], "_bh") := bh[idx]]
    }
  } else {
    for (j in seq_along(pairsCol)) tab[i, paste0(pairsCol[j], "_bh") := NA_real_]
  }
}

for (pc in pairsCol) {
  bc <- paste0(pc, "_bh")
  tab[, (paste0(pc, "_txt")) := sapply(seq_len(nrow(tab)), function(i) {
    if (!tab$kw_sig[i] || is.na(tab[[bc]][i])) "na"
    else if (tab[[bc]][i] < 0.05) fmt2(tab[[bc]][i])
    else "ns"
  })]
  tab[, (paste0(pc, "_raw_txt")) := sapply(tab[[pc]], fmt2)]
}
tab[, kw_p_txt := sapply(kw_p, fmt2)]
tab[, st_idx := match(state12, states)]
setorder(tab, st_idx)
tab[, st_idx := NULL]
fwrite(tab, file.path(res, "Figure6B_lineage_tests.csv"))

for (st in states) {
  r <- tab[state12 == st]
  cat(sprintf("%-24s KW=%-9s %s | HvHi raw=%-9s bh=%-7s | HvLo raw=%-9s bh=%-7s | HiLo raw=%-9s bh=%-7s\n",
              lab12[[st]], r$kw_p_txt, ifelse(r$kw_sig, "sig", "ns "),
              r$p_HvHi_raw_txt, ifelse(is.na(r$p_HvHi_bh), "na", ifelse(r$p_HvHi_bh<0.05, r$p_HvHi_txt, "ns")),
              r$p_HvLo_raw_txt, ifelse(is.na(r$p_HvLo_bh), "na", ifelse(r$p_HvLo_bh<0.05, r$p_HvLo_txt, "ns")),
              r$p_HiLo_raw_txt, ifelse(is.na(r$p_HiLo_bh), "na", ifelse(r$p_HiLo_bh<0.05, r$p_HiLo_txt, "ns"))))
}
cat("[out] Figure6B_lineage_group_summary.csv\n[out] Figure6B_lineage_tests.csv\n")
