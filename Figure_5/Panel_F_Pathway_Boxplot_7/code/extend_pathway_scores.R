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

suppressPackageStartupMessages({library(qs); library(Seurat); library(GSVA); library(msigdbr);
                                library(data.table)})
options(warn = -1)
W <- "output/figure_5/pathways"
dir.create(W, recursive = TRUE, showWarnings = FALSE)
TMP <- "output/figure_5/pathways"
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
amld <- "data/figure_5/spatial"
NEW  <- file.path(TMP, "aml_scores_add3.csv")

hm <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- split(hm$gene_symbol, hm$gs_name)
ADD <- c("HALLMARK_GLYCOLYSIS","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_P53_PATHWAY")

s2c <- setNames(as.character(1:4), as.character(1:4))
if (!file.exists(NEW) || Sys.getenv("FORCE") == "1") {
  out <- list()
  for (sid in names(s2c)) {
    o  <- qs::qread(file.path(amld, paste0(sid, ".qs")), nthreads = 4)
    md <- o[[]]
    us <- if ("RCTD_usable_for_composition" %in% colnames(md)) as.character(md$RCTD_usable_for_composition) else "TRUE"
    cnt <- as.matrix(SeuratObject::LayerData(o, assay = "Spatial", layer = "counts"))
    keep <- colnames(cnt)[us %in% c("TRUE","True")]
    cnt <- cnt[, keep, drop = FALSE]
    nc <- colSums(cnt); nf <- colSums(cnt > 0)
    cpm <- t(t(cnt) / nc) * 1e4; lcpm <- log1p(cpm)
    gs <- hallmark[names(hallmark) %in% ADD]
    gs <- gs[vapply(gs, function(g) sum(g %in% rownames(lcpm)) >= 5, logical(1))]
    es <- GSVA::gsva(GSVA::ssgseaParam(lcpm, gs), verbose = FALSE)
    sc <- data.table(cell_id = anonymous_ids(colnames(lcpm)), code = s2c[sid], nCount = nc, nFeature = nf)
    for (p in rownames(es)) sc[[p]] <- as.numeric(es[p, ])
    out[[s2c[sid]]] <- sc
    cat("[add3]", s2c[sid], "spots", ncol(lcpm), "sets", length(gs), "\n")
    rm(o, cnt, cpm, lcpm, es); gc()
  }
  A3 <- rbindlist(out, fill = TRUE); fwrite(A3, NEW)
} else { A3 <- fread(NEW); cat("Output completed.\n") }

P0 <- fread(file.path(W, "aml_pathway_scores_cpm.csv"))
SIXD <- c("HALLMARK_E2F_TARGETS","HALLMARK_GLYCOLYSIS","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_P53_PATHWAY",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1","HALLMARK_APOPTOSIS")
P <- merge(P0, A3[, c("cell_id","code", setdiff(names(A3), names(P0))), with = FALSE],
           by = c("cell_id","code"), all.x = TRUE)
miss <- setdiff(SIXD, names(P)); stopifnot(length(miss) == 0)
cat("[13H] missing columns:", length(miss), "| spots", nrow(P), "\n")
KEEP_COLS <- c("cell_id","code","nCount","nFeature", SIXD, "HALLMARK_HYPOXIA")
fwrite(P[, ..KEEP_COLS], file.path(TMP, "aml_pathway_scores_cpm_13H.csv"))

LEG <- fread("data/figure_5/spatial/aml_niche_mapping_srnhi_legacy.csv")[, .(cell_id, code, z_legacy = sRNHI_z)]
GLB <- fread("data/figure_5/spatial/aml_niche_mapping_srnhi.csv")[
        , .(cell_id, code, z_global = sRNHI_z_global)]
D <- merge(merge(P, LEG, by = c("cell_id","code")), GLB, by = c("cell_id","code"))
cat(sprintf("[tier] rho(legacy per-niche z, global z) = %.4f\n", cor(D$z_legacy, D$z_global, method = "spearman")))

tier_of <- function(z) ifelse(z <= quantile(z, .25), "Low", ifelse(z >= quantile(z, .75), "High", "Middle"))
rows <- list()
for (cd in names(s2c)) for (zz in c("z_legacy","z_global")) {
  d <- D[code == cd]
  d[, tier := factor(tier_of(get(zz)), levels = c("Low","Middle","High"))]
  for (p in SIXD) {
    sub <- d[tier %in% c("Low","High")]
    vhi <- sub[[p]][sub$tier == "High"]; vlo <- sub[[p]][sub$tier == "Low"]
    w <- suppressWarnings(wilcox.test(vhi, vlo, exact = FALSE))
    hi <- median(vhi); lo <- median(vlo)
    rows[[paste(cd, zz, p)]] <- data.table(code = cd, zcol = zz, pathway = p,
      n_Low = length(vlo), n_High = length(vhi),
      median_Low = lo, median_High = hi, delta_High_minus_Low = hi - lo, wilcox_p = w$p.value)
  }
}
PW <- rbindlist(rows); PW[, fdr := p.adjust(wilcox_p, "BH"), by = .(zcol)]
fwrite(PW, file.path(W, "figF_spot_CPM_aligned_wilcoxon.csv"))

cat("\n=== Direction: High - Low (positive = higher with high sRNHI); global-z tier ===\n")
piv <- dcast(PW[zcol=="z_global"], pathway ~ code, value.var = "delta_High_minus_Low")
piv[, n_up := rowSums(.SD > 0), .SDcols = names(s2c)]
piv[, n_dn := rowSums(.SD < 0), .SDcols = names(s2c)]
print(piv[order(-n_up)])
cat("\nDONE\n")
