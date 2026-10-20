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
                                library(data.table); library(ggplot2)})
options(warn = -1)

W <- "output/figure_5/pathways"
dir.create(W, recursive = TRUE, showWarnings = FALSE)
amld <- "data/figure_5/spatial"
SC0  <- as.data.frame(fread("output/figure_5/pathways/aml_pathway_scores.csv"))
C    <- as.data.frame(fread("data/figure_5/spatial/aml_niche_mapping_srnhi_legacy.csv"))

hm <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- split(hm$gene_symbol, hm$gs_name)
wanted <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
            "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1",
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
hallmark <- hallmark[names(hallmark) %in% wanted]
support_genes <- c("CXCL12","KITL","JAG1","ANGPT1","VCAM1","SCF","CSF1","IL7","KITLG","CDH2",
                   "POSTN","SPP1","LEPR","VEGFA","KDR","FLT1","NOTCH1","DLL4","PDGFRB",
                   "COL1A1","COL1A2","FN1")
gsets <- c(hallmark, list(HALLMARK_CXCL12_HEMATOPOIETIC_SUPPORT = support_genes))
cat("[gene sets]", length(gsets), "\n")

s2c <- setNames(as.character(1:4), as.character(1:4))
CACHE <- file.path(W, "aml_pathway_scores_cpm.csv")
if (file.exists(CACHE) && Sys.getenv("FORCE") != "1") {
  P <- fread(CACHE); cat("Output completed.\n"); SKIP <- TRUE
} else { SKIP <- FALSE; out <- list() }
if (!SKIP) for (sid in names(s2c)) {
  code <- s2c[sid]
  o <- qs::qread(file.path(amld, paste0(sid, ".qs")), nthreads = 4)
  md <- o[[]]
  us <- if ("RCTD_usable_for_composition" %in% colnames(md)) as.character(md$RCTD_usable_for_composition) else "TRUE"
  cnt <- as.matrix(SeuratObject::LayerData(o, assay = "Spatial", layer = "counts"))
  keep <- colnames(cnt)[us %in% c("TRUE","True")]
  cnt <- cnt[, keep, drop = FALSE]
  nc <- colSums(cnt); nf <- colSums(cnt > 0)

  cpm <- t(t(cnt) / nc) * 1e4
  lcpm <- log1p(cpm)
  gs <- gsets[vapply(gsets, function(g) sum(g %in% rownames(lcpm)) >= 5, logical(1))]
  es <- GSVA::gsva(GSVA::ssgseaParam(lcpm, gs), verbose = FALSE)
  sc <- data.table(cell_id = anonymous_ids(colnames(lcpm)), code = code, nCount = nc, nFeature = nf)
  for (p in rownames(es)) sc[[p]] <- as.numeric(es[p, ])
  out[[code]] <- sc
  cat("[CPM ssGSEA]", code, "spots", ncol(lcpm), "gene sets", length(gs), "\n")
  rm(o, cnt, cpm, lcpm, es); gc()
}
if (!SKIP) { P <- rbindlist(out, fill = TRUE)
fwrite(P, file.path(W, "aml_pathway_scores_cpm.csv")) }

rows <- list()
for (cd in names(s2c)) {
  d <- P[code == cd]
  for (p in names(gsets)) {
    if (!p %in% names(d)) next
    rho_cpm <- suppressWarnings(cor(d[["nCount"]], d[[p]], method = "spearman"))
    r <- merge(SC0[SC0$code == cd, ], P[code == cd, .(cell_id, code, nCount)], by = c("cell_id","code"))
    rho_raw <- if (p %in% names(r)) suppressWarnings(cor(as.numeric(r[["nCount"]]), as.numeric(r[[p]]), method = "spearman")) else NA_real_
    rows[[paste(cd, p)]] <- data.table(code = cd, pathway = p,
                                       rho_rawcount = round(rho_raw, 3), rho_CPM = round(rho_cpm, 3))
  }
}
DR <- rbindlist(rows); fwrite(DR, file.path(W, "figF_spot_cpm_depth_rho.csv"))
cat("\n[depth confounding] median absolute rho: raw-count =", round(median(abs(DR$rho_rawcount), na.rm = TRUE), 3),
    "  CPM =", round(median(abs(DR$rho_CPM), na.rm = TRUE), 3), "\n")

s2c2 <- setNames(as.character(1:4), as.character(1:4))
dat <- merge(P, C[, c("cell_id", "code", "sRNHI_z")], by = c("cell_id", "code"))
stopifnot(nrow(dat) == 5333)
dat$code <- factor(dat$code, levels = names(s2c2))
dat$tier <- NA_character_
for (cd in names(s2c2)) {
  z <- dat$sRNHI_z[dat$code == cd]; q <- quantile(z, probs = c(.25, .75))
  dat$tier[dat$code == cd] <- ifelse(z <= q[1], "Low", ifelse(z >= q[2], "High", "Middle"))
}
dat$tier <- factor(dat$tier, levels = c("Low", "Middle", "High"))

KEEP_up <- c("HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA","HALLMARK_INTERFERON_ALPHA_RESPONSE",
             "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB")
LAB <- c(HALLMARK_APOPTOSIS = "Apoptosis", HALLMARK_HYPOXIA = "Hypoxia",
         HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-\u03b1 resp.",
         HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-\u03b3 resp.",
         HALLMARK_TNFA_SIGNALING_VIA_NFKB  = "TNF\u03b1\u2013NF-\u03baB")
stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns")))
fmtp  <- function(p) ifelse(p < 2.2e-16, "p<2.2e-16", sprintf("p=%.2e", p))

L <- list(); PW <- list()
for (cd in names(s2c2)) {
  d <- dat[dat$code == cd]
  for (p in KEEP_up) {
    L[[paste(cd,p)]] <- data.table(code = cd, tier = d$tier,
                                   pathway = factor(unname(LAB[p]), levels = unname(LAB[KEEP_up])),
                                   score = d[[p]])
    sub <- d[tier %in% c("Low","High")]
    w <- suppressWarnings(wilcox.test(sub[[p]] ~ droplevels(sub$tier), exact = FALSE))
    hi <- median(sub[[p]][sub$tier == "High"]); lo <- median(sub[[p]][sub$tier == "Low"])
    PW[[paste(cd,p)]] <- data.table(code = cd, pathway = unname(LAB[p]),
                                    n_Low = sum(sub$tier == "Low"), n_High = sum(sub$tier == "High"),
                                    median_Low = lo, median_High = hi,
                                    delta_High_minus_Low = hi - lo,
                                    wilcox_p = w$p.value, stars = stars(w$p.value))
  }
}
ES <- rbindlist(L); PW <- rbindlist(PW)
PW[, fdr := p.adjust(wilcox_p, method = "BH")]
fwrite(PW, file.path(W, "figF3_spot_CPM_wilcoxon.csv"))

ann <- ES[tier %in% c("Low","High"), .(y = max(score, na.rm = TRUE)), by = .(code, pathway)]
ann <- merge(ann, PW[, .(code, pathway, lab = paste0(fmtp(wilcox_p), " ", stars))], by = c("code","pathway"))

g <- ggplot(ES, aes(tier, score, fill = tier)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.35, width = 0.68) +
  geom_text(data = ann, aes(x = 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.6, size = 2.6) +
  facet_grid(pathway ~ code, scales = "free_y") +
  scale_fill_manual(values = c(Low = "#74ADD1", Middle = "#F0F0F0", High = "#D73027"), guide = "none") +
  labs(x = "sRNHI tier (per-patient 25/50/75 quantiles)",
       y = "ssGSEA score (log1p-CPM input, kcdf=Gaussian)",
       title = "Figure F supplement | spot-level ssGSEA recomputed on log1p-CPM (depth-corrected at input)",
       subtitle = paste0("p = unpaired Wilcoxon Low vs High (exploratory only \u2014 spots are spatially autocorrelated, not independent replicates). ",
                         "Median |rho(score, nCount)| across pathways: raw-count input ",
                         round(median(abs(DR$rho_rawcount), na.rm = TRUE), 2), " -> CPM input ",
                         round(median(abs(DR$rho_CPM), na.rm = TRUE), 2), ".")) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 8.5, face = "bold"),
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(size = 8.5),
        plot.title = element_text(size = 11.5, face = "bold"),
        plot.subtitle = element_text(size = 7.4, colour = "grey25"),
        panel.grid.minor = element_blank())
ggsave(file.path(W, "FigF_spot_CPM_supplement.png"), g, width = 11, height = 11, dpi = 300, device = png, type = "cairo")
ggsave(file.path(W, "FigF_spot_CPM_supplement.pdf"), g, width = 11, height = 11, device = cairo_pdf)

cat("\n[done] CPM spot scores, supplementary plots, and p-value tables\n\n")
print(PW[order(pathway, code), .(code, pathway, delta = round(delta_High_minus_Low, 4),
                                 wilcox_p = signif(wilcox_p, 3), stars, fdr = signif(fdr, 3))])
cat("\n[pseudobulk direction comparison] (pseudobulk High-Low)\n")
dir <- fread(file.path(W, "high_minus_low_direction.csv"))
print(dir[pathway %in% KEEP_up, .(pathway, "1" = HL_1, "2" = HL_2, "3" = HL_3, "4" = HL_4)])
