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

suppressPackageStartupMessages({library(qs); library(Seurat); library(GSVA); library(msigdbr); library(data.table)})
options(warn=-1)
W <- "output/figure_5/pathways"
dir.create(W, recursive = TRUE, showWarnings = FALSE)

hm <- msigdbr(species="Homo sapiens", category="H")
hallmark <- split(hm$gene_symbol, hm$gs_name)
wanted_hm <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_APOPTOSIS","HALLMARK_HYPOXIA","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
hallmark <- hallmark[names(hallmark) %in% c(wanted_hm)]

support_genes <- c("CXCL12","KITL","JAG1","ANGPT1","VCAM1","SCF","CSF1","IL7","KITLG","CDH2","POSTN","SPP1","LEPR","VEGFA","KDR","FLT1","NOTCH1","DLL4","PDGFRB","COL1A1","COL1A2","FN1")
extra <- list(HALLMARK_CXCL12_HEMATOPOIETIC_SUPPORT=support_genes)
gsets <- c(hallmark, extra)
cat("[gene sets]", length(gsets), "\n")
s2c <- setNames(as.character(1:4), as.character(1:4))
amld <- "data/figure_5/spatial"
out_rows <- list()
for(sid in names(s2c)){
  code <- s2c[sid]
  o <- qs::qread(file.path(amld,paste0(sid,".qs")), nthreads=4)
  cnt <- as.matrix(SeuratObject::LayerData(o, assay="Spatial", layer="counts"))

  md <- o[[]]; us <- if("RCTD_usable_for_composition" %in% colnames(md)) as.character(md$RCTD_usable_for_composition) else "TRUE"
  keep <- colnames(cnt)[us %in% c("TRUE","True")]
  cnt <- cnt[, keep, drop=FALSE]

  gs <- gsets[ vapply(gsets, function(g) sum(g %in% rownames(cnt)) >= 5, logical(1)) ]
  cat("[gsva]", code, "spots", ncol(cnt), "usable gene sets", length(gs), "\n")
  es <- GSVA::gsva(GSVA::ssgseaParam(cnt, gs), verbose=FALSE)
  sc <- data.table(cell_id=anonymous_ids(colnames(cnt))); for(p in rownames(es)) sc[[p]] <- as.numeric(es[p,])
  sc$code <- code
  out_rows[[code]] <- sc
  rm(o,cnt,es); gc()
}
P <- rbindlist(out_rows, fill=TRUE)
fwrite(P, file.path(W,"aml_pathway_scores.csv"))
cat("DONE step1; total", nrow(P), "pathways", ncol(P)-2, "\n")
