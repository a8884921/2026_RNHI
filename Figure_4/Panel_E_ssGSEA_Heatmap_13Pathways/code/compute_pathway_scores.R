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

suppressPackageStartupMessages({library(data.table); library(GSVA); library(edgeR); library(Matrix)})
R6D <- "../../../.."
pbdir <- "data/figure_4/pseudobulk"; merged <- "output/figure_4/ssgsea/merged"; dir.create(merged, recursive=TRUE, showWarnings=FALSE)
res <- "output/figure_4/ssgsea"; dir.create(res, recursive=TRUE, showWarnings=FALSE)
COMP <- c("Normal_HSC_MPP","Normal_MyeloidProgenitor","Normal_EarlyErythroid","Normal_EarlyMyeloidProg")
compL <- c("Normal_HSC_MPP"="HSC/MPP","Normal_MyeloidProgenitor"="MyeloidProgenitor","Normal_EarlyErythroid"="EarlyErythroid","Normal_EarlyMyeloidProg"="EarlyMyeloidProg")
HALL <- c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
          "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
          "HALLMARK_P53_PATHWAY","HALLMARK_APOPTOSIS","HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_GLYCOLYSIS","HALLMARK_MYC_TARGETS_V1")
hl <- readRDS("data/figure_4/hallmark_list.rds")

read_gz <- function(tag, cmp) {
  p <- file.path(pbdir, sprintf("%s_%s.csv.gz", tag, cmp))
  if (!file.exists(p)) return(NULL)
  m <- fread(p)
  genes <- m$gene; m[, gene := NULL]
  as.matrix(m) |> {function(x){rownames(x)<-genes; storage.mode(x)<-"integer"; x}}()
}
excl_all <- list()
for (cmp in COMP) {
  cat("== ", cmp, " ==\n")
  mats <- list(reference=read_gz("reference",cmp), validation=read_gz("validation",cmp), disease=read_gz("disease",cmp))
  mats <- mats[!vapply(mats, is.null, logical(1))]; if (!length(mats)) next
  gset <- Reduce(intersect, lapply(mats, rownames)); cat("common genes:", length(gset), "\n")
  cnt <- do.call(cbind, lapply(seq_along(mats), function(i) {
    x <- mats[[i]][gset, , drop=FALSE]
    colnames(x) <- anonymous_ids(colnames(x))
    x
  }))
  stopifnot(!anyDuplicated(colnames(cnt)))
  all_meta <- rbindlist(lapply(list.files(pbdir, pattern="_meta.csv$", full.names=TRUE), fread), fill=TRUE)
  stopifnot(all(c("sample", "group", "cells", "comp") %in% names(all_meta)))
  all_meta$sample <- anonymous_ids(all_meta$sample)
  allowed_groups <- c("Healthy", "RNHI-high AML", "RNHI-low AML")
  stopifnot(all(all_meta$group %in% allowed_groups))
  sample_groups <- unique(all_meta[, .(sample, group)])
  stopifnot(!anyDuplicated(sample_groups$sample))
  fwrite(sample_groups, file.path(res, "sample_groups.csv"))
  meta <- all_meta[comp == cmp, .(sample, group, cells)]
  stopifnot(!anyDuplicated(meta$sample))
  si <- merge(data.table(sample=colnames(cnt)), meta, by="sample", all.x=TRUE, sort=FALSE)
  stopifnot(!anyNA(si$group))
  MIN_CELLS <- as.integer(Sys.getenv("F6D_MIN_CELLS", "0"))
  excluded <- si[is.na(cells) | cells < MIN_CELLS]
  if (nrow(excluded)) excl_all[[length(excl_all)+1L]] <- cbind(data.table(comp=cmp), excluded)
  si <- si[!is.na(cells) & cells >= MIN_CELLS]
  cnt <- cnt[, si$sample, drop=FALSE]
  cat("samples kept:", nrow(si), " per group:", paste(si[, .N, by=group][, paste(group,N,sep="=")], collapse=", "), "\n")
  saveRDS(list(counts=cnt, info=si), file.path(merged, sprintf("merged_%s.rds", cmp)))

  dge <- edgeR::DGEList(cnt)
  cpm <- edgeR::cpm(dge, log=TRUE, prior.count=1)
  gsets <- hl[intersect(names(hl), HALL)]
  par <- GSVA::ssgseaParam(cpm, gsets)
  esc <- GSVA::gsva(par)
  out <- as.data.table(esc, keep.rownames="pathway")
  fwrite(out, file.path(res, sprintf("Figure6D_ssgsea_%s.csv", cmp)))
  cat("ssGSEA pathways:", nrow(out), "\n")

  si[, grp := factor(group, levels=c("Healthy","RNHI-high AML","RNHI-low AML"))]
  design <- model.matrix(~0+grp, data=si); colnames(design) <- make.names(gsub("grp","",colnames(design)))
  y <- DGEList(cnt); keep <- filterByExpr(y, design); y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y); y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, contrast=makeContrasts("RNHI.low.AML - Healthy", levels=design))
  tt <- topTags(qlf, n=200)
  t1 <- as.data.table(tt$table, keep.rownames="gene")
  fwrite(t1, file.path(res, sprintf("Figure6D_edgeR_lowVShealthy_top200_%s.csv", cmp)))
  cat("edgeR top200 saved for", cmp, "\n")
}

e <- rbindlist(lapply(excl_all, function(z) if (nrow(z)) z else NULL), fill = TRUE)
if (!nrow(e)) e <- data.table(comp=character(), sample=character(), group=character(),
                              cells=integer())
fwrite(e, file.path(res,"Figure6D_pseudobulk_excluded_samples.csv"))
cat("exclusions:", nrow(e), "\n")
cat("[done]\n")
