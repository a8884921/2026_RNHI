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

suppressPackageStartupMessages({
  library(qs); library(Seurat); library(Matrix); library(BoneMarrowMap)
})

src_qs <- "data/figure_4/projection/query.qs"
ref_file <- "data/figure_4/projection/reference.rds"
uwot_file <- "data/figure_4/projection/umap_model.uwot"
out_dir <- "output/figure_4/projection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(out_dir, "run.log")
logit <- function(...) { cat(..., "\n"); cat(paste0(..., "\n"), file = logf, append = TRUE) }
logit("Projection run log")
logit("start:", Sys.time())

logit("loading reference ...")
ref <- readRDS(ref_file)
ref$save_uwot_path <- uwot_file
logit("reference class:", class(ref)[1], "cells:", if (inherits(ref,"Seurat")) ncol(ref) else "n/a")

logit("loading query ...")
ob <- qread(src_qs, nthreads = 4)
stopifnot(inherits(ob, "Seurat"))
counts <- Seurat::GetAssayData(ob, assay = "RNA", layer = "counts")
if (is.null(counts) || nrow(counts) == 0) counts <- Seurat::GetAssayData(ob, assay = "RNA", slot = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(as(counts, "Matrix"), "dgCMatrix")
meta <- ob@meta.data
required_fields <- c("sample_id", "subject_id", "copykat", "anno_Symphony", "azimuth")
stopifnot(all(required_fields %in% names(meta)))
meta <- meta[, required_fields, drop = FALSE]
meta$sample_id <- anonymous_ids(meta$sample_id)
meta$subject_id <- anonymous_ids(meta$subject_id)
colnames(counts) <- anonymous_ids(colnames(counts))
rownames(meta) <- anonymous_ids(rownames(meta))
stopifnot(nrow(meta) == ncol(counts))
if (is.null(rownames(meta)) || !all(colnames(counts) %in% rownames(meta))) rownames(meta) <- colnames(counts)
meta <- meta[colnames(counts), , drop = FALSE]
meta$cell_id <- rownames(meta)

meta$sampleID <- meta$sample_id
batch_key <- "sampleID"
logit(sprintf("query: %d cells x %d genes; patients=%d", ncol(counts), nrow(counts), length(unique(meta$subject_id))))

n_cells <- ncol(counts)
block_size <- 20000L
block_ids <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / block_size))
logit(sprintf("blocks: %d (size<=%d)", length(block_ids), block_size))
parts <- vector("list", length(block_ids))
for (i in seq_along(block_ids)) {
  ii <- block_ids[[i]]
  logit(sprintf("map block %d/%d (%d cells)", i, length(block_ids), length(ii)), Sys.time())
  parts[[i]] <- BoneMarrowMap::map_Query(exp_query = counts[, ii, drop=FALSE],
                                         metadata_query = meta[ii, , drop=FALSE],
                                         ref_obj = ref, vars = batch_key, do_normalize = TRUE)
  parts[[i]] <- BoneMarrowMap::calculate_MappingError(parts[[i]], reference = ref, MAD_threshold = 2.5)
  gc(full = TRUE)
}

cm <- list()
cm$meta <- do.call(rbind, lapply(parts, function(x) x@meta.data))
cm$umap <- do.call(rbind, lapply(parts, function(x) Seurat::Embeddings(x, reduction="umap_projected")))
rm(parts, block_ids, counts); gc(full=TRUE)

md <- cm$meta; md$cell_id <- rownames(md)
md$umap_projected_1 <- cm$umap[,1]; md$umap_projected_2 <- cm$umap[,2]
md$subject_id <- anonymous_ids(md$subject_id)

oq <- meta
oi <- match(md$cell_id, rownames(oq))
md$copykat <- oq$copykat[oi]
md$anno_Symphony <- oq$anno_Symphony[oi]
md$azimuth <- oq$azimuth[oi]

tab <- data.frame(metric=c("cells","pass","fail"),
                  value=c(nrow(md), sum(md$mapping_error_QC=="Pass",na.rm=TRUE), sum(md$mapping_error_QC=="Fail",na.rm=TRUE)))
logit("QC:", paste(apply(tab,1,paste,collapse="="), collapse=" ; "))

keep_fields <- c("cell_id", "subject_id", "sample_id", "copykat", "anno_Symphony", "azimuth", "umap_projected_1", "umap_projected_2", "mapping_error_QC")
md <- md[, keep_fields, drop = FALSE]
md$cell_id <- anonymous_ids(md$cell_id)
rownames(md) <- md$cell_id
rownames(cm$umap) <- anonymous_ids(rownames(cm$umap))
saveRDS(list(meta=md, umap=cm$umap), file.path(out_dir, "projected_coords.rds"))
write.csv(tab, file.path(out_dir, "projection_summary.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(md[, c("cell_id","subject_id","sample_id","copykat","anno_Symphony","azimuth","umap_projected_1","umap_projected_2","mapping_error_QC")],
          file.path(out_dir, "projected_cells.csv"), row.names=FALSE, fileEncoding="UTF-8")
logit("saved rds + csvs. done:", Sys.time())
cat("[complete]\n")
