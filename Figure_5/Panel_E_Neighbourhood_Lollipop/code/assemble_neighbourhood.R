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

suppressPackageStartupMessages({library(qs); library(Seurat); library(data.table); library(spdep)})
options(warn=-1)
W <- "output/figure_5/neighbourhood"
dir.create(W, recursive = TRUE, showWarnings = FALSE)
frozen <- fread("data/figure_5/spatial/aml_niche_mapping_srnhi_legacy.csv")
THRESH <- -1.256
frozen[, state := ifelse(sRNHI_z >= THRESH, "preserved","impaired")]
niche_cols <- c(MSC="RCTD_weight_Niche_MSC_CARlike", Endo="RCTD_weight_Niche_Endothelial", Adipo="RCTD_weight_Niche_Adipocyte")
code_sid <- setNames(as.character(1:4), as.character(1:4))
hdir <- "data/figure_5/spatial"
read_niche <- function(path, code){
  o <- qs::qread(path, nthreads=4); md <- o[[]]
  us <- if("RCTD_usable_for_composition" %in% colnames(md)) as.character(md$RCTD_usable_for_composition) else "TRUE"
  d <- data.table(sample_id=code, cell_id=anonymous_ids(rownames(md)))
  for(short in names(niche_cols)) d[[short]] <- as.numeric(md[[niche_cols[short]]])
  d <- d[us %in% c("TRUE","True")]
}
H <- rbind(read_niche(file.path(hdir,"5.qs"),"5"),
           read_niche(file.path(hdir,"6.qs"),"6"),
           read_niche(file.path(hdir,"7.qs"),"7"))
p75 <- H[, lapply(.SD, quantile, .75), .SDcols=c("MSC","Endo","Adipo")]
cat("[Healthy P75] MSC=", round(p75$MSC,4), " Endo=", round(p75$Endo,4), " Adipo=", round(p75$Adipo,4), "\n")
fwrite(p75, file.path(W,"healthy_p75_thresholds.csv"))
amld <- "data/figure_5/spatial"
rows <- list()
for(code in names(code_sid)){
  sid <- code_sid[code]
  o <- qs::qread(file.path(amld,paste0(sid,".qs")), nthreads=4); md <- o[[]]
  us <- if("RCTD_usable_for_composition" %in% colnames(md)) as.character(md$RCTD_usable_for_composition) else "TRUE"
  d <- data.table(code=code, cell_id=anonymous_ids(rownames(md)))
  for(short in names(niche_cols)) d[[short]] <- as.numeric(md[[niche_cols[short]]])
  d <- d[us %in% c("TRUE","True")]
  fr <- frozen[, .(code,cell_id,x,y,sRNHI_z,state,Normal_abundance)]
  m <- merge(d, fr, by=c("code","cell_id")); m <- m[is.finite(x) & is.finite(y)]
  xy <- as.matrix(m[,.(x,y)])
  nb <- spdep::tri2nb(xy)
  adj <- lapply(seq_len(length(nb)), function(i) as.integer(nb[[i]]))
  for(short in names(niche_cols)){
    v <- m[[short]]
    m[[paste0("neighborhood_",short)]] <- vapply(seq_len(nrow(m)), function(i){
      j <- adj[[i]]; if(!length(j)) return(NA_real_); mean(v[j]) }, numeric(1))
    m[[paste0(short,"_rich")]] <- v >= p75[[short]]
  }
  m$code <- code; rows[[code]] <- m
  cat("[aml]", code, nrow(m), "usable\n")
}
A <- rbindlist(rows)
fwrite(A, file.path(W,"aml_niche_microenv_input.csv"))
cat("DONE STEP1b; total", nrow(A), "\n")
