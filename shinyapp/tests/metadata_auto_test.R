# Regression test for the compact Seurat/Visium metadata input path.
args <- commandArgs(FALSE)
file_arg <- args[grep('^--file=', args)]
file_arg <- if (length(file_arg)) sub('^--file=', '', file_arg[[1]]) else 'tests/metadata_auto_test.R'
setwd(normalizePath(file.path(dirname(file_arg), '..')))
source('rnhi_core.R', local=TRUE)
ref <- parse_table(read_input('composition_A4.csv'))
health <- ref$x[ref$healthy,,drop=FALSE]

sc_meta <- data.frame(
  orig.ident=rep(c('Van_meta_01','Van_meta_02'), each=length(LINEAGES)),
  celltype=rep(LINEAGES, 2),
  umap_1=rep(seq_along(LINEAGES), 2),
  umap_2=rep(rev(seq_along(LINEAGES)), 2),
  stringsAsFactors=FALSE)
sc_comp <- prepare_modality_table(sc_meta, 'sc')
sc_r <- calculate_rnhi(sc_comp, health)
stopifnot(nrow(sc_comp)==18L, all(sc_comp$abundance==1), sc_r$n==2L)
sc_upper <- sc_meta
names(sc_upper)[names(sc_upper)=='umap_1'] <- 'UMAP_1'
names(sc_upper)[names(sc_upper)=='umap_2'] <- 'UMAP_2'
sc_upper_comp <- prepare_modality_table(sc_upper, 'sc')
sc_upper_r <- calculate_rnhi(sc_upper_comp, health)
sc_upper_r$plot_data <- make_sc_umap(sc_upper_r, sc_upper)
stopifnot(nrow(sc_upper_r$plot_data)==nrow(sc_upper),
          all.equal(sort(sc_upper_r$plot_data$umap_1), sort(as.numeric(sc_upper$UMAP_1))),
          all.equal(sort(sc_upper_r$plot_data$umap_2), sort(as.numeric(sc_upper$UMAP_2))))

sp_meta <- data.frame(
  spot_id=rep(sprintf('spot_%02d',1:4), each=length(LINEAGES)),
  celltype=rep(LINEAGES, 4),
  x=rep(rep(1:2, each=2), each=length(LINEAGES)),
  y=rep(rep(1:2, 2), each=length(LINEAGES)),
  stringsAsFactors=FALSE)
sp_comp <- prepare_modality_table(sp_meta, 'spatial')
sp_r <- calculate_rnhi(sp_comp, health)
sp_r$plot_data <- make_spatial_data(sp_r, sp_meta)
stopifnot(nrow(sp_comp)==36L, all(sp_comp$abundance==1), sp_r$n==4L, nrow(sp_r$plot_data)==4L)
sp_pxl <- sp_meta
names(sp_pxl)[names(sp_pxl)=='x'] <- 'pxl_col_in_fullres'
names(sp_pxl)[names(sp_pxl)=='y'] <- 'pxl_row_in_fullres'
sp_pxl_comp <- prepare_modality_table(sp_pxl, 'spatial')
sp_pxl_r <- calculate_rnhi(sp_pxl_comp, health)
sp_pxl_r$plot_data <- make_spatial_data(sp_pxl_r, sp_pxl)
stopifnot(all.equal(sort(unique(sp_pxl_r$plot_data$x)), sort(unique(as.numeric(sp_pxl$pxl_col_in_fullres)))),
          all.equal(sort(unique(sp_pxl_r$plot_data$y)), sort(unique(as.numeric(sp_pxl$pxl_row_in_fullres)))))
conflict <- sc_meta
conflict$healthy_reference <- 0L
conflict$healthy_reference[1] <- 1L
conflict_error <- tryCatch(prepare_modality_table(conflict, 'sc'), error=function(e) e)
stopifnot(inherits(conflict_error, 'error'))
cat('METADATA AUTO-ABUNDANCE PASS — sc samples:',sc_r$n,'; Visium spots:',sp_r$n,'\n')
