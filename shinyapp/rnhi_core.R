# Frozen A4: nine axes, per-group detection limits, frac=0.65.
LINEAGES <- c('MyeProg','Mono','Gran','Mega','EarlyEry','LateEry','B','NK','DC')
MUSIC_NORMAL_MAP <- c(
  Normal_MyeloidProgenitor='MyeProg', Normal_Monocyte='Mono',
  Normal_Granulocytic='Gran', Normal_Megakaryocytic='Mega',
  Normal_EarlyErythroid='EarlyEry', Normal_LateErythroid='LateEry',
  Normal_B='B', Normal_NK='NK', Normal_Dendritic='DC')

validate_frac <- function(frac) {
  if (!is.numeric(frac) || length(frac) != 1L || !is.finite(frac) || frac <= 0 || frac >= 1)
    stop('Zero replacement fraction must be one finite number between 0 and 1.')
  frac
}
numeric_matrix <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0))
    stop('Abundances must be finite, numeric and non-negative; blanks/NA are not accepted.')
  if (!nrow(x) || ncol(x) != length(LINEAGES)) stop('A non-empty matrix with all nine lineages is required.')
  x
}
closure <- function(x) {
  x <- numeric_matrix(x)
  if (any(rowSums(x) <= 0)) stop('At least one sample has zero total across the nine lineages; RNHI is undefined for that sample.')
  x / rowSums(x)
}
multRepl <- function(x, frac=.65) {
  frac <- validate_frac(frac)
  x <- numeric_matrix(x)
  dl <- vapply(seq_len(ncol(x)), function(j) {
    positive <- x[x[,j] > 0,j]
    if (length(positive)) min(positive) else NA_real_
  }, numeric(1))
  if (anyNA(dl)) stop('No detection limit for an all-zero lineage: ', paste(LINEAGES[is.na(dl)], collapse=', '))
  out <- x
  for (i in seq_len(nrow(x))) {
    zero <- which(x[i,] == 0)
    if (!length(zero)) next
    nonzero <- which(x[i,] > 0)
    replacement <- frac * dl[zero]
    scale <- (sum(x[i,]) - sum(replacement)) / sum(x[i,nonzero])
    if (!is.finite(scale) || scale <= 0) stop('Zero replacement has no positive remaining mass for row ',i,'.')
    out[i,zero] <- replacement
    out[i,nonzero] <- x[i,nonzero] * scale
  }
  out
}
clr <- function(x) { lx <- log(x); sweep(lx, 1, rowMeans(lx), '-') }
center_from <- function(h, frac=.65) colMeans(clr(multRepl(closure(h),frac)))
score_rnhi <- function(x, health, frac=.65) {
  z <- clr(multRepl(closure(x),frac))
  delta <- sweep(z, 2, center_from(health,frac), '-')
  -sqrt(rowSums(delta^2))
}
read_input <- function(path, name=basename(path)) {
  # Base connections avoid readr path re-encoding on macOS with a C locale.
  sep <- if (grepl('\\.(tsv|txt)$',name,ignore.case=TRUE)) '\t' else ','
  dat <- read.table(path, header=TRUE, sep=sep, check.names=FALSE,
                    stringsAsFactors=FALSE, comment.char='', quote='"', fileEncoding='UTF-8-BOM')
  if (!nrow(dat)) stop('The uploaded table has no data rows.')
  dat
}
expression_gene_column <- function(dat) {
  hits <- intersect(c('gene','genes','gene_id','geneid','feature','features'), names(dat))
  if (length(hits)) hits[[1]] else names(dat)[[1]]
}
is_expression_matrix <- function(dat) {
  if (any(LINEAGES %in% names(dat)) || all(c('cell_type','abundance') %in% names(dat))) return(FALSE)
  if (!ncol(dat) || ncol(dat) < 3) return(FALSE)
  gcol <- expression_gene_column(dat); sample_cols <- setdiff(names(dat),gcol)
  if (length(sample_cols) < 2) return(FALSE)
  ok <- vapply(dat[sample_cols], function(v) { z <- suppressWarnings(as.numeric(v)); any(is.finite(z)) && !all(is.na(z)) }, logical(1))
  sum(ok) >= 2
}
parse_expression_matrix <- function(dat) {
  if (!is_expression_matrix(dat)) stop('This table is not a recognizable gene-by-sample expression matrix.')
  gcol <- expression_gene_column(dat); sample_cols <- setdiff(names(dat),gcol)
  genes <- trimws(as.character(dat[[gcol]]))
  if (anyNA(genes) || any(!nzchar(genes))) stop('The expression matrix gene column contains blank IDs.')
  if (anyDuplicated(sample_cols)) stop('Expression sample columns must have unique names.')
  mat <- do.call(cbind,lapply(dat[sample_cols], function(v) suppressWarnings(as.numeric(v))))
  colnames(mat) <- sample_cols
  if (any(!is.finite(mat)) || any(mat < 0)) stop('Expression values must be finite and non-negative.')
  if (anyDuplicated(genes)) { mat <- rowsum(mat, group=genes, reorder=FALSE); genes <- rownames(mat) } else rownames(mat) <- genes
  if (!nrow(mat) || !ncol(mat)) stop('The expression matrix has no usable genes or samples.')
  list(matrix=mat, sample_ids=colnames(mat), gene_column=gcol)
}
music_fallback_fit <- function(y, D, iter.max=1000, eps=.01) {
  # Deterministic non-negative, sum-to-one projection used only when the
  # Bioconductor MuSiC package is unavailable on a deployment host.
  k <- ncol(D); p <- rep(1/k, k)
  DtD <- crossprod(D); Dty <- as.numeric(crossprod(D, y))
  for (iter in seq_len(iter.max)) {
    den <- as.numeric(DtD %*% p) + 1e-12
    p_new <- p * pmax(Dty, 0) / den
    if (!all(is.finite(p_new)) || sum(p_new) <= 0) p_new <- rep(1/k, k)
    p_new <- p_new / sum(p_new)
    if (max(abs(p_new - p)) < eps) { p <- p_new; break }
    p <- p_new
  }
  fitted <- as.numeric(D %*% p); ss_res <- sum((y - fitted)^2); ss_tot <- sum((y - mean(y))^2)
  list(p.weight=p, R.squared=if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_, converge='fallback_projected_gradient')
}
first_column <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit)) hit[[1]] else NULL
}
metadata_to_composition <- function(dat, mode=c('sc','spatial')) {
  mode <- match.arg(mode)
  cell_col <- first_column(dat, c('celltype','cell_type','CellType','cell.type'))
  id_candidates <- if (mode == 'sc') c('sample_id','orig.ident','sample','patient_id') else c('spot_id','sample_id','barcode','orig.ident')
  id_col <- first_column(dat, id_candidates)
  if (is.null(cell_col) || is.null(id_col)) {
    req <- if (mode == 'sc') 'celltype + sample_id (or orig.ident)' else 'celltype + spot_id (or sample_id)'
    stop(mode, ' input requires ', req, '. Export Seurat/Visium meta.data with one row per annotated cell/spot.')
  }
  sid <- trimws(as.character(dat[[id_col]])); ct <- trimws(as.character(dat[[cell_col]]))
  if (anyNA(sid) || any(!nzchar(sid))) stop(id_col, ' contains blank IDs.')
  if (anyNA(ct) || any(!nzchar(ct))) stop(cell_col, ' contains blank cell-type labels.')
  if (any(!ct %in% LINEAGES)) {
    bad <- unique(ct[!ct %in% LINEAGES])
    stop('Unknown cell-type labels: ', paste(utils::head(bad, 8L), collapse=', '),
         '. Use the nine frozen names: ', paste(LINEAGES, collapse=', '), '.')
  }
  d <- data.frame(sample_id=sid, cell_type=ct, abundance=1, stringsAsFactors=FALSE)
  out <- aggregate(abundance ~ sample_id + cell_type, d, sum)
  # Preserve a stable healthy flag when it is present in exported metadata.
  hcol <- first_column(dat, c('healthy_reference','healthy','is_healthy'))
  if (!is.null(hcol)) {
    hv <- dat[[hcol]]
    if (is.logical(hv)) hv <- as.integer(hv)
    if (is.character(hv)) hv <- suppressWarnings(as.numeric(trimws(hv)))
    if (!is.numeric(hv) || anyNA(hv) || any(!hv %in% c(0,1))) stop(hcol, ' must contain only 0 or 1.')
    hgroups <- split(as.integer(hv), sid)
    if (any(vapply(hgroups, function(v) length(unique(v)) != 1L, logical(1))))
      stop('Conflicting ', hcol, ' values for one sample/spot.')
    hf <- data.frame(sample_id=names(hgroups), healthy_reference=vapply(hgroups, `[`, integer(1), 1L), stringsAsFactors=FALSE)
    out <- merge(out, hf, by='sample_id', all.x=TRUE, sort=FALSE)
  }
  # Coordinates are carried through for spatial plots. Keep the first coordinate
  # per spot/sample so repeated cell annotations do not create duplicate points.
  coord_candidates <- list(x=c('x','X','pxl_col_in_fullres'), y=c('y','Y','pxl_row_in_fullres'),
                           imagecol=c('imagecol','imagerow_x'), imagerow=c('imagerow','imagerow_y'),
                           umap_1=c('umap_1','UMAP_1'), umap_2=c('umap_2','UMAP_2'))
  coord_names <- names(coord_candidates)[vapply(coord_candidates, function(x) !is.null(first_column(dat,x)), logical(1))]
  if (length(coord_names)) {
    keep <- data.frame(sample_id=sid, stringsAsFactors=FALSE)
    for (nm in coord_names) keep[[nm]] <- as.numeric(dat[[first_column(dat,coord_candidates[[nm]])]])
    keep <- keep[!duplicated(keep$sample_id),,drop=FALSE]
    out <- merge(out, keep, by='sample_id', all.x=TRUE, sort=FALSE)
  }
  out
}
prepare_modality_table <- function(dat, mode=c('sc','spatial')) {
  mode <- match.arg(mode)
  if (all(LINEAGES %in% names(dat)) || all(c('cell_type','abundance') %in% names(dat))) return(dat)
  metadata_to_composition(dat, mode)
}
music_deconvolve <- function(dat, basis_path=file.path('data','music_basis.rds'), iter.max=1000, eps=.01) {
  # Keep MuSiC optional at deployment time: resolve the package name lazily so
  # rsconnect does not force a site-wide Bioconductor snapshot just to run the
  # bundled deterministic fallback.
  music_pkg <- paste0('Mu','SiC')
  use_music <- requireNamespace(music_pkg,quietly=TRUE)
  music_iter <- if (use_music) getExportedValue(music_pkg,'music.iter') else NULL
  if (!file.exists(basis_path)) stop('Bulk expression was detected, but the bundled locked MuSiC basis is missing.')
  expr <- parse_expression_matrix(dat); basis <- readRDS(basis_path)
  common <- intersect(rownames(basis$D),rownames(expr$matrix)); min_common <- max(20L,ceiling(.10*min(length(common),nrow(basis$D))))
  if (length(common) < min_common) stop('Too few genes overlap the locked MuSiC basis (',length(common),' common genes; need at least ',min_common,').')
  D <- basis$D[common,,drop=FALSE]; Sigma <- basis$Sigma[common,,drop=FALSE]; S <- basis$S
  rows <- vector('list',ncol(expr$matrix)); qc <- vector('list',ncol(expr$matrix))
  for (i in seq_len(ncol(expr$matrix))) {
    y <- expr$matrix[common,i]; pos <- is.finite(y) & y > 0; sid <- colnames(expr$matrix)[i]
    if (sum(pos) < max(20L,ceiling(.20*length(common)))) { qc[[i]] <- data.frame(sample_id=sid,common_genes=length(common),positive_genes=sum(pos),status='failed',reason='too_few_positive_genes',stringsAsFactors=FALSE); next }
    yy <- y[pos]; names(yy) <- common[pos]
    fit <- tryCatch(if (use_music) music_iter(Y=yy,D=D[pos,,drop=FALSE],S=S,Sigma=Sigma[pos,,drop=FALSE],iter.max=iter.max,nu=1e-4,eps=eps) else music_fallback_fit(yy,D[pos,,drop=FALSE],iter.max=iter.max,eps=eps), error=function(e)e)
    if (inherits(fit,'error')) { qc[[i]] <- data.frame(sample_id=sid,common_genes=length(common),positive_genes=sum(pos),status='failed',reason=conditionMessage(fit),stringsAsFactors=FALSE); next }
    p <- fit$p.weight
    if (is.null(p) || any(!is.finite(p)) || sum(p) <= 0) { qc[[i]] <- data.frame(sample_id=sid,common_genes=length(common),positive_genes=sum(pos),status='failed',reason='invalid_MuSiC_proportion',stringsAsFactors=FALSE); next }
    p <- p/sum(p); names(p) <- colnames(D)
    normal <- vapply(LINEAGES,function(ct) sum(p[names(MUSIC_NORMAL_MAP)[MUSIC_NORMAL_MAP==ct]],na.rm=TRUE),numeric(1))
    rows[[i]] <- data.frame(sample_id=sid,as.list(normal),healthy_reference=0,stringsAsFactors=FALSE,check.names=FALSE)
    qc[[i]] <- data.frame(sample_id=sid,common_genes=length(common),positive_genes=sum(pos),status='success',R_squared=as.numeric(fit$R.squared),converge=as.character(fit$converge),stringsAsFactors=FALSE)
  }
  good <- !vapply(rows,is.null,logical(1))
  if (!all(good)) {
    failed_qc <- do.call(rbind,qc[!good])
    detail <- paste0(failed_qc$sample_id, ' (', failed_qc$reason, ')')
    detail <- paste(utils::head(detail, 6L), collapse='; ')
    more <- if (length(detail) && nrow(failed_qc) > 6L) paste0('; +', nrow(failed_qc)-6L, ' more') else ''
    stop('MuSiC failed for ', nrow(failed_qc), ' sample(s); no partial RNHI result was returned. ', detail, more, '.')
  }
  comp <- do.call(rbind,rows); rownames(comp) <- NULL
  list(data=comp,qc=do.call(rbind,qc),expression=expr,basis_version=basis$version,common_genes=length(common),method=if (use_music) 'MuSiC deconvolution' else 'Locked-basis fallback projection (MuSiC unavailable)')
}
parse_table <- function(dat) {
  if (anyDuplicated(names(dat))) stop('Duplicate column names are not accepted.')
  if (!'sample_id' %in% names(dat)) {
    if ('spot_id' %in% names(dat)) dat$sample_id <- dat$spot_id
    else stop('Add a sample_id column (or spot_id for spot-level compositions).')
  }
  dat$sample_id <- as.character(dat$sample_id)
  if (anyNA(dat$sample_id) || any(!nzchar(trimws(dat$sample_id)))) stop('Sample IDs must not be blank.')
  ids <- unique(dat$sample_id)
  flags <- rep(FALSE,length(ids))
  if ('healthy_reference' %in% names(dat)) {
    flag_values <- dat$healthy_reference
    if (is.logical(flag_values)) flag_values <- as.integer(flag_values)
    if (is.character(flag_values)) flag_values <- suppressWarnings(as.numeric(trimws(flag_values)))
    if (!is.numeric(flag_values) || anyNA(flag_values) || any(!flag_values %in% c(0,1))) stop('healthy_reference must contain only 0 or 1.')
    groups <- split(flag_values,dat$sample_id)
    if (any(vapply(groups,function(v) length(unique(v)) != 1L,logical(1)))) stop('Conflicting healthy_reference values for one sample.')
    flags <- as.logical(flag_values[match(ids,dat$sample_id)])
  }
  if (all(LINEAGES %in% names(dat))) {
    if (anyDuplicated(dat$sample_id)) stop('Wide input requires a unique sample_id per row.')
    x <- as.matrix(dat[,LINEAGES,drop=FALSE])
  } else if (all(c('cell_type','abundance') %in% names(dat))) {
    if (anyNA(dat$cell_type) || any(!dat$cell_type %in% LINEAGES)) stop('Long input cell_type must use only the nine frozen lineage names.')
    if (!is.numeric(dat$abundance) || any(!is.finite(dat$abundance)) || any(dat$abundance<0)) stop('abundance must be finite, numeric and non-negative.')
    x <- matrix(0, length(ids), length(LINEAGES), dimnames=list(ids,LINEAGES))
    # Aggregate by keys, never by the incidental order of a metadata table.
    z <- aggregate(dat$abundance, list(sample_id=dat$sample_id,cell_type=dat$cell_type),sum)
    x[cbind(match(z$sample_id,ids), match(z$cell_type,LINEAGES))] <- z$x
  } else stop('Provide nine lineage columns, or the long columns sample_id, cell_type, abundance.')
  list(ids=ids,x=numeric_matrix(x),healthy=flags)
}
calculate_rnhi <- function(dat, default_health, custom_health=NULL, frac=.65) {
  frac <- validate_frac(frac)
  q <- parse_table(dat)
  composition <- closure(q$x)
  if (!is.null(custom_health) && any(q$healthy)) stop('Use either embedded healthy_reference rows or a separate healthy reference, not both.')
  health <- if (any(q$healthy)) q$x[q$healthy,,drop=FALSE] else if (!is.null(custom_health)) parse_table(custom_health)$x else default_health
  health <- numeric_matrix(health)
  center <- center_from(health,frac)
  positive <- composition
  # Same frozen per_group rule as Figure 2: healthy and patient DL separately.
  for (mask in list(q$healthy, !q$healthy)) {
    if (any(mask)) positive[mask,] <- multRepl(composition[mask,,drop=FALSE],frac)
  }
  delta <- sweep(clr(positive),2,center,'-')
  distance <- sqrt(rowSums(delta^2))
  out <- data.frame(sample_id=q$ids,healthy_reference=as.integer(q$healthy),RNHI=-distance,
                    AitchisonDistance=distance,selected_lineage_total=rowSums(q$x),
                    detected_lineages=rowSums(q$x>0),stringsAsFactors=FALSE)
  out$quality_flag <- ifelse(out$detected_lineages==1,'single_lineage_degenerate','ok')
  out <- cbind(out,as.data.frame(composition))
  list(data=out,composition=composition,center=center,n=nrow(out),nhealth=nrow(health),frac=frac)
}
distribution_plot <- function(r) {
  d <- r$data
  d$population <- ifelse(d$healthy_reference==1,'Healthy reference','AML / query')
  p <- ggplot2::ggplot(d,ggplot2::aes(RNHI,fill=population,color=population))
  if (nrow(d)>1 && length(unique(d$RNHI))>1)
    p <- p+ggplot2::geom_density(alpha=.18,linewidth=.9)
  p + ggplot2::geom_rug(alpha=.2)+ggplot2::scale_color_manual(values=c('AML / query'='#2454a6','Healthy reference'='#008574'))+
    ggplot2::scale_fill_manual(values=c('AML / query'='#2454a6','Healthy reference'='#008574'))+
    ggplot2::theme_classic(base_size=13)+ggplot2::theme(legend.position='top')+
    ggplot2::labs(x='RNHI (negative Aitchison distance)',y='Density',fill=NULL,color=NULL)
}
composition_plot <- function(r) {
  ord <- order(r$data$RNHI,decreasing=TRUE)
  x <- r$composition[ord,,drop=FALSE]
  z <- data.frame(rank=rep(seq_len(nrow(x)),ncol(x)),lineage=factor(rep(colnames(x),each=nrow(x)),levels=rev(LINEAGES)),proportion=as.vector(x))
  ggplot2::ggplot(z,ggplot2::aes(rank,lineage,fill=proportion))+ggplot2::geom_raster()+
    ggplot2::scale_fill_gradientn(colours=c('#F7FBFF','#DEEBF7','#9ECAE1','#4292C6','#08519C'),limits=c(0,1),oob=scales::squish,labels=scales::percent)+
    ggplot2::theme_classic(base_size=13)+ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank(),legend.position='right')+
    ggplot2::labs(x='Samples ordered by RNHI (high to low)',y=NULL,fill='Within-sample\nproportion')
}

make_sc_umap <- function(r, dat=NULL, n_per_sample=60) {
  # If uploaded data already contain coordinates, preserve them; otherwise create a
  # deterministic composition-space embedding with cell-level jitter for a useful demo.
  score <- r$data[,c('sample_id','RNHI','healthy_reference'),drop=FALSE]
  plot_dat <- dat
  if (!is.null(plot_dat) && !'sample_id' %in% names(plot_dat) && 'spot_id' %in% names(plot_dat)) plot_dat$sample_id <- plot_dat$spot_id
  if (!is.null(plot_dat) && !'sample_id' %in% names(plot_dat) && 'orig.ident' %in% names(plot_dat)) plot_dat$sample_id <- plot_dat$orig.ident
  if (!is.null(plot_dat) && !'x' %in% names(plot_dat) && 'pxl_col_in_fullres' %in% names(plot_dat)) plot_dat$x <- plot_dat$pxl_col_in_fullres
  if (!is.null(plot_dat) && !'y' %in% names(plot_dat) && 'pxl_row_in_fullres' %in% names(plot_dat)) plot_dat$y <- plot_dat$pxl_row_in_fullres
  if (!is.null(plot_dat) && !'umap_1' %in% names(plot_dat) && 'UMAP_1' %in% names(plot_dat)) plot_dat$umap_1 <- plot_dat$UMAP_1
  if (!is.null(plot_dat) && !'umap_2' %in% names(plot_dat) && 'UMAP_2' %in% names(plot_dat)) plot_dat$umap_2 <- plot_dat$UMAP_2
  if (!is.null(plot_dat) && all(c('umap_1','umap_2','sample_id') %in% names(plot_dat))) {
    z <- plot_dat[,c('sample_id','umap_1','umap_2'),drop=FALSE]
    z$umap_1 <- as.numeric(z$umap_1); z$umap_2 <- as.numeric(z$umap_2)
    if (any(!is.finite(z$umap_1)) || any(!is.finite(z$umap_2)))
      stop('UMAP coordinates must be finite numeric values (umap_1/umap_2).')
    z <- merge(z,score,by='sample_id',all.x=TRUE,sort=FALSE)
    z$cell_index <- seq_len(nrow(z)); return(z)
  }
  set.seed(20260913)
  X <- log(r$composition + 1e-6)
  if (nrow(X) >= 2) {
    # prcomp(scale.=TRUE) rejects a zero-variance lineage; composition previews
    # remain deterministic and useful when a user uploads a constant axis.
    sds <- apply(X,2,stats::sd)
    pc <- stats::prcomp(X,center=TRUE,scale.=all(is.finite(sds) & sds > 0))$x
    centers <- pc[,seq_len(min(2,ncol(pc))),drop=FALSE]
  } else centers <- matrix(c(0,0),nrow=1,ncol=2)
  if (ncol(centers)<2) centers <- cbind(centers,0)
  z <- do.call(rbind,lapply(seq_len(nrow(r$data)),function(i) {
    prop <- r$composition[i,]; cells <- sample(LINEAGES,n_per_sample,replace=TRUE,prob=prop)
    data.frame(sample_id=r$data$sample_id[i],cell_type=cells,umap_1=centers[i,1]+rnorm(n_per_sample,sd=.22),umap_2=centers[i,2]+rnorm(n_per_sample,sd=.22),stringsAsFactors=FALSE)
  }))
  merge(z,score,by='sample_id',all.x=TRUE,sort=FALSE)
}
sc_umap_plot <- function(r) {
  z <- r$plot_data
  ggplot2::ggplot(z,ggplot2::aes(umap_1,umap_2,color=RNHI))+ggplot2::geom_point(size=.9,alpha=.62)+
    ggplot2::scale_color_gradientn(colours=c('#A8C7DD','#C8DCEA','#F0E0C2','#F2B18F','#D36B67'),name='RNHI')+
    ggplot2::coord_equal()+ggplot2::theme_classic(base_size=13)+ggplot2::labs(x='UMAP 1',y='UMAP 2',title='scRNA-seq cells coloured by sample RNHI')+
    ggplot2::theme(plot.background=ggplot2::element_rect(fill='#faf8f3',colour=NA),panel.background=ggplot2::element_rect(fill='#faf8f3',colour=NA))
}
make_spot_umap <- function(r) {
  X <- log(r$composition + 1e-6)
  if (nrow(X) >= 2) {
    sds <- apply(X,2,stats::sd)
    pc <- stats::prcomp(X,center=TRUE,scale.=all(is.finite(sds) & sds > 0))$x
    xy <- pc[,seq_len(min(2,ncol(pc))),drop=FALSE]
  } else xy <- matrix(c(0,0),nrow=1,ncol=2)
  if (ncol(xy)<2) xy <- cbind(xy,0)
  data.frame(sample_id=r$data$sample_id,umap_1=xy[,1],umap_2=xy[,2],RNHI=r$data$RNHI,stringsAsFactors=FALSE)
}
spot_umap_plot <- function(r) {
  z <- r$plot_data_umap
  ggplot2::ggplot(z,ggplot2::aes(umap_1,umap_2,color=RNHI))+ggplot2::geom_point(size=3.6,alpha=.86)+
    ggplot2::scale_color_gradientn(colours=c('#A8C7DD','#C8DCEA','#F0E0C2','#F2B18F','#D36B67'),name='RNHI')+
    ggplot2::coord_equal()+ggplot2::theme_classic(base_size=13)+ggplot2::labs(x='UMAP 1',y='UMAP 2',title='Visium spots — UMAP view')+
    ggplot2::theme(plot.background=ggplot2::element_rect(fill='#faf8f3',colour=NA),panel.background=ggplot2::element_rect(fill='#faf8f3',colour=NA))
}
make_spatial_data <- function(r, dat=NULL) {
  score <- r$data[,c('sample_id','RNHI','healthy_reference'),drop=FALSE]
  plot_dat <- dat
  if (!is.null(plot_dat) && !'sample_id' %in% names(plot_dat) && 'spot_id' %in% names(plot_dat)) plot_dat$sample_id <- plot_dat$spot_id
  if (!is.null(plot_dat) && !'sample_id' %in% names(plot_dat) && 'orig.ident' %in% names(plot_dat)) plot_dat$sample_id <- plot_dat$orig.ident
  image_coordinate_input <- !is.null(plot_dat) && ((all(c('imagecol','imagerow') %in% names(plot_dat))) || all(c('pxl_col_in_fullres','pxl_row_in_fullres') %in% names(plot_dat)))
  if (!is.null(plot_dat) && !'x' %in% names(plot_dat) && 'pxl_col_in_fullres' %in% names(plot_dat)) plot_dat$x <- plot_dat$pxl_col_in_fullres
  if (!is.null(plot_dat) && !'y' %in% names(plot_dat) && 'pxl_row_in_fullres' %in% names(plot_dat)) plot_dat$y <- plot_dat$pxl_row_in_fullres
  coord_names <- if (!is.null(plot_dat) && all(c('x','y') %in% names(plot_dat))) c('x','y') else if (!is.null(plot_dat) && all(c('imagecol','imagerow') %in% names(plot_dat))) c('imagecol','imagerow') else NULL
  if (!is.null(coord_names)) {
    z <- plot_dat[,c('sample_id',coord_names),drop=FALSE]; names(z)[2:3] <- c('x','y'); z <- z[!duplicated(z$sample_id),,drop=FALSE]
    z$coord_type <- if (image_coordinate_input) 'image' else 'cartesian'
    z$x <- as.numeric(z$x); z$y <- as.numeric(z$y)
    if (any(!is.finite(z$x)) || any(!is.finite(z$y)))
      stop('Spatial coordinates must be finite numeric values (x/y or imagecol/imagerow).')
    return(merge(z,score,by='sample_id',all.x=TRUE,sort=FALSE))
  }
  n <- nrow(r$data); side <- ceiling(sqrt(n)); grid <- expand.grid(x=seq_len(side),y=seq_len(side))[seq_len(n),]
  grid$coord_type <- 'cartesian'
  cbind(grid,score)
}
spatial_plot <- function(r) {
  z <- r$plot_data
  rng <- range(z$RNHI, finite=TRUE)
  if (!length(rng) || any(!is.finite(rng))) stop('Spatial RNHI map has no finite scores.')
  if (diff(rng) == 0) rng <- rng + c(-.5,.5)
  map_limits <- c(rng[1], max(0, rng[2]))
  if (diff(map_limits) == 0) map_limits <- map_limits + c(-.5,.5)
  # A restrained editorial palette keeps the map legible in print and on screen:
  # soft blue = more negative imbalance, terracotta = closer to the healthy centre.
  pal <- c('#A8C7DD','#C8DCEA','#F0E0C2','#F2B18F','#D36B67')
  span <- max(diff(range(z$x,finite=TRUE)),diff(range(z$y,finite=TRUE)),na.rm=TRUE)
  spot_size <- if (is.finite(span) && span <= 20) 8.2 else 5.4
  y_scale <- if ('coord_type' %in% names(z) && any(z$coord_type == 'image')) ggplot2::scale_y_reverse(expand=ggplot2::expansion(mult=.13)) else ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=.13))
  ggplot2::ggplot(z,ggplot2::aes(x,y,fill=RNHI))+
    ggplot2::geom_point(shape=21,size=spot_size,alpha=.96,colour='#fffdf8',stroke=.32)+
    ggplot2::scale_fill_gradientn(colours=pal,limits=map_limits,oob=scales::squish,name='RNHI\n(0 = healthy-like)',guide=ggplot2::guide_colorbar(frame.colour='#d6d0c5',ticks.colour='#8b8379',barheight=grid::unit(62,'mm'),barwidth=grid::unit(5,'mm')))+
  ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.13)) + y_scale +
    ggplot2::coord_equal()+ggplot2::theme_minimal(base_size=13)+
    ggplot2::labs(title='Visium spot spatial RNHI map',subtitle='Each spot is coloured by its compositional distance from the healthy centre',x=NULL,y=NULL)+
    ggplot2::theme(plot.background=ggplot2::element_rect(fill='#faf8f3',colour=NA),panel.background=ggplot2::element_rect(fill='#faf8f3',colour=NA),panel.grid=ggplot2::element_blank(),axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),plot.title=ggplot2::element_text(hjust=.5,face='bold',colour='#334947',margin=ggplot2::margin(b=3)),plot.subtitle=ggplot2::element_text(hjust=.5,colour='#6a756d',size=10,margin=ggplot2::margin(b=10)),legend.position='right',legend.title=ggplot2::element_text(colour='#455b56',face='bold'),legend.text=ggplot2::element_text(colour='#5b6862'))
}
