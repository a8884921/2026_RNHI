# Run the bundled MuSiC path on the real BeatAML2 matrix (first eight samples).
# The subset keeps this reproducibility check quick while using the project file.
source('rnhi_core.R', local=TRUE)
src <- file.path('..','02_Figure2','01_RNHI构建与定义','A_输入数据','02_真实bulk_12队列','bulk','BeatAML2.csv')
if (!file.exists(src)) stop('Real BeatAML2 file not found at ', normalizePath(src, mustWork=FALSE))
full <- read_input(src)
dat <- full[, c(names(full)[1], names(full)[seq_len(min(8L, ncol(full)-1L)) + 1L]), drop=FALSE]
mu <- music_deconvolve(dat, iter.max=300, eps=.02)
ref <- parse_table(read_input('composition_A4.csv'))
r <- calculate_rnhi(mu$data, ref$x[ref$healthy,,drop=FALSE])
stopifnot(nrow(mu$data)==8L, all(mu$qc$status=='success'), nrow(r$data)==8L, all(is.finite(r$data$RNHI)))
bad <- dat
bad[[2]] <- 0
bad_error <- tryCatch(music_deconvolve(bad, iter.max=50, eps=.02), error=function(e) e)
stopifnot(inherits(bad_error, 'error'), grepl('no partial RNHI result', conditionMessage(bad_error), fixed=TRUE))
cat('REAL BEATAML2 MUSIC PASS —',nrow(r$data),'samples;',mu$common_genes,'common genes; median RNHI',median(r$data$RNHI),'\n')
