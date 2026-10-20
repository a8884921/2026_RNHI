# Server-side regression check for the same path used by the browser upload.
file_arg <- commandArgs(FALSE)[grep('^--file=', commandArgs(FALSE))]
file_arg <- if(length(file_arg)) sub('^--file=', '', file_arg[[1]]) else 'tests/local_real_data_test.R'
setwd(normalizePath(file.path(dirname(file_arg), '..')))
source('rnhi_core.R', local=TRUE)
dat <- read_input('composition_A4.csv')
parsed <- parse_table(dat)
calc <- calculate_rnhi(dat, parsed$x[parsed$healthy,,drop=FALSE])
formal <- read.csv('../02_Figure2/01_RNHI构建与定义/复现包_RNHI_A4_20260907/03_输出/a4_rnhi_scores.csv', check.names=FALSE)
keep <- formal[,c('sample_id','RNHI')]
joined <- merge(calc$data[calc$data$healthy_reference==0,c('sample_id','RNHI')], keep, by='sample_id', suffixes=c('_shiny','_formal'))
delta <- joined$RNHI_shiny - joined$RNHI_formal
report <- data.frame(n_samples=nrow(joined), max_abs_delta=max(abs(delta)), exact_fraction=mean(abs(delta)<1e-10), median_rnhi=median(joined$RNHI_shiny))
write.csv(report, 'data/local_validation_report.csv', row.names=FALSE)
stopifnot(nrow(joined)==3310, max(abs(delta)) < 1e-10, mean(abs(delta)<1e-10)==1)
cat('LOCAL REAL DATA PASS\n'); print(report)
