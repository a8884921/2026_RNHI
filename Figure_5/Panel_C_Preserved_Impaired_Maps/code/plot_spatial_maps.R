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

suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table); library(ComplexHeatmap); library(circlize)})
options(warn=-1)
W <- "output/figure_5/fragmentation"
O <- "output/figure_5/fragmentation"
dir.create(O, recursive = TRUE, showWarnings = FALSE)
dir.create(O, showWarnings=FALSE, recursive=TRUE)
ann <- fread(file.path(W,"preserved_impaired_annotated.csv"))
ann$code <- anonymous_ids(ann$code)
fr  <- fread(file.path(W,"per_sample_fragmentation_final.csv"))
adj <- fread(file.path(W,"aml_niche_adjacency_preservation.csv"))
prof<- as.matrix(fread(file.path(W,"healthy_niche_adjacency_profile.csv")), rownames=1)
label <- function(s) anonymous_ids(s)
pal_state <- c(preserved="#2166AC", impaired="#B2182B")
journal <- theme_void(base_size=12)+theme(plot.title=element_text(face="bold",size=13,hjust=.5,margin=margin(0,0,4,0)),
  plot.subtitle=element_text(size=9.5,hjust=.5,colour="grey25"),legend.title=element_text(face="bold",size=10),
  legend.text=element_text(size=9),legend.key.height=unit(.4,"cm"),legend.key.width=unit(.35,"cm"),plot.margin=margin(4,4,4,4))

pl <- list()
for(g in unique(ann$code)){
  d <- ann[code==g]
  p <- ggplot(d, aes(x,y,colour=state)) + geom_point(size=1.6,stroke=0,alpha=.85) +
    scale_colour_manual(values=pal_state,name="Normal haematopoiesis",drop=FALSE) +
    scale_y_reverse()+coord_fixed()+journal+
    labs(title=label(g), subtitle=paste0("n=",nrow(d)," | preserved=",round(100*mean(d$state=="preserved"),1),"%"))
  pl[[g]] <- p
}
g1 <- wrap_plots(pl, ncol=4)
ggsave(file.path(O,"FigD1_preserved_impaired_maps.pdf"), g1, width=4.6*4, height=4.6*2, device=cairo_pdf, limitsize=FALSE)
ggsave(file.path(O,"FigD1_preserved_impaired_maps.png"), g1, width=4.6*4, height=4.6*2, dpi=300, device=png, type="cairo", limitsize=FALSE)

fr$lab <- vapply(fr$sample_id,label,character(1))
fr$grp_col <- ifelse(fr$group=="Healthy","#00A087","#E64B35")
mkbar <- function(var, ylab){
  fr$v <- fr[[var]]
  ggplot(fr, aes(reorder(lab,-v), v, fill=group)) + geom_col(width=.7) +
    scale_fill_manual(values=c(Healthy="#00A087",AML="#E64B35")) +
    coord_flip() + theme_minimal(base_size=12) +
    labs(x="",y=ylab) + theme(legend.position="none")
}
m1 <- mkbar("preserved_area_fraction","Preserved area fraction")
m2 <- mkbar("largest_cc_fraction","Largest component / preserved")
m3 <- mkbar("mean_component_size","Mean island size")
m4 <- mkbar("boundary_density","Boundary density")
g2 <- (m1+m2)/(m3+m4)
ggsave(file.path(O,"FigD2_fragmentation_metrics.pdf"), g2, width=9, height=9, device=cairo_pdf)
ggsave(file.path(O,"FigD2_fragmentation_metrics.png"), g2, width=9, height=9, dpi=300, device=png, type="cairo")

M <- as.matrix(prof); colnames(M)<-rownames(M)
pdf(file.path(O,"FigD3_healthy_niche_adjacency.pdf"), width=4.5,height=4.5)
h <- Heatmap(M, name="norm adjacency", col=colorRamp2(c(0,.35,.7),c("white","#FDD49E","#B30000")),
             cluster_rows=FALSE,cluster_columns=FALSE,row_names_gp=gpar(fontsize=8),column_names_gp=gpar(fontsize=7,rot=40),border=TRUE)
draw(h, padding=unit(c(2,2,14,2),"mm")); dev.off()
png(file.path(O,"FigD3_healthy_niche_adjacency.png"), width=4.5,height=4.5,units="in",res=300)
draw(h, padding=unit(c(2,2,14,2),"mm")); dev.off()

write.csv(adj, file.path(O,"figD3_adjacency_preservation.csv"), row.names=FALSE)
cat("FIGD DONE\n")
