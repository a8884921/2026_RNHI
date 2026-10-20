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

suppressPackageStartupMessages({library(ggplot2); library(data.table); library(patchwork); library(ggpubr); library(gtable); library(grid)})

.fmt2 <- function(x) formatC(x, format="f", digits=0)

healthy_ellipses <- function(hfull, level=0.90, npts=100) {
  lvs <- unique(hfull$state12)
  out <- list()
  qr  <- sqrt(stats::qchisq(level, df=2))
  for (st in lvs) {
    d <- hfull[state12==st]
    X <- as.matrix(d[, .(u1,u2)])
    if (nrow(X) < 5) next
    mu <- colMeans(X); S <- stats::cov(X)
    e  <- eigen(S, symmetric=TRUE)
    a  <- qr*sqrt(max(e$values, 1e-9)); b <- qr*sqrt(max(e$values[2], 1e-9))
    th <- atan2(e$vectors[2,1], e$vectors[1,1])
    tt <- seq(0, 2*pi, length.out=npts)
    x <- mu[1] + a*cos(th)*cos(tt) - b*sin(th)*sin(tt)
    y <- mu[2] + a*sin(th)*cos(tt) + b*cos(th)*sin(tt)
    out[[length(out)+1]] <- data.frame(state12=st, x=x, y=y)
  }
  rbindlist(out)
}

draw_f6a <- function(sub, hfull=NULL, out_dir, wA=4.5, hA=4.5) {
  cols12 <- c("Normal_HSC_MPP"="#E64B35","Normal_LMPP"="#F39B7F","Normal_MyeloidProgenitor"="#4DBBD5",
   "Normal_Monocyte"="#00A087","Normal_Granulocytic"="#3C5488","Normal_Megakaryocytic"="#7E6148",
   "Normal_EarlyErythroid"="#B2182B","Normal_LateErythroid"="#D6604D","Normal_B"="#4393C3",
   "Normal_T"="#92C5DE","Normal_NK"="#F4A582","Normal_Dendritic"="#FDDBC7")
  stateLv <- names(cols12)
  lab12 <- c("Normal_HSC_MPP"="HSC/MPP","Normal_LMPP"="LMPP","Normal_MyeloidProgenitor"="Myeloid prog.",
   "Normal_Monocyte"="Monocyte","Normal_Granulocytic"="Granulocytic","Normal_Megakaryocytic"="Megakaryocytic",
   "Normal_EarlyErythroid"="Early erythroid","Normal_LateErythroid"="Late erythroid","Normal_B"="B",
   "Normal_T"="T","Normal_NK"="NK","Normal_Dendritic"="Dendritic")
  sub$st_lab <- factor(lab12[as.character(sub$state12)], levels=lab12[stateLv])

  n_H <- sum(sub$grp3=="Healthy"); n_hi <- sum(sub$grp3=="AML_RNHI-high"); n_lo <- sum(sub$grp3=="AML_RNHI-low")

  pA <- ggplot() +
    geom_point(data=sub, aes(u1,u2,colour=state12), size=0.12, stroke=0, alpha=0.7) +
    scale_colour_manual(values=cols12, name=NULL, breaks=stateLv, labels=lab12[stateLv]) +
    guides(colour=guide_legend(override.aes=list(size=2.4, alpha=1), ncol=4)) +
    labs(x="BoneMarrowMap UMAP 1", y="BoneMarrowMap UMAP 2") +
    theme_minimal(base_size=13) +
    theme(panel.grid=element_blank(), axis.line=element_line(colour="grey40", linewidth=0.4),
          axis.text=element_blank(), axis.ticks=element_blank(),
          legend.position="bottom", legend.text=element_text(size=8.5),
          legend.key.size=unit(0.34,"cm"), legend.key.spacing.x=unit(0.14,"cm"),
          legend.margin=margin(0,0,0,0), plot.margin=margin(4,6,2,4)) +
    ggtitle(sprintf("Normal-like haematopoietic cells\n12 states (n = %s)", .fmt2(nrow(sub)))) +
    theme(plot.title=element_text(size=11.5, face="bold", hjust=0.5, lineheight=1.05))

  grp_cols <- c("Healthy"="#009E73", "AML_RNHI-high"="#D55E00", "AML_RNHI-low"="#CC79A7")
  grp_lab  <- c("Healthy"="Healthy\nn = 54 donors",
                "AML_RNHI-high"="RNHI-high\nn = 5 patients",
                "AML_RNHI-low"="RNHI-low\nn = 5 patients")
  sub$grp3f <- factor(sub$grp3, levels=names(grp_cols))
  pB <- ggplot() +
    geom_point(data=sub, aes(u1,u2,colour=grp3f), size=0.12, stroke=0, alpha=0.7) +
    scale_colour_manual(values=grp_cols, name=NULL, labels=grp_lab) +
    guides(colour=guide_legend(override.aes=list(size=2.6, alpha=1), ncol=3)) +
    labs(x="BoneMarrowMap UMAP 1", y="BoneMarrowMap UMAP 2") +
    theme_minimal(base_size=13) +
    theme(panel.grid=element_blank(), axis.line=element_line(colour="grey40", linewidth=0.4),
          axis.text=element_blank(), axis.ticks=element_blank(),
          legend.position="bottom", legend.text=element_text(size=8.6, lineheight=0.95),
          legend.key.size=unit(0.30,"cm"), legend.key.spacing.x=unit(0.16,"cm"),
          legend.margin=margin(2,4,0,4), plot.margin=margin(4,6,2,4)) +
    ggtitle("Normal-like cells by source\nHealthy donors vs AML RNHI-high/low") +
    theme(plot.title=element_text(size=11.5, face="bold", hjust=0.5, lineheight=1.05))

  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  emit <- function(p, base, w, h) {
    pdf_f<-file.path(out_dir,paste0(base,".pdf")); png_f<-file.path(out_dir,paste0(base,".png"))
    grDevices::cairo_pdf(pdf_f,width=w,height=h,family="Helvetica"); print(p); grDevices::dev.off()
    grDevices::png(png_f,width=w,height=h,units="in",res=300,type="cairo"); print(p); grDevices::dev.off()
    cat("Output completed.\n")
  }
  emit(pA,"Figure6A_panelA_12states_UMAP",wA,hA)
  emit(pB,"Figure6A_panelB_threegroup_UMAP",wA,hA)

  comp<-pA + pB + plot_layout(ncol=2, widths=c(1,1)) +
    plot_annotation(title="Residual normal haematopoiesis in Healthy vs RNHI-high/low AML (BoneMarrowMap)",
                    subtitle="Healthy = green (45 ref + 9 ind donors); RNHI-high = orange; RNHI-low = purple. Same BoneMarrowMap coordinates across all cells.",
                    tag_levels="A", theme=theme(plot.title=element_text(size=13,face="bold",hjust=0.5),
                                                plot.subtitle=element_text(size=9,hjust=0.5,colour="grey25")))
  pdf_f<-file.path(out_dir,"Figure6A_composite.pdf"); png_f<-file.path(out_dir,"Figure6A_composite.png")
  grDevices::cairo_pdf(pdf_f,width=10.2,height=5.2,family="Helvetica"); print(comp); grDevices::dev.off()
  grDevices::png(png_f,width=10.2,height=5.2,units="in",res=300,type="cairo"); print(comp); grDevices::dev.off()
  cat("Output completed.\n")
  invisible(list(A=pA,B=pB))
}
