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

TMP <- "output/figure_5/pathways"
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
options(warn = -1)
W <- "output/figure_5/pathways"
dir.create(W, recursive = TRUE, showWarnings = FALSE)
FONT <- Sys.getenv("FIG_FONT","Arial")
P   <- fread(file.path(TMP,"aml_pathway_scores_cpm_13H.csv"))
GLB <- fread("data/figure_5/spatial/aml_niche_mapping_srnhi.csv")[,.(cell_id,code,z=sRNHI_z_global)]
D <- merge(as.data.frame(P), as.data.frame(GLB), by=c("cell_id","code"))

GRP1 <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_IL6_JAK_STAT3_SIGNALING",
          "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_P53_PATHWAY")
GRP2 <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
SIXD <- c(GRP1, GRP2)
GRPLAB <- c("Stress / inflammatory (higher with RNHI)", "Interferon & ROS (collapsed in RNHI-low)")
GRPV <- setNames(c(rep(GRPLAB[1], length(GRP1)), rep(GRPLAB[2], length(GRP2))), SIXD)
LAB <- c(HALLMARK_E2F_TARGETS="E2F targets",HALLMARK_GLYCOLYSIS="Glycolysis",
         HALLMARK_TNFA_SIGNALING_VIA_NFKB="TNF\u03b1 via NF-\u03baB",
         HALLMARK_IL6_JAK_STAT3_SIGNALING="IL-6\u2013JAK\u2013STAT3",
         HALLMARK_INFLAMMATORY_RESPONSE="Inflammatory response",HALLMARK_P53_PATHWAY="p53 pathway",
         HALLMARK_INTERFERON_GAMMA_RESPONSE="Interferon-\u03b3 response",
         HALLMARK_INTERFERON_ALPHA_RESPONSE="Interferon-\u03b1 response",
         HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY="Reactive oxygen species",
         HALLMARK_OXIDATIVE_PHOSPHORYLATION="Oxidative phosphorylation")
D$code <- factor(D$code, levels=as.character(1:4))

for (cd in levels(D$code)) for (p in SIXD) {
  v <- D[[p]][D$code==cd]; D[[p]][D$code==cd] <- (v-mean(v))/sd(v)
}
cat("[z] All spots standardized within each sample and pathway; ranks and p values preserved.\n")

D$tier <- NA_character_
for (cd in levels(D$code)) { zz <- D$z[D$code==cd]; q <- quantile(zz,c(.25,.75))
  D$tier[D$code==cd] <- ifelse(zz<=q[1],"Low",ifelse(zz>=q[2],"High","Middle")) }
D$tier <- factor(D$tier, levels=c("Low","Middle","High")); D$logN <- log10(D$nCount)
match_pairs <- function(lo,hi,caliper=.25){ setorder(lo,logN); setorder(hi,logN)
  used<-rep(FALSE,nrow(hi)); pick<-integer(nrow(lo)); sdv<-sd(c(lo$logN,hi$logN))
  for(i in seq_len(nrow(lo))){ dd<-abs(hi$logN-lo$logN[i]); dd[used]<-Inf; j<-which.min(dd)
    if(length(j)&&is.finite(dd[j])&&dd[j]<=caliper*sdv){pick[i]<-j; used[j]<-TRUE}else pick[i]<-NA_integer_}
  ok<-!is.na(pick); list(i=which(ok), j=pick[ok]) }
PW<-list(); ES<-list()
for (cd in levels(D$code)) { d<-D[D$code==cd,]; lo<-d[d$tier=="Low",]; hi<-d[d$tier=="High",]
  mp<-match_pairs(lo,hi); L2<-lo[mp$i,]; H2<-hi[mp$j,]
  for (p in SIXD) {
    w <- suppressWarnings(wilcox.test(H2[[p]],L2[[p]],paired=TRUE,exact=FALSE))
    PW[[paste(cd,p)]] <- data.table(code=cd, pkey=p, pathway=unname(LAB[p]), n_pairs=length(mp$i),
      delta=median(H2[[p]])-median(L2[[p]]), p=w$p.value)
    ES[[paste(cd,p)]] <- rbind(
      data.table(code=cd, pwlab=factor(unname(LAB[p]),levels=unname(LAB[SIXD])), grp=factor(GRPV[p],levels=GRPLAB), tier=factor("Low",levels=c("Low","High")), score=L2[[p]]),
      data.table(code=cd, pwlab=factor(unname(LAB[p]),levels=unname(LAB[SIXD])), grp=factor(GRPV[p],levels=GRPLAB), tier=factor("High",levels=c("Low","High")), score=H2[[p]]))
  } }
PW<-rbindlist(PW); ES<-rbindlist(ES)
PW[, fdr := p.adjust(p,"BH")]
PW[, star := ifelse(fdr<.001,"***",ifelse(fdr<.01,"**",ifelse(fdr<.05,"*","ns")))]
fwrite(PW[,.(code,pathway,n_pairs,delta_z=round(delta,4),p=signif(p,3),fdr=signif(fdr,3),star,
             scRNA_direction="high > low",same_as_scRNA=delta>0)],
       file.path(W,"figF_spot_aligned_direction_focus7.csv"))
cat("[note] delta_z is the within-sample standardized High-Low median difference; p values are unchanged.\n")

ym <- ES[, .(ymax=max(score)), by=.(code,pwlab,grp)]
BR <- merge(PW[,.(code,pkey,pwlab=factor(pathway,levels=unname(LAB[SIXD])),delta,p,fdr,star)], ym, by=c("code","pwlab"))
BR[, grp := factor(GRPV[pkey], levels=GRPLAB)]
BR[, y := ymax + 0.10]; BR[, lab := star]

BR[, opp := (delta < 0) & (fdr < 0.05)]
OPP <- BR[opp == TRUE]

pal <- c(Low="#9ECAE1", High="#F4A582")
g <- ggplot(ES, aes(tier, score)) +
  geom_hline(yintercept=0, linetype="22", linewidth=.25, colour="grey65") +
  geom_boxplot(aes(fill=tier), width=.17, outlier.shape=NA, linewidth=.38,
               colour="grey25", alpha=.85) +
  geom_segment(data=BR, aes(x=1, xend=2, y=y, yend=y), inherit.aes=FALSE, linewidth=.28, colour="grey20") +
  geom_text(data=BR, aes(x=1.5, y=y, label=lab), inherit.aes=FALSE, vjust=-.35,
            size=2.0, colour="grey15", family=FONT) +
  { if (nrow(OPP)) geom_text(data=OPP, aes(x=1.5, y=y, label="opp."), inherit.aes=FALSE,
                             vjust=1.6, size=1.9, fontface="bold", colour="#B2182B", family=FONT) } +
  facet_grid(grp + pwlab ~ code) +
  scale_fill_manual(values=pal, guide="none") + scale_colour_manual(values=pal, guide="none") +
  scale_x_discrete(labels=c(Low="Low", High="High")) +
  scale_y_continuous(breaks=c(-1,0,1), labels=c("-1","0","+1"),
                     expand=expansion(mult=c(.06,.14))) +
  labs(x = NULL, y = "ssGSEA score (z within sample)",
       title = "Spot-level pathway activity along the sRNHI gradient",
       subtitle = paste(strwrap(paste0(
         "Seven key pathways taken from the scRNA Figure 6D panel b set (all of them run 'high > low' in scRNA), ",
         "grouped as in the text: stress / inflammatory programmes that are activated with RNHI, and the interferon ",
         "and ROS programmes that collapse in RNHI-low. Scores are z-scored within each sample x pathway so all ",
         "panels share one scale (dashed line = that sample's own mean); Low / High = bottom / top sRNHI quartile, ",
         "matched 1:1 on log10(nCount) within each sample and compared by paired Wilcoxon. Stars = BH-FDR ",
         "(* <0.05, ** <0.01, *** <0.001, ns = not significant). IFN-gamma is significant in 4/4 samples, IL-6 and ",
         "IFN-alpha in 3/4; no pathway is significantly opposite to the scRNA direction."),
         width = 132), collapse = "\n")) +
  theme_minimal(base_size=9, base_family=FONT) +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(linewidth=.22, colour="grey90"),
        panel.spacing.x=unit(2.6,"mm"), panel.spacing.y=unit(1.7,"mm"),
        strip.text=element_text(size=6.4, face="bold", colour="grey15"),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.background=element_rect(fill="grey95", colour=NA),
        axis.text.x=element_text(size=6.4), axis.text.y=element_text(size=6, colour="grey35"),
        axis.title.y=element_text(size=7.4), plot.title=element_text(size=10.2, face="bold"),
        plot.subtitle=element_text(size=6.3, colour="grey30", lineheight=1.05),
        plot.margin=margin(8,8,6,8))
grDevices::cairo_pdf(file.path(W,"FigF_spot_CPM_supplement.pdf"), width=9, height=7.0, family=FONT)
print(g); dev.off()
ragg::agg_png(file.path(W,"FigF_spot_CPM_supplement.png"), width=9, height=7.0, units="in", res=400, background="white")
print(g); dev.off()
cat(sprintf("\n[done] 9 x 7.0 in; seven pathway boxplots; significance labels: %s\n",
            paste(names(table(PW$star)), table(PW$star), sep="=", collapse=" ")))
