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

options(warn=-1)
suppressMessages({library(ggplot2); library(ggpubr); library(png); library(magrittr)})

FROZEN <- "data/figure_2/Figure2C_data_frozenUMAP_withRNHI.csv"
ANNOT <- "data/figure_2/annot_master_keep7.csv"
OUTD <- "output/figure_2/molecular_umap"
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

um <- data.table::fread(FROZEN, data.table=FALSE, check.names=FALSE, encoding="UTF-8")
an <- data.table::fread(ANNOT,  data.table=FALSE, check.names=FALSE, encoding="UTF-8")
um$sample_id <- anonymous_ids(um$sample_id)
an$sample_id <- anonymous_ids(an$sample_id)
an1 <- an[!duplicated(an$sample_id), c("sample_id","TP53","NPM1","FLT3_ITD","IDH1_2","ELN2017")]
d <- merge(um[,c("sample_id","UMAP1","UMAP2")], an1, by="sample_id", all.x=TRUE)
stopifnot(nrow(d)==nrow(um))

na_lab <- "NA"
lab <- list(TP53=c("Mut"="TP53-mutated","WT"="TP53-WT"),
            NPM1=c("Mut"="NPM1-mutated","WT"="NPM1-WT"),
            FLT3_ITD=c("Positive"="FLT3-ITD+","Negative"="FLT3-ITD-"),
            IDH1_2=c("Mut"="IDH1/2-mutated","WT"="IDH1/2-WT"),
            ELN2017=c("Favorable"="ELN-favorable","Intermediate"="ELN-intermediate","Adverse"="ELN-adverse"))
to_fct <- function(x,lv){ x[is.na(x)|x==""|x=="NA"|x=="nan"]<-na_lab; factor(x,levels=c(lv,na_lab)) }
for(v in names(lab)) d[[v]] <- to_fct(d[[v]], names(lab[[v]]))

xmid <- mean(range(d$UMAP1)); ymid <- mean(range(d$UMAP2))
span <- max(diff(range(d$UMAP1)), diff(range(d$UMAP2)))*1.03
XLIM <- c(xmid-span/2, xmid+span/2); YLIM <- c(ymid-span/2, ymid+span/2)

base_panel <- function(var){
  lv <- c(names(lab[[var]]), na_lab)
  cols <- setNames(c("#D55E00","#0072B2","#009E73"), setdiff(lv,na_lab)); cols[na_lab]<-"#C7C7C7"; cols<-cols[lv]
  ntab<-table(d[[var]])
  labs_txt<-sapply(lv, function(z) sprintf("%s (n=%d)", if (is.na(lab[[var]][z])) na_lab else lab[[var]][z], as.integer(ntab[z]))); names(labs_txt)<-lv
  ggplot(d, aes(UMAP1, UMAP2, colour=.data[[var]])) +
    geom_point(size=1.4, alpha=0.9, stroke=0) +
    scale_colour_manual(values=cols, labels=labs_txt, breaks=lv, drop=FALSE, name=NULL,
                        guide=guide_legend(override.aes=list(size=3.2, alpha=1), ncol=1)) +
    labs(x="UMAP1", y="UMAP2") +
    coord_equal(xlim=XLIM, ylim=YLIM) +
    theme_bw(base_size=16) +
    theme(text=element_text(colour="#1F1F1F"),
          axis.title=element_text(size=16), axis.text=element_text(size=12.5),
          axis.ticks=element_line(colour="#1F1F1F", linewidth=0.5),
          panel.border=element_rect(colour="#1F1F1F", fill=NA, linewidth=0.7),
          panel.grid=element_blank(),
          legend.position="right", legend.text=element_text(size=11),
          legend.key.size=unit(0.30,"cm"), legend.key=element_blank(),
          legend.background=element_blank(), legend.margin=margin(0,0,0,4),
          plot.margin=margin(2,2,1,2))
}
pl <- lapply(names(lab), base_panel)

for(i in seq_along(pl)){
  v<-names(lab)[i]
  ggsave(file.path(OUTD, sprintf("Figure2D_UMAP_overlay_%s.pdf", v)), pl[[i]], width=5.6, height=5.6, device=cairo_pdf, dpi=300)
  ggsave(file.path(OUTD, sprintf("Figure2D_UMAP_overlay_%s.png", v)), pl[[i]], width=5.6, height=5.6, dpi=300)
}

comb <- ggarrange(plotlist=pl, ncol=3, nrow=2, labels=LETTERS[1:5],
                  font.label=list(size=20, face="bold"), align="hv",
                  widths=rep(1,3), heights=c(1,1)) %>%
  annotate_figure(top=text_grob("Known AML molecular features on the frozen transcriptomic landscape",
                                face="bold", size=22),
                  bottom=text_grob("Same frozen Figure 2B UMAP (2,575 patients) in every panel; grey = not reported (NA); coordinates never re-computed; PCA/UMAP not re-run",
                                   color="grey25", size=12.5))
ggsave(file.path(OUTD,"Figure2D_main_overlay_composite_final.pdf"), comb, width=16.8, height=12.2, device=cairo_pdf, dpi=300, limitsize=FALSE)
ggsave(file.path(OUTD,"Figure2D_main_overlay_composite_final.png"), comb, width=16.8, height=12.2, dpi=300, limitsize=FALSE)
img <- png::readPNG(file.path(OUTD,"Figure2D_main_overlay_composite_final.png"))
png::writePNG(img, file.path(OUTD,"Figure2D_main_overlay_composite_final_preview.png"))
cat("Output completed.\n")
