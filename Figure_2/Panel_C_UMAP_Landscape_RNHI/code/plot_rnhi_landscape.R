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

options(warn = -1)
suppressMessages({library(ggplot2)})

U <- "data/figure_2/11_UMAP_coordinates.csv"
P <- "data/figure_2/fig2b_keep7_patients.csv"
OUTD <- "output/figure_2/rnhi_umap"

um <- read.csv(U, check.names=FALSE, fileEncoding="UTF-8-BOM")
pt <- read.csv(P, check.names=FALSE, fileEncoding="UTF-8-BOM")[, c("sample_id","RNHI")]
um$sample_id <- anonymous_ids(um$sample_id)
pt$sample_id <- anonymous_ids(pt$sample_id)
d <- merge(um[, c("sample_id","UMAP1","UMAP2")], pt, by="sample_id", all.x=TRUE)
d$is_na <- is.na(d$RNHI); d <- d[order(d$is_na), ]
rng <- range(d$RNHI, na.rm=TRUE); n <- nrow(d)
cat("mapped:",n,"NA:",sum(d$is_na),"range:",round(rng,2),"\n")

p <- ggplot(d, aes(UMAP1, UMAP2)) +
  geom_point(data=subset(d,is_na), colour="#C7C7C7", size=1.5, alpha=0.9) +
  geom_point(data=subset(d,!is_na), aes(colour=RNHI), size=1.7, alpha=1.0, stroke=0) +
  scale_colour_gradientn(colours=c("#0D0887","#5B02A3","#9C179E","#CC4778","#ED7953","#FDB42F","#F0F921"),
                         limits=rng, breaks=seq(floor(rng[1]),ceiling(rng[2]),by=2),
                         labels=function(x)sprintf("%.0f",x), name="RNHI") +
  guides(colour=guide_colourbar(barwidth=1.1,barheight=6.2,title.position="top",
                                title.hjust=0.5,ticks.colour="black",frame.colour="black")) +
  labs(x="UMAP1", y="UMAP2",
       title="RNHI in the AML landscape",
       subtitle=sprintf("n = %d | 7 cohorts | frozen UMAP", n),
       caption=NULL) +
  theme_bw(base_size=19) +
  theme(text=element_text(colour="#1F1F1F"),
        plot.title=element_text(size=23, face="bold", hjust=0.5, margin=margin(b=6)),
        plot.subtitle=element_text(size=15, colour="#4A4A4A", hjust=0.5, margin=margin(b=4)),
        plot.caption=element_text(size=10.5, colour="#4A4A4A", hjust=0.5, margin=margin(t=6)),
        axis.title=element_text(size=20), axis.text=element_text(size=17),
        axis.ticks=element_line(colour="#1F1F1F", linewidth=0.5),
        panel.border=element_rect(colour="#1F1F1F", fill=NA, linewidth=0.7),
        panel.grid.major=element_line(colour="#E6EAEE", linewidth=0.3),
        panel.grid.minor=element_blank(), panel.background=element_rect(fill="white"),
        plot.background=element_rect(fill="white"),
        legend.position="right",
        legend.title=element_text(size=17, face="bold"),
        legend.text=element_text(size=15), legend.key=element_blank(),
        legend.background=element_blank(), plot.margin=margin(10,10,8,10))

ggsave(file.path(OUTD,"Figure2C_RNHI_plasma_main.pdf"), p, width=6.4, height=5.2, device=cairo_pdf, dpi=300)
ggsave(file.path(OUTD,"Figure2C_RNHI_plasma_main.png"), p, width=6.4, height=5.2, dpi=300)
cat("saved plasma main map\n")
