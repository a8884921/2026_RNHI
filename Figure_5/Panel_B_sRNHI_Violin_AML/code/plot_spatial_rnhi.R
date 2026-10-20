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

MM      <- function(mm) mm / 25.4
W_FULL  <- MM(183)
Z_LIM   <- 3.0

SPOT_UM <- 100
BAR_UM  <- 1000
ALPHA   <- 0.05
EXPORT_PANELS <- TRUE

FONT <- local({
  cand <- c("Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "sans")
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    for (f in cand) {
      hit <- try(systemfonts::match_font(f), silent = TRUE)
      if (!inherits(hit, "try-error") && !is.null(hit$path) &&
          grepl(gsub("[^a-z]", "", tolower(f)),
                gsub("[^a-z]", "", tolower(basename(hit$path))), fixed = TRUE)) return(f)
    }
    return("sans")
  }
  cand[1]
})

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(grid); library(scales)
})

IN <- "data/figure_5/spatial"
OUT <- "output/figure_5/spatial"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(IN))

open_pdf <- function(fp, w, h) grDevices::cairo_pdf(fp, width = w, height = h, family = FONT)
open_png <- function(fp, w, h, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(fp, width = w, height = h, units = "in", res = dpi)
  else grDevices::png(fp, width = w, height = h, units = "in", res = dpi, type = "cairo")
}
quiet_mbcs <- function(expr) withCallingHandlers(expr, warning = function(w)
  if (grepl("mbcsToSbcs|conversion failure", conditionMessage(w))) invokeRestart("muffleWarning"))

need <- c("aml_niche_mapping_srnhi.csv", "healthy_heldout_validation.csv",
          "niche_healthy_srnhi_reference.csv", "content_validation.csv")
miss <- need[!file.exists(file.path(IN, need))]
if (length(miss)) stop("Missing input: ", paste(miss, collapse = ", "))

res <- fread(file.path(IN, "niche_healthy_srnhi_reference.csv"))
map <- fread(file.path(IN, "aml_niche_mapping_srnhi.csv"))
hld <- fread(file.path(IN, "healthy_heldout_validation.csv"))
val <- fread(file.path(IN, "content_validation.csv"))
map$sample_id <- anonymous_ids(map$sample_id)
hld$donor <- anonymous_ids(hld$donor)
val$sample <- anonymous_ids(val$sample)

SAMP <- as.character(1:4)

SAMP_LAB <- setNames(SAMP, SAMP)
DL <- function(x) unname(SAMP_LAB[as.character(x)])
stopifnot(all(SAMP %in% map$sample_id))
map[, sid := factor(sample_id, levels = SAMP)]

nmap <- unique(map[!is.na(assigned_healthy_niche), .(assigned_healthy_niche, assigned_niche_name)])
setorder(nmap, assigned_healthy_niche)
NICHE <- nmap$assigned_niche_name
cat("[niche] observed categories:", length(NICHE), " categories:",
    paste(sprintf("%s=%s", nmap$assigned_healthy_niche, NICHE), collapse = "; "), "\n")

NC_ALL <- c("Granulocytic-erythroid"      = "#D55E00",
            "Erythroid-rich"              = "#CC79A7",
            "Megakaryocyte-rich"          = "#E69F00",
            "Stromal/vascular-osteogenic" = "#0072B2",
            "Low-density mixed"           = "#009E73",

            "HSPC-enriched niche"              = "#0072B2",
            "Erythroid-myeloid niche"          = "#D55E00",
            "Erythroid-megakaryocytic niche"   = "#CC79A7")
unknown <- setdiff(NICHE, names(NC_ALL))
if (length(unknown)) stop("Niche missing from palette: ", paste(unknown, collapse = ", "))
NC <- NC_ALL[NICHE]
NICHE_ORD <- names(NC_ALL)[names(NC_ALL) %in% NICHE]
UNMAP_COL <- "grey82"

nn_px <- function(d, k = 400) {
  i <- seq_len(min(k, nrow(d)))
  median(vapply(i, function(j) {
    dd <- sqrt((d$x - d$x[j])^2 + (d$y - d$y[j])^2); min(dd[dd > 0])
  }, 0))
}
scl <- map[, .(nn = nn_px(.SD)), by = sid, .SDcols = c("x", "y")]
scl[, um_per_px := SPOT_UM / nn]
map <- merge(map, scl[, .(sid, nn, um_per_px)], by = "sid")
map[, `:=`(ux = x * um_per_px, uy = -y * um_per_px)]
map[, `:=`(ux = ux - min(ux), uy = uy - min(uy)), by = sid]
cat("[scale] "); cat(sprintf("%s %.1f um/px  ", scl$sid, scl$um_per_px), "\n")

PAD_UM <- SPOT_UM * 2.2
ext <- map[, .(xr = diff(range(ux)), yr = diff(range(uy)),
               xc = mean(range(ux)), yc = mean(range(uy))), by = sid]
Y_COM <- max(ext$yr) + 2 * PAD_UM + SPOT_UM * 3.2
ext[, `:=`(x0 = xc - (xr / 2 + PAD_UM), x1 = xc + (xr / 2 + PAD_UM),
           y0 = yc - Y_COM / 2,         y1 = yc + Y_COM / 2)]
ext[, wgt := (xr + 2 * PAD_UM)]
EXT <- split(ext, ext$sid)
W_ROW <- setNames(ext$wgt, as.character(ext$sid))[SAMP]

ROW_MM  <- 42
PT_SIZE <- max(0.22, ROW_MM * (SPOT_UM / Y_COM) / 0.75 * 0.95)
PS <- setNames(rep(PT_SIZE, length(SAMP)), SAMP)
cat(sprintf("[scale] shared scale: Y_COM=%.0f um, point size=%.2f, width weights %s\n",
            Y_COM, PT_SIZE, paste(round(W_ROW), collapse = "/")))

base_thm <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.line = element_line(linewidth = .3, colour = "black"),
        axis.ticks = element_line(linewidth = .3, colour = "black"),
        axis.ticks.length = unit(.9, "mm"),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        plot.title = element_text(size = 7, face = "plain", hjust = .5,
                                  margin = margin(b = 1.2)),
        plot.margin = margin(1, 2, 1, 2),
        legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0))
map_thm <- theme_void(base_size = 7, base_family = FONT) +
  theme(plot.margin = margin(1, 1, 1, 1),
        legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(3, "mm"))

fix_coord <- function(s) {
  e <- EXT[[s]]
  coord_fixed(xlim = c(e$x0, e$x1), ylim = c(e$y0, e$y1), expand = FALSE, clip = "off")
}
scalebar <- function(s, show) {
  if (!show) return(NULL)
  e <- EXT[[s]]; y <- e$y0 + SPOT_UM * 1.6; x <- e$x0 + SPOT_UM * 1.2
  list(annotate("segment", x = x, xend = x + BAR_UM, y = y, yend = y,
                linewidth = .6, colour = "black", lineend = "butt"),
       annotate("text", x = x + BAR_UM / 2, y = y + SPOT_UM * 1.0,
                label = "1 mm", size = 1.75, family = FONT, vjust = 0))
}

inlab <- function(s, n) {
  e <- EXT[[s]]
  annotate("text", x = e$x0 + SPOT_UM * 1.2, y = e$y1 - SPOT_UM * 1.2,
           label = sprintf("%s\nn = %s", DL(s), format(n, big.mark = ",")),
           hjust = 0, vjust = 1, size = 2.05, family = FONT, lineheight = .95)
}

UNMAP_LAB <- sprintf("Unmapped (n = %d)", sum(is.na(map$assigned_healthy_niche)))
A_COL <- c(NC, setNames(UNMAP_COL, UNMAP_LAB))
panel_a <- function(s, show_leg = FALSE, show_bar = FALSE) {
  d <- copy(map[sid == s])
  d[, lab := fifelse(is.na(assigned_niche_name), UNMAP_LAB, assigned_niche_name)]
  d[, lab := factor(lab, levels = c(NICHE_ORD, UNMAP_LAB))]
  setorder(d, lab)
  ggplot(d, aes(ux, uy, colour = lab)) +
    geom_point(size = PS[[s]], shape = 16) +
    scale_colour_manual(values = A_COL, breaks = c(NICHE_ORD, UNMAP_LAB),
                        name = "Assigned healthy niche", drop = FALSE,
                        guide = guide_legend(override.aes = list(size = 1.7), ncol = 1,
                                             keyheight = unit(3.1, "mm"))) +
    fix_coord(s) + scalebar(s, show_bar) + inlab(s, nrow(d)) + map_thm +
    theme(legend.position = if (show_leg) "right" else "none")
}

SR_COL <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B")
panel_b <- function(s, show_leg = FALSE, show_bar = FALSE) {
  d  <- map[sid == s]
  dm <- d[is.finite(sRNHI_z_global)]
  du <- d[!is.finite(sRNHI_z_global)]
  ggplot() +
    { if (nrow(du)) geom_point(data = du, aes(ux, uy), colour = UNMAP_COL,
                               size = PS[[s]], shape = 16) } +
    geom_point(data = dm, aes(ux, uy, colour = sRNHI_z_global),
               size = PS[[s]], shape = 16) +
    scale_colour_gradientn(
      colours = SR_COL, limits = c(-Z_LIM, Z_LIM), oob = scales::squish,
      breaks = c(-Z_LIM, 0, Z_LIM),
      labels = c(paste0("\u2264 ", -Z_LIM), "0", paste0("\u2265 ", Z_LIM)),
      name = "sRNHI (z, global)",
      guide = guide_colourbar(raster = FALSE, barheight = unit(20, "mm"), barwidth = unit(2.4, "mm"),
                              ticks.colour = "white", frame.colour = "grey30",
                              frame.linewidth = .3)) +
    fix_coord(s) + scalebar(s, show_bar) + inlab(s, nrow(dm)) + map_thm +
    theme(legend.position = if (show_leg) "right" else "none")
}

gmu <- mean(hld$sRNHI_raw); gsd <- sd(hld$sRNHI_raw)
hz <- data.table(donor = hld$donor, zp = (hld$sRNHI_raw - gmu) / gsd)
H_MED <- median(hz$zp); H_Q <- quantile(hz$zp, c(.25, .75))
cat(sprintf("[healthy ref] n=%d  median=%.2f  IQR=[%.2f, %.2f]\n",
            nrow(hz), H_MED, H_Q[1], H_Q[2]))

dz <- map[is.finite(sRNHI_z_global)]
stat_c <- dz[, {
  w <- suppressWarnings(wilcox.test(sRNHI_z_global, hz$zp))
  .(n = .N, med = median(sRNHI_z_global), p = w$p.value)
}, by = sid]
stat_c[, p_bh := p.adjust(p, "BH")]
stat_c[, star := fifelse(p_bh < 1e-4, "****", fifelse(p_bh < 1e-3, "***",
               fifelse(p_bh < 0.01, "**", fifelse(p_bh < ALPHA, "*", "ns"))))]
print(stat_c[, .(sid, n, med = round(med, 2), p_bh = signif(p_bh, 2), star)])
dz <- merge(dz, stat_c[, .(sid, med_fill = med)], by = "sid")

stat_c[, delta := med - H_MED]
ytop <- max(dz$sRNHI_z_global)
panel_c <- ggplot(dz, aes(sid, sRNHI_z_global)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = H_Q[1], ymax = H_Q[2],
           fill = "grey86", alpha = .8) +
  annotate("segment", x = -Inf, xend = Inf, y = H_MED, yend = H_MED,
           linetype = "22", linewidth = .35, colour = "grey25") +
  geom_violin(aes(fill = med_fill), scale = "width", width = .86,
              linewidth = .22, colour = "grey45") +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white",
               linewidth = .32, fatten = 1.5) +
  geom_text(data = stat_c, aes(sid, ytop + .22, label = star),
            size = 2.0, family = FONT, inherit.aes = FALSE) +
  geom_text(data = stat_c, aes(sid, ytop + .62, label = sprintf("%+.2f", delta)),
            size = 1.85, family = FONT, colour = "grey20", inherit.aes = FALSE) +
  annotate("text", x = 4.58, y = H_MED, label = "healthy\nreference",
           size = 1.75, family = FONT, colour = "grey30", hjust = 0, lineheight = .95) +
  scale_fill_gradientn(colours = SR_COL, limits = c(-Z_LIM, Z_LIM),
                       oob = scales::squish, guide = "none",
                       aesthetics = "fill") +
  scale_x_discrete(labels = SAMP_LAB) +
  scale_y_continuous(expand = expansion(mult = c(.04, .10))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "sRNHI (z, global)") + base_thm +
  theme(plot.margin = margin(1, 26, 1, 2))

dq <- map[is.finite(sRNHI_z_global) & is.finite(Local_AML_fraction)]
stat_d <- dq[, {
  ct <- suppressWarnings(cor.test(sRNHI_z_global, Local_AML_fraction, method = "spearman"))
  .(rho = unname(ct$estimate), p = ct$p.value, n = .N)
}, by = sid]
stat_d[, p_bh := p.adjust(p, "BH")]

set.seed(20260101)
stat_d[, `:=`(ci_lo = NA_real_, ci_hi = NA_real_)]
for (k in seq_len(nrow(stat_d))) {
  d  <- dq[sid == stat_d$sid[k]]
  nb <- max(20L, round(nrow(d) / 50))
  cl <- try(stats::kmeans(cbind(d$ux, d$uy), centers = nb, iter.max = 50)$cluster, silent = TRUE)
  if (inherits(cl, "try-error")) cl <- rep(1L, nrow(d))
  bl <- split(seq_len(nrow(d)), cl)
  bs <- vapply(seq_len(400), function(i) {
    idx <- unlist(bl[sample.int(length(bl), length(bl), replace = TRUE)], use.names = FALSE)
    suppressWarnings(cor(d$sRNHI_z_global[idx], d$Local_AML_fraction[idx], method = "spearman"))
  }, 0)
  stat_d$ci_lo[k] <- quantile(bs, .025, na.rm = TRUE)
  stat_d$ci_hi[k] <- quantile(bs, .975, na.rm = TRUE)
}
stat_d[, lab := sprintf("\u03c1 = %+.2f [%.2f, %.2f]", rho, ci_lo, ci_hi)]
print(stat_d[, .(sid, n, rho = round(rho, 3), ci = sprintf("[%+.2f, %+.2f]", ci_lo, ci_hi), p_bh = signif(p_bh, 2))])
print(stat_d[, .(sid, n, rho = round(rho, 3), p_bh = signif(p_bh, 2))])

panel_d <- ggplot(dq, aes(Local_AML_fraction, sRNHI_z_global)) +
  annotate("segment", x = -Inf, xend = Inf, y = H_MED, yend = H_MED,
           linetype = "22", linewidth = .3, colour = "grey55") +
  geom_point(size = .18, alpha = .22, colour = "grey25", stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey15",
              fill = "grey55", alpha = .30, linewidth = .4) +
  geom_text(data = stat_d, aes(0.02, 2.30, label = lab), parse = FALSE,
            hjust = 0, size = 1.62, family = FONT, inherit.aes = FALSE) +
  facet_wrap(~sid, nrow = 1, labeller = as_labeller(SAMP_LAB)) +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(-3.1, 2.6)) +
  labs(x = "Local AML fraction", y = "sRNHI (z, global)") +
  base_thm + theme(strip.background = element_blank(),
                   strip.text = element_text(size = 6.5, margin = margin(b = 1)),
                   panel.spacing.x = unit(2.4, "mm"))

val[, rho := as.numeric(spearman_rho)]
val[, p   := as.numeric(spearman_p)]
val[, p_bh := p.adjust(p, "BH")]
val[, star := fifelse(p_bh < 1e-4, "****", fifelse(p_bh < 1e-3, "***",
              fifelse(p_bh < 0.01, "**", fifelse(p_bh < ALPHA, "*", ""))))]
val[, sid := factor(sample, levels = SAMP)]
ord <- val[, .(m = mean(rho, na.rm = TRUE)), by = program][order(-m), program]
val[, program := factor(program, levels = ord)]
PROG_LAB <- c(HSPC = "HSPC", Erythroid = "Erythroid", Granulo = "Granulocytic",
              MonoMacro = "Mono/macro", LymphB = "B cell", LymphTNK = "T/NK",
              EndoVasc = "Endothelial", StromalMSC = "Stromal/MSC", Osteo = "Osteogenic",
              Adipo = "Adipocyte", Megakaryo = "Megakaryocytic", Inflammation = "Inflammation")
RL <- max(abs(val$rho), na.rm = TRUE)

panel_e <- ggplot(val, aes(program, sid, fill = rho)) +
  geom_tile(colour = "white", linewidth = .45) +
  geom_text(aes(label = star), size = 1.9, family = FONT, vjust = .78, colour = "grey10") +
  scale_fill_gradientn(colours = SR_COL, limits = c(-RL, RL), oob = scales::squish,
                       breaks = c(-.4, -.2, 0, .2, .4),
                       name = expression(bold(Spearman~rho)),
                       guide = guide_colourbar(raster = FALSE, barheight = unit(2.4, "mm"),
                                               barwidth = unit(24, "mm"),
                                               title.position = "left", title.vjust = 1,
                                               ticks.colour = "white",
                                               frame.colour = "grey30", frame.linewidth = .3)) +
  scale_x_discrete(labels = PROG_LAB, expand = c(0, 0)) +
  scale_y_discrete(limits = rev(SAMP), labels = SAMP_LAB, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  base_thm +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 6),
        legend.position = "bottom", legend.title = element_text(size = 6.5, face = "bold"),
        legend.margin = margin(t = -3), plot.margin = margin(1, 2, 0, 2))

pad <- theme(plot.margin = margin(t = 4.5, r = 1, b = 0, l = 1))
lab_thm <- theme(plot.tag = element_text(size = 8, face = "bold", family = FONT,
                                         hjust = 0, vjust = 1))

LEG_W <- sum(W_ROW) * 0.34
rowA <- wrap_plots(list(panel_a(SAMP[1], FALSE, TRUE), panel_a(SAMP[2]),
                        panel_a(SAMP[3]), panel_a(SAMP[4], TRUE)),
                   nrow = 1, widths = as.numeric(W_ROW)) +
        plot_layout(guides = "collect") & theme(legend.position = "right")
rowB <- wrap_plots(list(panel_b(SAMP[1], FALSE, TRUE), panel_b(SAMP[2]),
                        panel_b(SAMP[3]), panel_b(SAMP[4], TRUE)),
                   nrow = 1, widths = as.numeric(W_ROW)) +
        plot_layout(guides = "collect") & theme(legend.position = "right")

fig <- quiet_mbcs((wrap_elements(full = patchworkGrob(rowA)) + pad) /
       (wrap_elements(full = patchworkGrob(rowB)) + pad)) /
       (panel_c + panel_d + plot_layout(widths = c(1, 2.35))) /
       panel_e +
  plot_layout(heights = c(1.00, 1.00, 1.02, 0.86)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))

H_FIG <- MM(205)
quiet_mbcs({
  open_pdf(file.path(OUT, "FigureC_main.pdf"), W_FULL, H_FIG); print(fig); dev.off()
  open_png(file.path(OUT, "FigureC_main.png"), W_FULL, H_FIG); print(fig); dev.off()
})
cat("[fig] FigureC_main  183 x", round(H_FIG * 25.4), "mm\n")

FR_COL <- viridisLite::mako(9, direction = -1)

panel_frac <- function(s, show_leg = FALSE, show_bar = FALSE) {
  d  <- map[sid == s]; dm <- d[is.finite(Local_AML_fraction)]
  du <- d[!is.finite(Local_AML_fraction)]
  ggplot() +
    { if (nrow(du)) geom_point(data = du, aes(ux, uy), colour = UNMAP_COL,
                               size = PS[[s]], shape = 16) } +
    geom_point(data = dm, aes(ux, uy, colour = Local_AML_fraction),
               size = PS[[s]], shape = 16) +
    scale_colour_gradientn(colours = FR_COL, limits = c(0, 1), oob = scales::squish,
      breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), name = "Local AML fraction",
      guide = guide_colourbar(raster = FALSE, barheight = unit(20, "mm"), barwidth = unit(2.4, "mm"),
                              ticks.colour = "white", frame.colour = "grey30",
                              frame.linewidth = .3)) +
    fix_coord(s) + scalebar(s, show_bar) + inlab(s, nrow(dm)) + map_thm +
    theme(legend.position = if (show_leg) "right" else "none")
}

panel_zg <- function(s, show_leg = FALSE, show_bar = FALSE) {
  d  <- map[sid == s]; dm <- d[is.finite(sRNHI_z_perniche)]
  du <- d[!is.finite(sRNHI_z_perniche)]
  ggplot() +
    { if (nrow(du)) geom_point(data = du, aes(ux, uy), colour = UNMAP_COL,
                               size = PS[[s]], shape = 16) } +
    geom_point(data = dm, aes(ux, uy, colour = sRNHI_z_perniche),
               size = PS[[s]], shape = 16) +
    scale_colour_gradientn(colours = SR_COL, limits = c(-Z_LIM, Z_LIM),
      oob = scales::squish, breaks = c(-Z_LIM, 0, Z_LIM),
      labels = c(paste0("\u2264 ", -Z_LIM), "0", paste0("\u2265 ", Z_LIM)),
      name = "sRNHI (z, per niche)",
      guide = guide_colourbar(raster = FALSE, barheight = unit(20, "mm"), barwidth = unit(2.4, "mm"),
                              ticks.colour = "white", frame.colour = "grey30",
                              frame.linewidth = .3)) +
    fix_coord(s) + scalebar(s, show_bar) + inlab(s, nrow(dm)) + map_thm +
    theme(legend.position = if (show_leg) "right" else "none")
}

mkrow <- function(fn) wrap_plots(list(fn(SAMP[1], FALSE, TRUE), fn(SAMP[2]),
                                      fn(SAMP[3]), fn(SAMP[4], TRUE)),
                                 nrow = 1, widths = as.numeric(W_ROW)) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

figS2 <- quiet_mbcs((wrap_elements(full = patchworkGrob(mkrow(panel_frac))) + pad) /
         (wrap_elements(full = patchworkGrob(mkrow(panel_zg))) + pad)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))

H_S2 <- MM(96)
quiet_mbcs({
  open_pdf(file.path(OUT, "FigureC_S2_maps.pdf"), W_FULL, H_S2); print(figS2); dev.off()
  open_png(file.path(OUT, "FigureC_S2_maps.png"), W_FULL, H_S2); print(figS2); dev.off()
})
cat("[fig] FigureC_S2_maps  183 x", round(H_S2 * 25.4), "mm\n")

hz[, dlab := anonymous_ids(donor)]
hz[, dlab := anonymous_ids(donor)]
s1a <- ggplot(hz, aes(dlab, zp)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = H_Q[1], ymax = H_Q[2],
           fill = "grey88", alpha = .8) +
  annotate("segment", x = -Inf, xend = Inf, y = H_MED, yend = H_MED,
           linetype = "22", linewidth = .32, colour = "grey25") +
  geom_boxplot(width = .6, outlier.size = .12, outlier.alpha = .2,
               fill = "white", linewidth = .3, fatten = 1.4) +
  labs(x = NULL, y = "sRNHI (z, global)",
       title = sprintf("Healthy hold-out calibration (%d donors, n = %s)",
                       uniqueN(hz$donor), format(nrow(hz), big.mark = ","))) +
  coord_cartesian(ylim = quantile(hz$zp, c(.001, .999))) +
  base_thm + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 5.4))

comp <- map[!is.na(assigned_niche_name), .N, by = .(sid, assigned_niche_name)]
comp[, frac := N / sum(N), by = sid]
comp[, assigned_niche_name := factor(assigned_niche_name, levels = rev(NICHE_ORD))]
s1b <- ggplot(comp, aes(sid, frac, fill = assigned_niche_name)) +
  geom_col(width = .72, colour = "white", linewidth = .3) +
  scale_fill_manual(values = NC, breaks = NICHE_ORD, name = NULL,
                    guide = guide_legend(ncol = 1, keyheight = unit(3, "mm"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, .02))) +
  scale_x_discrete(labels = SAMP_LAB) +
  labs(x = NULL, y = "Spots assigned", title = "Healthy-niche composition per sample") +
  base_thm + theme(legend.position = "right", legend.text = element_text(size = 5.8))

dn <- map[is.finite(sRNHI_z_global) & is.finite(Normal_abundance)]
dn[, lowconf := Normal_abundance < 0.05]
st_n <- dn[, {ct <- suppressWarnings(cor.test(Normal_abundance, sRNHI_z_global, method = "spearman"))
             .(lab = sprintf('rho == "%+.2f"', unname(ct$estimate)),
               nl = sprintf("low-confidence: %.1f%%", 100 * mean(lowconf)))}, by = sid]
s1c <- ggplot(dn, aes(Normal_abundance, sRNHI_z_global)) +
  geom_point(aes(colour = lowconf), size = .18, alpha = .25, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey15",
              fill = "grey55", alpha = .3, linewidth = .4) +
  geom_text(data = st_n, aes(0.02, 2.3, label = lab), parse = TRUE, hjust = 0,
            size = 1.85, family = FONT, inherit.aes = FALSE) +
  geom_text(data = st_n, aes(0.02, 1.75, label = nl), hjust = 0, size = 1.7,
            family = FONT, colour = "#B2182B", inherit.aes = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B2182B"),
                      labels = c("Normal abundance \u2265 0.05", "< 0.05 (low confidence)"),
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.4, alpha = 1))) +
  facet_wrap(~sid, nrow = 1, labeller = as_labeller(SAMP_LAB)) +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  coord_cartesian(xlim = c(0, 1.08), ylim = c(-3.1, 2.6)) +
  labs(x = "Normal-cell abundance (RCTD, absolute)", y = "sRNHI (z, global)") +
  base_thm + theme(strip.background = element_blank(),
                   strip.text = element_text(size = 6.5, margin = margin(b = 1)),
                   panel.spacing.x = unit(2.4, "mm"), legend.position = "bottom")

s1d <- ggplot(map[is.finite(assignment_confidence)], aes(sid, assignment_confidence)) +
  geom_violin(scale = "width", width = .8, fill = "grey90",
              colour = "grey45", linewidth = .22) +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white",
               linewidth = .3, fatten = 1.4) +
  scale_x_discrete(labels = SAMP_LAB) +
  labs(x = NULL, y = "Niche-assignment confidence",
       title = "Assignment confidence") + base_thm

s1e <- ggplot(map[is.finite(Local_AML_fraction)], aes(sid, Local_AML_fraction)) +
  geom_violin(scale = "width", width = .8, fill = "grey90",
              colour = "grey45", linewidth = .22) +
  geom_boxplot(width = .16, outlier.shape = NA, fill = "white",
               linewidth = .3, fatten = 1.4) +
  scale_x_discrete(labels = SAMP_LAB) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Local AML fraction",
       title = "Leukaemic burden per spot") + base_thm

figS <- (s1a | s1b) / s1c / (s1d | s1e) +
  plot_layout(heights = c(1, 1.05, 0.85), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))

H_S <- MM(170)
quiet_mbcs({
  open_pdf(file.path(OUT, "FigureC_S1_QC.pdf"), W_FULL, H_S); print(figS); dev.off()
  open_png(file.path(OUT, "FigureC_S1_QC.png"), W_FULL, H_S); print(figS); dev.off()
})
cat("[fig] FigureC_S1_QC  183 x", round(H_S * 25.4), "mm\n")

if (EXPORT_PANELS) {
  pdir <- file.path(OUT, "panels"); dir.create(pdir, showWarnings = FALSE)
  one <- function(nm, obj, w, h) {
    quiet_mbcs({ open_pdf(file.path(pdir, paste0(nm, ".pdf")), w, h); print(obj); dev.off() })
  }
  one("C_a_niche_maps",   quiet_mbcs(wrap_elements(full = patchworkGrob(rowA))), W_FULL, MM(46))
  one("C_b_sRNHI_maps",   quiet_mbcs(wrap_elements(full = patchworkGrob(rowB))), W_FULL, MM(46))
  one("C_c_distribution", panel_c,                                           MM(62), MM(52))
  one("C_d_orthogonality", panel_d,                                          MM(120), MM(52))
  one("C_e_modules",      panel_e,                                           W_FULL, MM(46))
  one("S2_a_AMLfraction", quiet_mbcs(wrap_elements(full = patchworkGrob(mkrow(panel_frac)))), W_FULL, MM(46))
  one("S2_b_sRNHI_perniche", quiet_mbcs(wrap_elements(full = patchworkGrob(mkrow(panel_zg)))),  W_FULL, MM(46))
  cat("Output completed.\n")
}

calib <- map[is.finite(sRNHI_z_global) & is.finite(sRNHI_z_perniche), {
  ct <- suppressWarnings(cor.test(sRNHI_z_global, sRNHI_z_perniche, method = "spearman"))
  hpg <- (hld$sRNHI_raw - gmu) / gsd
  hp  <- (hld$sRNHI_raw - res$mean[match(hld$assigned_domain, res$domain)]) / res$sd[match(hld$assigned_domain, res$domain)]
  .(rho = unname(ct$estimate), n = .N,
    pres_global = mean(sRNHI_z_global >= quantile(hpg, .10)),
    pres_perniche = mean(sRNHI_z_perniche >= quantile(hp, .10)))
}, by = sid]
print(calib[, .(sid, n, rho = round(rho, 3),
                pres_global = round(pres_global, 3), pres_perniche = round(pres_perniche, 3))])

out_stats <- rbindlist(list(
  stat_c[, .(panel = "c", key = unname(SAMP_LAB[as.character(sid)]), n, value = round(med, 3),
             metric = "median sRNHI(z, global)", p_bh = signif(p_bh, 3))],
  stat_d[, .(panel = "d", key = unname(SAMP_LAB[as.character(sid)]), n, value = round(rho, 3),
             metric = "Spearman rho vs Local AML fraction", p_bh = signif(p_bh, 3))],
  calib[, .(panel = "calib", key = unname(SAMP_LAB[as.character(sid)]), n, value = round(rho, 3),
            metric = "Spearman rho: global vs per-niche calibration", p_bh = NA_real_)],
  val[, .(panel = "e", key = paste(unname(SAMP_LAB[sample]), program, sep = "|"),
          n = NA_integer_, value = round(rho, 3),
          metric = "Spearman rho vs sRNHI", p_bh = signif(p_bh, 3))]))
fwrite(out_stats, file.path(OUT, "FigureC_source_statistics.csv"))
cat("[done]\n")

if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))
