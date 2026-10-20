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

suppressPackageStartupMessages({library(ggplot2); library(data.table); library(patchwork)})

res <- "output/figure_4/composition"
dir.create(res, recursive = TRUE, showWarnings = FALSE)
d   <- fread(file.path(res, "Figure6B_patient_lineage_prop.csv"))
d$unit <- anonymous_ids(d$unit)
tst <- fread(file.path(res, "Figure6B_lineage_tests.csv"))
uni <- fread(file.path(res, "Figure6B_unit_summary.csv"))
uni$unit <- anonymous_ids(uni$unit)

states <- c("Normal_HSC_MPP","Normal_LMPP","Normal_MyeloidProgenitor","Normal_Monocyte",
            "Normal_Granulocytic","Normal_Megakaryocytic","Normal_EarlyErythroid","Normal_LateErythroid",
            "Normal_B","Normal_T","Normal_NK","Normal_Dendritic")
cols12 <- c("Normal_HSC_MPP"="#E64B35","Normal_LMPP"="#F39B7F","Normal_MyeloidProgenitor"="#4DBBD5",
            "Normal_Monocyte"="#00A087","Normal_Granulocytic"="#3C5488","Normal_Megakaryocytic"="#7E6148",
            "Normal_EarlyErythroid"="#B2182B","Normal_LateErythroid"="#D6604D","Normal_B"="#4393C3",
            "Normal_T"="#92C5DE","Normal_NK"="#F4A582","Normal_Dendritic"="#FDDBC7")
lab12 <- c("Normal_HSC_MPP"="HSC/MPP","Normal_LMPP"="LMPP","Normal_MyeloidProgenitor"="Myeloid prog.",
           "Normal_Monocyte"="Monocyte","Normal_Granulocytic"="Granulocytic","Normal_Megakaryocytic"="Megakaryocytic",
           "Normal_EarlyErythroid"="Early erythroid","Normal_LateErythroid"="Late erythroid","Normal_B"="B",
           "Normal_T"="T","Normal_NK"="NK","Normal_Dendritic"="Dendritic")
grpCols <- c("Healthy"="#009E73","RNHI-high AML"="#D55E00","RNHI-low AML"="#CC79A7")
fmt2 <- function(x) formatC(x, format = "e", digits = 1)

d[, state12 := factor(state12, levels = states)]
d[, grpDisp := factor(grpDisp, levels = names(grpCols))]
uni[, grpDisp := factor(grpDisp, levels = names(grpCols))]
setorder(uni, unit_order)
uni[, pos := .I]
d[, unit := factor(unit, levels = uni$unit)]

seg <- d[, .(pct = 100 * prop), by = .(unit, grpDisp, state12, unit_order)]

gp <- uni[, .(x0 = min(pos) - 0.5, x1 = max(pos) + 0.5,
              n = .N, ncells = sum(n_cells)), by = grpDisp]
gp[, mid := (x0 + x1) / 2]
gp[, txt := ifelse(grpDisp == "Healthy", "Healthy donors (n=54)",
            ifelse(grpDisp == "RNHI-high AML", "Hi", "Lo"))]
gpCol <- grpCols[as.character(gp$grpDisp)]

pA <- ggplot(seg, aes(x = unit, y = pct, fill = state12)) +
  geom_col(width = 1, colour = "white", linewidth = 0.08) +
  geom_rect(data = gp, aes(xmin = x0, xmax = x1, ymin = -10.5, ymax = 0),
            inherit.aes = FALSE, fill = gpCol, alpha = 0.95) +
  geom_text(data = gp, aes(x = mid, y = -5.25, label = txt), inherit.aes = FALSE,
            colour = "white", size = 3.1, fontface = "bold") +
  scale_fill_manual(values = cols12, name = NULL, breaks = states, labels = lab12[states]) +
  scale_y_continuous(limits = c(-10.5, NA), breaks = seq(0, 100, 25),
                     labels = paste0(seq(0, 100, 25), "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(drop = FALSE) +
  guides(fill = guide_legend(ncol = 4, byrow = TRUE, override.aes = list(colour = NA))) +
  coord_cartesian(ylim = c(-10.5, 102), clip = "off") +
  labs(x = NULL, y = "Composition (% of residual normal cells)",
       title = "Patient-level residual normal composition",
       subtitle = "64 units: Healthy 54 donors -> RNHI-high AML (n=5) -> RNHI-low AML (n=5);\nAML ordered by RNHI high->low. Cells: Healthy 338,738 / Hi 5,835 / Lo 10,892.",
       caption = "Strip colours: green Healthy \u00b7 orange RNHI-high AML (Hi) \u00b7 purple RNHI-low AML (Lo)\nBars = 12 normal states (legend below); each column = one donor/patient, 100% stacked.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line.y = element_line(colour = "grey45", linewidth = 0.4),
        axis.text.y = element_text(size = 8.6, colour = "black"),
        axis.title = element_text(size = 9.6),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8.4, hjust = 0.5, colour = "grey25", lineheight = 1.05),
        plot.caption = element_text(size = 6.8, hjust = 0.5, colour = "grey35", lineheight = 1.15),
        legend.position = "bottom", legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.28, "cm"), legend.key.spacing.x = unit(0.08, "cm"),
        legend.margin = margin(1, 0, 0, 0), plot.margin = margin(6, 10, 3, 6))

sig <- list()
for (i in seq_len(nrow(tst))) {
  st <- tst$state12[i]
  lines <- character(0)
  if (isTRUE(tst$kw_sig[i])) {
    if (!is.na(tst$p_HvHi_bh[i]) && tst$p_HvHi_bh[i] < 0.05) lines <- c(lines, sprintf("Healthy\u2212Hi %s", fmt2(tst$p_HvHi_bh[i])))
    if (!is.na(tst$p_HvLo_bh[i]) && tst$p_HvLo_bh[i] < 0.05) lines <- c(lines, sprintf("Healthy\u2212Lo %s", fmt2(tst$p_HvLo_bh[i])))
    if (!is.na(tst$p_HiLo_bh[i]) && tst$p_HiLo_bh[i] < 0.05) lines <- c(lines, sprintf("Hi\u2212Lo %s", fmt2(tst$p_HiLo_bh[i])))
  }
  sig[[length(sig) + 1L]] <- data.table(state12 = st,
    kw_txt = ifelse(isTRUE(tst$kw_sig[i]), sprintf("KW p = %s", fmt2(tst$kw_p[i])), sprintf("KW p = %s (ns)", fmt2(tst$kw_p[i]))),
    pair_txt = if (length(lines)) paste(lines, collapse = "\n") else NA_character_,
    has_pair = length(lines) > 0L)
}
sa <- rbindlist(sig)

d2 <- d[, .(state12, grpDisp, pct = 100 * prop)]
mx <- d2[, .(mx = max(pct)), by = state12]
sa <- merge(sa, mx, by = "state12")
sa[, kw_y := mx * 1.06 + 0.5]
sa[, pair_y := kw_y + pmax(mx * 0.06, 3.5)]
sa[has_pair == FALSE, pair_y := NA_real_]

xshort <- c("Healthy"="Healthy","RNHI-high AML"="RNHI-high","RNHI-low AML"="RNHI-low")
pB <- ggplot(d2, aes(x = grpDisp, y = pct, fill = grpDisp, colour = grpDisp)) +
  geom_violin(scale = "width", width = 0.9, alpha = 0.30, linewidth = 0.35) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
  geom_jitter(width = 0.10, size = 0.55, alpha = 0.9, stroke = 0) +
  facet_wrap(~ state12, ncol = 4, scales = "free_y", labeller = labeller(state12 = lab12)) +
  geom_text(data = sa, aes(x = 2, y = kw_y, label = kw_txt), inherit.aes = FALSE,
            size = 2.5, hjust = 0.5, colour = "grey25") +
  geom_text(data = sa[has_pair == TRUE], aes(x = 2, y = pair_y, label = pair_txt), inherit.aes = FALSE,
            size = 2.6, hjust = 0.5, colour = "#A50F15", fontface = "bold", lineheight = 0.92) +
  scale_fill_manual(values = grpCols, name = NULL,
                    labels = c("Healthy"="Healthy (n=54)","RNHI-high AML"="RNHI-high AML (n=5)","RNHI-low AML"="RNHI-low AML (n=5)")) +
  scale_colour_manual(values = grpCols, guide = "none") +
  scale_x_discrete(labels = xshort) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Lineage proportion per patient (%)",
       title = "Residual normal lineage proportions by group (per patient)",
       subtitle = "Violin + box + per-patient dots. Panel text: Kruskal\u2013Wallis p; red lines: pairwise Wilcoxon p after global BH correction (<0.05).") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8.2, face = "bold"),
        axis.text.x = element_text(size = 7.6, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 9.4),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8.2, hjust = 0.5, colour = "grey25", lineheight = 1.05),
        legend.position = "bottom", legend.text = element_text(size = 8.4),
        legend.key.size = unit(0.42, "cm"), legend.margin = margin(2, 0, 0, 0),
        panel.spacing = unit(0.9, "lines"),
        plot.margin = margin(5, 10, 3, 6))
pB <- pB + guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 1.2, shape = 21)))

emit <- function(p, base, w, h) {
  pdf_f <- file.path(res, paste0(base, ".pdf")); png_f <- file.path(res, paste0(base, ".png"))
  grDevices::cairo_pdf(pdf_f, width = w, height = h, family = "Helvetica"); print(p); grDevices::dev.off()
  grDevices::png(png_f, width = w, height = h, units = "in", res = 300, type = "cairo"); print(p); grDevices::dev.off()
  cat("Output completed.\n")
}
emit(pA, "Figure6B_panelA_composition_stacked", 4.5, 4.5)
emit(pB, "Figure6B_panelB_lineage_violin", 9.6, 7.4)
comp <- pA + pB + plot_layout(ncol = 1, heights = c(0.78, 1.6)) +
  plot_annotation(title = "Figure 6B  Residual normal haematopoietic composition: Healthy vs RNHI-high/low AML",
                  theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))
emit(comp, "Figure6B_composite", 9.6, 11.8)
invisible(list(A = pA, B = pB))
