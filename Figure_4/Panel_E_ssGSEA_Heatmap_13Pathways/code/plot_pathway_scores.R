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

ROW_GROUPING <- "group"

FRAME_MODE <- "block"

DELTA      <- 0.05
SHOW_TITLE <- FALSE
NCOL_B     <- as.integer(Sys.getenv("F6D_NCOL_B", "4"))

LABEL_WRAP <- Sys.getenv("F6D_LABEL_WRAP", "word")

LABEL_LH   <- 1.12
LABEL_PT   <- 6.5
ROW_MIN_MM <- 5.0

FONT       <- Sys.getenv("F6D_FONT", "")
if (!nzchar(FONT)) FONT <- local({
  cand <- c("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "sans")
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
MM         <- function(mm) mm / 25.4
W_FULL     <- MM(183)

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(ComplexHeatmap)
  library(circlize); library(grid)
})

res <- Sys.getenv("F6D_RES", "output/figure_4/ssgsea")
res_in <- Sys.getenv("F6D_RES_IN", "output/figure_4/ssgsea")
stopifnot(dir.exists(res), dir.exists(res_in))

open_png <- function(fp, w, h, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE))
    ragg::agg_png(fp, width = w, height = h, units = "in", res = dpi)
  else grDevices::png(fp, width = w, height = h, units = "in", res = dpi,
                      type = if (capabilities("cairo")) "cairo" else "Xlib")
}

quiet_mbcs <- function(expr) withCallingHandlers(expr, warning = function(w)
  if (grepl("mbcsToSbcs|conversion failure", conditionMessage(w))) invokeRestart("muffleWarning"))

PDF_MODE <- Sys.getenv("F6D_PDF_MODE", "text")
if (PDF_MODE == "outline") {
  stopifnot(requireNamespace("showtext", quietly = TRUE),
            requireNamespace("sysfonts", quietly = TRUE))
  .ff <- function(bold = FALSE, italic = FALSE) {
    p <- tryCatch(systemfonts::match_font(FONT, bold = bold, italic = italic)$path,
                  error = function(e) NULL)
    if (is.null(p) || !file.exists(p)) systemfonts::match_font(FONT)$path else p
  }
  sysfonts::font_add("F6DOut", regular = .ff(), bold = .ff(bold = TRUE),
                     italic = .ff(italic = TRUE), bolditalic = .ff(bold = TRUE, italic = TRUE))
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 600)
  cat("Output completed.\n")
}
open_pdf <- function(fp, w, h) {

  if (PDF_MODE == "svg") {
    stopifnot(requireNamespace("svglite", quietly = TRUE))
    return(svglite::svglite(sub("\\.pdf$", ".svg", fp), width = w, height = h, bg = "white",
                             system_fonts = list(sans = FONT, serif = FONT,
                                                 mono = "Courier New", symbol = FONT)))
  }
  grDevices::cairo_pdf(fp, width = w, height = h,
                       family = if (PDF_MODE == "outline") "F6DOut" else FONT)
}

dev_null <- function(...) if (capabilities("cairo"))
  grDevices::cairo_pdf(filename = nullfile(), ...) else grDevices::pdf(nullfile(), ...)

pwL <- c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
         "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
         "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
         "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_E2F_TARGETS",
         "HALLMARK_G2M_CHECKPOINT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
         "HALLMARK_GLYCOLYSIS", "HALLMARK_MYC_TARGETS_V1")

pwLab <- c(
  HALLMARK_INFLAMMATORY_RESPONSE           = "Inflammatory response",
  HALLMARK_IL6_JAK_STAT3_SIGNALING         = "IL-6\u2013JAK\u2013STAT3 signalling",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB         = "TNF\u03B1 signalling via NF-\u03BAB",
  HALLMARK_INTERFERON_ALPHA_RESPONSE       = "Interferon-\u03B1 response",
  HALLMARK_INTERFERON_GAMMA_RESPONSE       = "Interferon-\u03B3 response",
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "Reactive oxygen species",
  HALLMARK_P53_PATHWAY                     = "p53 pathway",
  HALLMARK_APOPTOSIS                       = "Apoptosis",
  HALLMARK_E2F_TARGETS                     = "E2F targets",
  HALLMARK_G2M_CHECKPOINT                  = "G2/M checkpoint",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION       = "Oxidative phosphorylation",
  HALLMARK_GLYCOLYSIS                      = "Glycolysis",
  HALLMARK_MYC_TARGETS_V1                  = "MYC targets V1")

pwLabB <- c(
  HALLMARK_INFLAMMATORY_RESPONSE           = "Inflammatory",
  HALLMARK_IL6_JAK_STAT3_SIGNALING         = "IL-6\u2013JAK\u2013STAT3",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB         = "TNF\u03B1 via NF-\u03BAB",
  HALLMARK_INTERFERON_ALPHA_RESPONSE       = "Interferon-\u03B1",
  HALLMARK_INTERFERON_GAMMA_RESPONSE       = "Interferon-\u03B3",
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "ROS",
  HALLMARK_P53_PATHWAY                     = "p53 pathway",
  HALLMARK_APOPTOSIS                       = "Apoptosis",
  HALLMARK_E2F_TARGETS                     = "E2F targets",
  HALLMARK_G2M_CHECKPOINT                  = "G2/M checkpoint",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION       = "OXPHOS",
  HALLMARK_GLYCOLYSIS                      = "Glycolysis",
  HALLMARK_MYC_TARGETS_V1                  = "MYC targets V1")
stopifnot(all(pwL %in% names(pwLabB)))

split_words <- function(x) {
  w <- strsplit(x, "[[:space:]]+")[[1]]
  w[nzchar(w)]
}
wrap_label <- function(x) {
  w <- split_words(x)
  if (LABEL_WRAP == "none" || length(w) < 2L) return(paste(w, collapse = " "))
  if (LABEL_WRAP == "two") {
    k <- ceiling(length(w) / 2)
    return(paste(paste(w[seq_len(k)], collapse = " "),
                 paste(w[-seq_len(k)], collapse = " "), sep = "\n"))
  }
  paste(w, collapse = "\n")
}
n_lines <- function(x) {
  w <- split_words(x)
  if (LABEL_WRAP == "none" || length(w) < 2L) return(1L)
  if (LABEL_WRAP == "two") return(min(2L, length(w)))
  length(w)
}
LAB_FULL <- vapply(pwLab, wrap_label, character(1), USE.NAMES = TRUE)
cat(sprintf("[label] wrap=%s longest label %d lines total %d lines\n", LABEL_WRAP,
            max(vapply(pwLab, n_lines, 1L)),
            sum(vapply(pwLab, n_lines, 1L))))

fblocks <- setNames(c(rep("Inflammation", 3), rep("Interferon", 2),
                      rep("Stress/apoptosis", 3), rep("Proliferation", 3),
                      rep("Metabolism", 2)), pwL)
fcols <- c(Inflammation = "#E8D3C0", Interferon = "#CDE0F0",
           `Stress/apoptosis` = "#EDD3DA", Proliferation = "#D9E7D4",
           Metabolism = "#F1E7D5")

tint <- function(x, p = 0.20)
  grDevices::rgb(grDevices::colorRamp(c("white", x))(p), maxColorValue = 255)

gord <- c("Healthy", "RNHI-high AML", "RNHI-low AML")
gcol <- c("Healthy" = "#3B7A57", "RNHI-high AML" = "#C77C2F", "RNHI-low AML" = "#8E6FA3")
gshort <- c("Healthy" = "Healthy", "RNHI-high AML" = "RNHI-high", "RNHI-low AML" = "RNHI-low")

gshort2 <- c("Healthy" = "Healthy", "RNHI-high AML" = "RNHI-\nhigh", "RNHI-low AML" = "RNHI-\nlow")

read_sc <- function(cmp) {
  fp <- file.path(res_in, sprintf("Figure6D_ssgsea_%s.csv", cmp))
  if (!file.exists(fp)) stop("Missing input file: ", fp)
  x <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(x[[1]])) stop("Duplicate pathway names: ", fp)
  m <- as.matrix(x[, -1, drop = FALSE]); rownames(m) <- x[[1]]
  storage.mode(m) <- "double"
  m
}
sample_groups <- read.csv(file.path(res_in, "sample_groups.csv"), stringsAsFactors=FALSE, check.names=FALSE)
sample_groups$sample <- anonymous_ids(sample_groups$sample)
stopifnot(!anyDuplicated(sample_groups$sample), all(sample_groups$group %in% gord))
grp_of <- function(cols) {
  cols <- anonymous_ids(cols)
  idx <- match(cols, sample_groups$sample)
  stopifnot(!anyNA(idx))
  factor(sample_groups$group[idx], levels=gord)
}

row_z <- function(m) {
  mu <- rowMeans(m); sdv <- apply(m, 1, stats::sd)
  bad <- !is.finite(sdv) | sdv < .Machine$double.eps
  if (any(bad)) warning("Zero-variance pathways set to zero: ", paste(rownames(m)[bad], collapse = ", "))
  sdv[bad] <- 1
  z <- (m - mu) / sdv; z[bad, ] <- 0
  z
}

runs_of <- function(idx) {
  idx <- sort(unique(idx)); if (!length(idx)) return(list())
  b <- c(0L, which(diff(idx) != 1L), length(idx))
  lapply(seq_len(length(b) - 1L), function(k) idx[(b[k] + 1L):b[k + 1L]])
}

build_heat <- function(cmp) {
  m <- read_sc(cmp)
  miss <- setdiff(pwL, rownames(m))
  if (length(miss)) stop("Pathways missing from matrix: ", paste(miss, collapse = ", "))
  g <- grp_of(colnames(m))
  if (any(table(g) == 0)) stop("Cannot plot empty groups.")

  ord <- order(as.integer(g), -apply(m, 2, stats::median))
  m <- m[, ord, drop = FALSE]; g <- droplevels(g[ord])
  z <- row_z(m)[pwL, , drop = FALSE]

  gm <- function(rw) vapply(levels(g), function(k) mean(z[rw, g == k]), 0)

  up <- vapply(pwL, function(rw) {
    v <- sort(gm(rw), decreasing = TRUE)
    if (v[1] - v[2] >= DELTA) names(v)[1] else ""
  }, character(1))
  dn <- vapply(pwL, function(rw) {
    v <- sort(gm(rw))
    isTRUE(names(v)[1] == "RNHI-low AML" && v[2] - v[1] >= DELTA)
  }, logical(1))

  fLv  <- unique(unname(fblocks[pwL]))
  mrk  <- setNames(match(fblocks[pwL], fLv), pwL)
  eff  <- function(rr, k) vapply(rr, function(r) mean(z[r, g == k]), 0)
  rowseq <- character(0)

  if (ROW_GROUPING == "group") {

    for (k in levels(g)) {
      rr <- pwL[up[pwL] == k]
      if (length(rr)) rowseq <- c(rowseq, rr[order(mrk[rr], -eff(rr, k))])
    }
    rr <- pwL[up[pwL] == "" & dn[pwL]]
    if (length(rr)) rowseq <- c(rowseq, rr[order(mrk[rr], eff(rr, "RNHI-low AML"))])
    rr <- pwL[up[pwL] == "" & !dn[pwL]]
    if (length(rr)) rowseq <- c(rowseq, rr[order(mrk[rr])])
  } else {

    for (fb in fLv) {
      rws <- pwL[fblocks[pwL] == fb]
      for (k in levels(g)) {
        rr <- rws[up[rws] == k]
        if (length(rr)) rowseq <- c(rowseq, rr[order(-eff(rr, k))])
      }
      rr <- rws[up[rws] == "" & dn[rws]]
      if (length(rr)) rowseq <- c(rowseq, rr[order(eff(rr, "RNHI-low AML"))])
      rr <- rws[up[rws] == "" & !dn[rws]]
      if (length(rr)) rowseq <- c(rowseq, rr)
    }
  }
  rowseq <- unique(rowseq)
  stopifnot(setequal(rowseq, pwL), length(rowseq) == length(pwL))

  z <- z[rowseq, , drop = FALSE]
  up <- up[rowseq]; dn <- dn[rowseq]

  if (ROW_GROUPING == "group") {
    lv  <- c(sprintf("Higher in\n%s", gshort[levels(g)]), "Lower in\nRNHI-low", "Not\nspecific")
    fil <- c(vapply(gcol[levels(g)], tint, ""), tint(gcol[["RNHI-low AML"]], 0.10), "grey94")
    names(fil) <- lv
    lab <- ifelse(up != "", sprintf("Higher in\n%s", gshort[up]),
                  ifelse(dn, "Lower in\nRNHI-low", "Not\nspecific"))
    rsplit <- droplevels(factor(lab, levels = lv))
  } else {
    fil <- fcols; rsplit <- factor(fblocks[rowseq], levels = fLv)
  }
  rlv <- levels(rsplit)

  LINE_MM  <- LABEL_PT * LABEL_LH / 72 * 25.4
  max_ln   <- max(vapply(pwLab[rowseq], n_lines, 1L))
  row_mm   <- max(max_ln * LINE_MM, ROW_MIN_MM) + 0.5
  body_mm  <- length(rowseq) * row_mm + 1.4 * (length(rlv) - 1)
  need_mm  <- body_mm + 29
  cat(sprintf("[geom] %s: line %.2f mm x maximum %d lines; row %.2f mm, heatmap %.1f mm, panel %.1f mm\n",
              cmp, LINE_MM, max_ln, row_mm, body_mm, need_mm))

  list(z = z, g = g, up = up, dn = dn, rowseq = rowseq, cmp = cmp,
       rsplit = rsplit, rlv = rlv, rfill = fil[rlv],
       mod = factor(fblocks[rowseq], levels = fLv),
       lab = LAB_FULL[rowseq], need_mm = need_mm)
}

COL_FUN <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

make_legends <- function() {
  lgd_cb <- Legend(col_fun = COL_FUN, title = "ssGSEA z-score", at = c(-2, -1, 0, 1, 2),
                   title_gp = gpar(fontsize = 6.5, fontface = "bold"),
                   labels_gp = gpar(fontsize = 6), direction = "horizontal",
                   title_position = "lefttop", legend_width = unit(26, "mm"),
                   grid_height = unit(2.6, "mm"))
  lgd_fr <- Legend(title = "Outlined block", title_gp = gpar(fontsize = 6.5, fontface = "bold"),
                   labels = c("higher in that group", "lower in RNHI-low"),
                   labels_gp = gpar(fontsize = 6), type = "lines", nrow = 1,
                   grid_width = unit(7, "mm"), grid_height = unit(3, "mm"),
                   legend_gp = gpar(col = c("grey25", gcol[["RNHI-low AML"]]),
                                    lty = c("solid", "22"), lwd = 1.6))
  lst <- list(lgd_cb, lgd_fr)
  if (ROW_GROUPING == "group")
    lst <- c(lst, list(Legend(title = "Functional module",
                              title_gp = gpar(fontsize = 6.5, fontface = "bold"),
                              labels = names(fcols), labels_gp = gpar(fontsize = 6),
                              legend_gp = gpar(fill = fcols), nrow = 1,
                              grid_width = unit(3, "mm"), grid_height = unit(3, "mm"))))
  packLegend(list = lst, direction = "horizontal", column_gap = unit(6, "mm"),
             row_gap = unit(1.5, "mm"), max_width = unit(178, "mm"))
}

draw_panelA <- function(H) {
  z <- H$z; g <- H$g; rlv <- H$rlv

  nlab <- sprintf("%s\n(n=%d)", gshort[levels(g)], as.integer(table(g)))
  names(nlab) <- levels(g)

  ha <- HeatmapAnnotation(
    Group = g, col = list(Group = gcol[levels(g)]),
    simple_anno_size = unit(2.6, "mm"), show_annotation_name = FALSE,
    show_legend = FALSE, border = TRUE)

  la <- if (ROW_GROUPING == "group")
    rowAnnotation(Module = H$mod, col = list(Module = fcols),
                  simple_anno_size = unit(3, "mm"), show_annotation_name = FALSE,
                  show_legend = FALSE, border = TRUE) else NULL

  hm <- Heatmap(
    z, name = "z", col = COL_FUN, show_heatmap_legend = FALSE,
    na_col = "grey90", border = TRUE, border_gp = gpar(col = "black", lwd = 0.5),

    column_split = g, cluster_columns = FALSE, show_column_names = FALSE,
    column_title = nlab[levels(g)], column_title_gp = gpar(fontsize = 6.5),
    column_gap = unit(1.6, "mm"),

    row_split = H$rsplit, cluster_rows = FALSE, row_order = seq_len(nrow(z)),
    row_names_side = "left", row_labels = H$lab,
    row_names_gp = gpar(fontsize = LABEL_PT, lineheight = LABEL_LH),
    row_title = rlv, row_title_side = "right", row_title_rot = 0,
    row_title_gp = gpar(fontsize = 6, fill = H$rfill, col = "black",
                        border = "grey40"),
    row_gap = unit(1.4, "mm"),
    top_annotation = ha, left_annotation = la)

  ht <- draw(hm, annotation_legend_list = list(make_legends()),
             annotation_legend_side = "bottom",
             padding = unit(c(1, 2.5, 1, 2), "mm"))

  ro <- row_order(ht); co <- column_order(ht)
  if (!is.list(ro)) ro <- list(ro)
  if (!is.list(co)) co <- list(co)
  ncs <- length(co)
  cgrp <- vapply(co, function(ix) as.character(unique(g[ix]))[1], character(1))

  seg <- function(x0, y0, x1, y1, gp) grid::grid.segments(x0, y0, x1, y1, gp = gp)

  for (si in seq_along(ro)) {
    lab <- rownames(z)[ro[[si]]]
    n   <- length(lab)
    for (tg in c(levels(g), "lowdown")) {
      pos <- if (tg == "lowdown") {

        if (FRAME_MODE == "block") which(H$dn[lab]) else which(H$dn[lab] & H$up[lab] == "")
      } else which(H$up[lab] == tg)
      if (!length(pos)) next
      col <- if (tg == "lowdown") gcol[["RNHI-low AML"]] else gcol[[tg]]
      lty <- if (tg == "lowdown") "22" else "solid"
      lwd <- 1.5

      cj <- if (FRAME_MODE == "block")
              which(cgrp == if (tg == "lowdown") "RNHI-low AML" else tg) else seq_len(ncs)
      for (run in runs_of(pos)) {
        r1 <- min(run); r2 <- max(run)
        for (j in cj) {
          first <- j == min(cj); last <- j == max(cj)
          decorate_heatmap_body("z", row_slice = si, column_slice = j, {
            pad <- unit(lwd / 2 + 0.2, "pt")
            yb <- unit(1 - r2 / n, "npc") + pad
            yt <- unit(1 - (r1 - 1) / n, "npc") - pad
            xl <- unit(0, "npc") + pad; xr <- unit(1, "npc") - pad
            gp <- gpar(col = col, lwd = lwd, lty = lty, fill = NA, lineend = "square")
            seg(xl, yt, xr, yt, gp); seg(xl, yb, xr, yb, gp)
            if (first) seg(xl, yb, xl, yt, gp)
            if (last)  seg(xr, yb, xr, yt, gp)
          })
        }
      }
    }
  }
  invisible(ht)
}

heat_hm <- function(cmp, fname, wd = W_FULL, ht = NULL) {
  H <- build_heat(cmp)
  if (is.null(ht)) ht <- MM(H$need_mm)
  open_pdf(file.path(res, paste0(fname, ".pdf")), wd, ht); quiet_mbcs(draw_panelA(H)); dev.off()
  open_png(file.path(res, paste0(fname, ".png")), wd, ht); quiet_mbcs(draw_panelA(H)); dev.off()
  cat("[fig]", fname, sprintf("(%.1f x %.1f mm)\n", wd * 25.4, H$need_mm))
  for (k in levels(H$g))
    cat("   ", k, "UP:", paste(pwLab[names(H$up)[H$up == k]], collapse = ", "), "\n")
  cat("    RNHI-low DOWN:", paste(pwLab[names(H$dn)[H$dn]], collapse = ", "), "\n")
  invisible(H)
}

HA <- heat_hm("Normal_EarlyMyeloidProg", "Figure6D_panelA_heatmap_EarlyMyeloidProg")
heat_hm("Normal_MyeloidProgenitor", "Figure6D_heatmap_MyeloidProg")

sel <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
         "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_P53_PATHWAY",
         "HALLMARK_E2F_TARGETS", "HALLMARK_GLYCOLYSIS")

SEL_NEW <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
             "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
SEL_OLD <- sel
sel <- c(SEL_OLD, SEL_NEW)
stopifnot(all(sel %in% pwL), !anyDuplicated(sel))

order_sel_by_group <- function(m, g, sel) {
  z  <- row_z(m)[sel, , drop = FALSE]
  lv <- levels(g)
  gm <- t(vapply(sel, function(pw) vapply(lv, function(k) mean(z[pw, g == k]), 0), numeric(length(lv))))
  colnames(gm) <- lv
  grank   <- apply(gm, 1, function(v) { o <- order(-v); if (v[o[1]] - v[o[2]] >= DELTA) o[1] else length(lv) + 1L })
  modrank <- match(fblocks[sel], unique(unname(fblocks[pwL])))
  eff     <- vapply(seq_along(sel), function(i) if (grank[i] <= length(lv)) gm[i, grank[i]] else 0, 0)
  ord     <- order(grank, modrank, -eff)
  attr(ord, "grank") <- grank
  ord
}

m <- read_sc("Normal_EarlyMyeloidProg"); g <- grp_of(colnames(m))

.oA <- order_sel_by_group(m, g, SEL_OLD)
.oB <- order_sel_by_group(m, g, SEL_NEW)
SEL_A <- SEL_OLD[.oA]
SEL_B <- SEL_NEW[.oB]
sel   <- c(SEL_A, SEL_B)
.gA <- attr(.oA, "grank"); .gB <- attr(.oB, "grank")
.lv <- levels(g)
cat("[panelB] three-row layout: upper block", length(SEL_A), " pathways (", NCOL_B, " columns -> 4+2), lower block",
    length(SEL_B), " pathways (1 row x", NCOL_B, " columns)\n")
for (i in seq_along(SEL_A))
  cat(sprintf("   Upper block %d. %-38s <- %s\n", i, pwLab[SEL_A[i]],
              if (.gA[i] <= length(.lv)) .lv[.gA[i]] else "not specific"))
for (i in seq_along(SEL_B))
  cat(sprintf("   Lower block %d. %-38s <- %s\n", i, pwLab[SEL_B[i]],
              if (.gB[i] <= length(.lv)) .lv[.gB[i]] else "not specific"))

nvec <- table(droplevels(g))
leg_lab <- sprintf("%s (n = %d)", names(nvec), as.integer(nvec))

tt <- fread(file.path(res_in, "Figure6D_pathway_tests.csv"))
pv <- tt[comp == "Normal_EarlyMyeloidProg" & pathway %in% sel]

pexpr <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) "P < 0.001" else sprintf("P = %s", formatC(p, format = "f", digits = 3))
}

pairs_def <- list(list("p_HvHi_bh", 1, 2), list("p_HiLo_bh", 2, 3), list("p_HvLo_bh", 1, 3))

thm <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        strip.background = element_blank(),

        strip.text = element_text(size = 6.5, face = "plain", hjust = 0, margin = margin(b = 1.5)),
        panel.spacing.x = unit(3, "mm"), panel.spacing.y = unit(3.5, "mm"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6.5), legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = -2), plot.margin = margin(2, 3, 1, 2))

make_panelB <- function(pws, ncol = NCOL_B, ytitle = NULL) {
  dd <- rbindlist(lapply(pws, function(pw)
    data.table(pathway = pw, sample = colnames(m), score = m[pw, ], grp = g)))
  dd[, pathway := factor(pathway, levels = pws)]
  dd[, pwlab   := factor(pwLabB[as.character(pathway)], levels = pwLabB[pws])]

  bk <- rbindlist(lapply(pws, function(pw) {
    r <- pv[pathway == pw]; if (!nrow(r)) return(NULL)
    y <- dd[pathway == pw, range(score)]; rg <- diff(y); if (rg == 0) rg <- 1
    k <- 0L
    out <- rbindlist(lapply(pairs_def, function(p) {
      val <- suppressWarnings(as.numeric(r[[p[[1]]]][1]))
      if (!is.finite(val) || val >= 0.05) return(NULL)
      k <<- k + 1L
      data.table(pathway = pw, x = p[[2]], xend = p[[3]],
                 y = y[2] + rg * (0.10 + 0.135 * (k - 1)), lab = pexpr(val))
    }), fill = TRUE)
    out
  }), fill = TRUE)
  if (nrow(bk)) {
    bk[, pwlab := factor(pwLabB[as.character(pathway)], levels = pwLabB[pws])]
    bk[, tick := 0]
    bk[, tick := (max(y) - min(y)) * 0.03, by = pwlab]
    bk[tick == 0, tick := y * 0.012]
  }

  hr <- dd[, .(score = max(score)), by = pwlab]
  if (nrow(bk)) hr <- rbind(hr,
    bk[, .(score = max(y) + (max(y) - min(y)) * 0.10 + diff(range(dd$score)) * 0.02), by = pwlab])

  ns_path <- setdiff(pws, unique(bk$pathway))

  nsd <- if (length(ns_path)) dd[pathway %in% ns_path, .(y = max(score)), by = pwlab] else dd[0L]
  if (nrow(nsd)) nsd[, pwlab := factor(as.character(pwlab), levels = levels(dd$pwlab))]
  cat(sprintf("[panelB] block (%d pathways) BH significant: %s | n.s.: %s\n", length(pws),
              if (nrow(bk)) paste(pwLab[unique(bk$pathway)], collapse = " / ") else "(none)",
              if (length(ns_path)) paste(pwLab[ns_path], collapse = " / ") else "(none)"))

  p <- ggplot(dd, aes(grp, score)) +
    geom_violin(aes(fill = grp), scale = "width", width = 0.85,
                alpha = 0.25, linewidth = 0.25, colour = "grey45") +
    geom_boxplot(aes(colour = grp), width = 0.24, outlier.shape = NA,
                 fill = "white", linewidth = 0.35, fatten = 1.6) +
    geom_jitter(aes(colour = grp), position = position_jitter(width = 0.10, height = 0, seed = 7),
                size = 0.7, alpha = 0.85, stroke = 0) +
    geom_blank(data = hr, aes(x = 1, y = score)) +
    facet_wrap(~pwlab, ncol = ncol, scales = "free_y") +
    scale_fill_manual(values = gcol, labels = leg_lab, drop = FALSE) +
    scale_colour_manual(values = gcol, guide = "none") +
    scale_x_discrete(labels = gshort2) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.04))) +
    labs(x = NULL, y = ytitle) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.9, colour = NA))) +
    thm
  if (nrow(bk)) p <- p +
    geom_segment(data = bk, aes(x = x, xend = xend, y = y, yend = y),
                 linewidth = 0.28, inherit.aes = FALSE) +
    geom_segment(data = bk, aes(x = x, xend = x, y = y - tick, yend = y),
                 linewidth = 0.28, inherit.aes = FALSE) +
    geom_segment(data = bk, aes(x = xend, xend = xend, y = y - tick, yend = y),
                 linewidth = 0.28, inherit.aes = FALSE) +
    geom_text(data = bk, aes(x = (x + xend) / 2, y = y, label = lab),
              parse = FALSE, size = 1.85, vjust = -0.35, inherit.aes = FALSE)
  if (nrow(nsd)) p <- p +
    geom_text(data = nsd, aes(x = 2, y = y, label = "n.s."),
              size = 1.85, colour = "grey35", vjust = -0.30, inherit.aes = FALSE)
  p
}

pB_top <- make_panelB(SEL_A, ytitle = "ssGSEA score")
pB_bot <- make_panelB(SEL_B, ytitle = NULL)
pB <- (pB_top / pB_bot) + plot_layout(heights = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

if (SHOW_TITLE) pB <- pB +
  plot_annotation(title = "Functional state of residual normal HSC/MPP (patient pseudobulk)") &
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))

H_B <- MM(H_B_mm <- 168)
open_pdf(file.path(res, "Figure6D_panelB_pathway_violin.pdf"), W_FULL, H_B); quiet_mbcs(print(pB)); dev.off()
open_png(file.path(res, "Figure6D_panelB_pathway_violin.png"), W_FULL, H_B); quiet_mbcs(print(pB)); dev.off()
cat("[fig] Figure6D_panelB_pathway_violin  (n:",
    paste(names(nvec), as.integer(nvec), sep = "=", collapse = ", "), ")\n")

gA <- quiet_mbcs(grid.grabExpr(draw_panelA(HA), width = W_FULL, height = MM(HA$need_mm),
                               device = dev_null, warn = 0))

gB <- quiet_mbcs(patchworkGrob(pB))
panel_tag <- function(lab) textGrob(lab, x = 0, hjust = 0,
                                    gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT))

draw_composite <- function() {
  grid.newpage()

  lay <- grid.layout(nrow = 4, ncol = 1,
                     heights = unit.c(unit(4, "mm"), unit(HA$need_mm, "null"),
                                      unit(4, "mm"), unit(H_B_mm, "null")))
  pushViewport(viewport(layout = lay, width = 0.97, height = 0.98))
  pushViewport(viewport(layout.pos.row = 1)); grid.draw(panel_tag("a")); popViewport()
  pushViewport(viewport(layout.pos.row = 2)); grid.draw(gA); popViewport()
  pushViewport(viewport(layout.pos.row = 3)); grid.draw(panel_tag("b")); popViewport()
  pushViewport(viewport(layout.pos.row = 4)); grid.draw(gB); popViewport()
  popViewport()
}
H_C_mm <- (4 + HA$need_mm + 4 + H_B_mm) * 1.04
H_C <- MM(H_C_mm)
open_pdf(file.path(res, "Figure6D_composite.pdf"), W_FULL, H_C); quiet_mbcs(draw_composite()); dev.off()
open_png(file.path(res, "Figure6D_composite.png"), W_FULL, H_C); quiet_mbcs(draw_composite()); dev.off()
cat("[fig] Figure6D_composite  (", round(W_FULL * 25.4), "x", round(H_C_mm), "mm )\n")

cat("[done]\n")
