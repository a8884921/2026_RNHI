library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(nnls)
source('rnhi_core.R', local = TRUE)
options(shiny.maxRequestSize = 100 * 1024^2)

# Reference data: the bulk compendium file holds the nine-lineage compositions of the
# 45 healthy reference donors (BoneMarrowMap) and of the AML samples of the bulk compendium.
# It provides the default healthy center and the fixed detection limits for bulk samples.
COMPENDIUM_FILE <- file.path('data', 'bulk_compendium_composition.csv')
compendium <- parse_table(read_input(COMPENDIUM_FILE))
DEFAULT_HEALTH <- compendium$x[compendium$healthy, , drop = FALSE]
COMPENDIUM_DL <- detection_limits(closure(compendium$x[!compendium$healthy, , drop = FALSE]))
N_COMPENDIUM_AML <- sum(!compendium$healthy)
N_DEFAULT_HEALTH <- nrow(DEFAULT_HEALTH)
REFERENCE_BASIS <- readRDS(file.path('data', 'joint_reference_basis.rds'))

fmt <- function(n) format(n, big.mark = ',', trim = TRUE)
EXAMPLES <- list(
  compendium = list(path = COMPENDIUM_FILE, short = 'Bulk AML compendium',
                    label = sprintf('Bulk AML compendium (%s AML transcriptomes with a detected lineage and %d healthy donors)', fmt(N_COMPENDIUM_AML), N_DEFAULT_HEALTH)),
  simcomp    = list(path = file.path('data', 'simulated_composition_demo.csv'), short = 'Simulated compositions',
                    label = 'Simulated compositions (8 samples)'),
  expression = list(path = file.path('data', 'simulated_expression_demo.csv'), short = 'Simulated expression matrix',
                    label = 'Simulated expression matrix (8 samples)'),
  sc         = list(path = file.path('data', 'simulated_scrnaseq_demo.csv'), label = 'Simulated scRNA-seq counts (6 samples)'),
  spatial    = list(path = file.path('data', 'simulated_visium_demo.csv'), label = 'Simulated Visium spots (144 spots)'))
example_key <- function(mode, bulk_example) if (identical(mode, 'bulk')) bulk_example else mode

# Input guide ------------------------------------------------------------------
format_table <- function(headers, rows, caption) {
  tags$div(class = 'example-table-wrap',
    tags$div(class = 'example-caption', caption),
    tags$table(class = 'table table-sm table-bordered example-table',
      tags$thead(tags$tr(lapply(headers, tags$th))),
      tags$tbody(lapply(rows, function(v) tags$tr(lapply(v, tags$td))))))
}

format_examples <- function(mode) {
  switch(mode,
    bulk = tagList(
      format_table(c('sample_id', LINEAGES), list(c('AML_001', '120', '80', '450', '15', '210', '95', '40', '30', '20')),
                   'Composition table: one row per sample (abundances or proportions).'),
      format_table(c('gene', 'Sample_01', 'Sample_02', 'Sample_03'),
                   list(c('A1BG', '12.6', '9.8', '15.1'), c('AAAS', '31.2', '27.5', '35.0'), c('ABL1', '18.7', '21.3', '16.9')),
                   'Expression matrix: gene symbols, one column per sample (linear scale, e.g., counts or TPM).')),
    sc = format_table(c('sample_id', 'celltype', 'umap_1', 'umap_2'),
                      list(c('S01', 'Mono', '1.42', '-0.73'), c('S01', 'B', '-2.10', '3.05'), c('S02', 'NK', '0.37', '2.88')),
                      'Cell metadata: one row per annotated cell (umap_1 and umap_2 optional).'),
    spatial = format_table(c('spot_id', 'cell_type', 'abundance', 'x', 'y'),
                           list(c('spot_001', 'MyeProg', '25', '1', '1'), c('spot_001', 'Mono', '15', '1', '1'),
                                c('spot_002', 'Gran', '45', '2', '1')),
                           'Spot composition: one row per spot and lineage (x and y optional).'))
}

mode_text <- c(
  bulk = 'Upload a lineage composition table (wide or long) or a gene-by-sample expression matrix. Expression data are deconvolved automatically (MuSiC estimator).',
  sc = 'Upload cell metadata (sample_id or orig.ident, celltype) or a long composition table. Cells are counted per sample; supplied UMAP coordinates are plotted.',
  spatial = 'Upload spot metadata (spot_id or barcode, celltype) or a long composition table (spot_id, cell_type, abundance). Supplied coordinates (x/y, imagecol/imagerow or pxl_col_in_fullres/pxl_row_in_fullres) are used for the spatial map.')

instructions <- function(mode) {
  title <- c(bulk = 'Bulk RNA / pseudobulk', sc = 'scRNA-seq', spatial = 'Visium')[[mode]]
  tagList(
    card(class = 'hero-card', card_header(paste(title, 'input')),
      div(class = 'rnhi-hero',
        div(class = 'eyebrow', 'AML · residual normal hematopoiesis'),
        h1('RNHI'),
        h3('Residual normal hematopoiesis integrity'),
        p('Measures how closely the lineage composition of residual normal hematopoiesis resembles healthy bone marrow. Only the lineage composition is used.')),
      div(class = 'workflow',
        div('01 · Prepare', p('Composition or expression table')),
        div('02 · Upload', p('CSV or TSV')),
        div('03 · Calculate', p('RNHI for each sample')),
        div('04 · Explore', p('Scores, figures, downloads'))),
      p(mode_text[[mode]]),
      format_examples(mode),
      tags$ul(class = 'guide-notes',
        tags$li(HTML(paste0('<strong>Lineages:</strong> <code>', paste(LINEAGES, collapse = ', '), '</code>'))),
        tags$li(HTML(paste0('<strong>Healthy reference:</strong> ', N_DEFAULT_HEALTH, ' BoneMarrowMap donors by default; mark your own healthy samples with <code>healthy_reference = 1</code> or upload them separately.'))),
        tags$li(HTML('<strong>Zeros:</strong> replaced with 0.65 × the lineage detection limit; by default, fixed compendium limits for bulk data and limits estimated from the upload for scRNA-seq and Visium.')))),
    card(card_header('Output'),
      p('RNHI is the negative Aitchison distance to the healthy center (RNHI ≤ 0; values closer to 0 indicate a composition closer to healthy bone marrow).'),
      p('Score table, RNHI distribution, lineage composition, lineage profile versus healthy and, when coordinates are supplied, a UMAP or spatial map; CSV, PDF and PNG downloads.')))
}

method_panel <- function() {
  card(card_header('Method'),
    tags$ol(
      tags$li(HTML('<strong>Closure:</strong> the nine lineage abundances are rescaled to sum to 1.')),
      tags$li(HTML('<strong>Zero replacement:</strong> zeros are replaced with 0.65 × the detection limit of the lineage (its minimum positive proportion), and non-zero parts are rescaled (following the multiplicative replacement rule of zCompositions::multRepl).')),
      tags$li(HTML('<strong>CLR transformation:</strong> clr(x) = log(x) − mean(log(x)).')),
      tags$li(HTML(sprintf('<strong>Healthy center c:</strong> mean CLR vector of the healthy reference donors (default: %d BoneMarrowMap donors).', N_DEFAULT_HEALTH))),
      tags$li(HTML('<strong>Score:</strong> RNHI(x) = −‖clr(x) − c‖<sub>2</sub>, the negative Aitchison distance (RNHI ≤ 0).'))),
    p(HTML('<strong>Lineages:</strong> MyeProg, Mono, Gran, Mega, EarlyEry, LateEry, B, NK and DC. HSC/MPP, LMPP and T cells are excluded because they are detected in fewer than 4% of bulk samples.')),
    p(HTML('<strong>Detection limits:</strong> <em>Fixed</em> (default for bulk data) uses the limits of the bulk AML compendium (3,377 transcriptomes, 11 public datasets), so each sample is scored independently. <em>Estimated</em> (default for scRNA-seq and Visium data) uses the minimum positive proportion of each lineage among the uploaded query samples, so scores depend on the uploaded set. Healthy reference samples use limits estimated within the reference.')),
    p(HTML('<strong>Expression input:</strong> each profile is scaled to relative abundance over the reference genes, genes with zero counts are removed, and cell-state proportions are estimated with the MuSiC weighted estimator against the basis matrix of the 19-state joint healthy/AML reference (7,000 genes). At least 10% of the reference genes must be present, and at least 20% of these must be expressed in each sample.')),
    p(HTML('<strong>Quality flag:</strong> samples with a single detected lineage are flagged (qc_flag).')),
    p('No survival, risk group or mutation data are used.'))
}

about_panel <- function() {
  card(card_header('About'),
    p('The RNHI Calculator computes residual normal hematopoiesis integrity (RNHI) for AML samples from bulk expression or composition data, scRNA-seq and Visium data.'),
    p(HTML('Source code: <a href="https://github.com/a8884921/2026_RNHI" target="_blank">https://github.com/a8884921/2026_RNHI</a>')),
    p('Uploaded files are processed only within the current session and are not stored. For research use only.'))
}

# User interface ---------------------------------------------------------------
ui <- page_navbar(
  title = span(class = 'brand-copy', span('RNHI Calculator', class = 'brand-title'),
               span('Residual normal hematopoiesis integrity', class = 'brand-subtitle')),
  window_title = 'RNHI Calculator',
  fillable = FALSE,
  header = tags$head(includeCSS(file.path('www', 'style.css'))),
  theme = bs_theme(version = 5, bootswatch = 'flatly', primary = '#176b87', secondary = '#ee6c4d'),
  nav_panel('Calculate', layout_sidebar(fillable = FALSE,
    sidebar = sidebar(width = 320,
      radioButtons('mode', 'Input type', c('Bulk RNA / pseudobulk' = 'bulk', 'scRNA-seq' = 'sc', 'Visium' = 'spatial')),
      conditionalPanel("input.mode == 'bulk'",
        selectInput('bulk_example', 'Example data',
                    choices = setNames(c('compendium', 'simcomp', 'expression'),
                                       c(EXAMPLES$compendium$short, EXAMPLES$simcomp$short, EXAMPLES$expression$short)))),
      actionButton('use_example', 'Load example', class = 'btn-outline-primary w-100'),
      fileInput('file', 'Or upload a table', accept = c('.csv', '.tsv', '.txt')),
      textOutput('source_label'),
      radioButtons('dl_mode', 'Detection limits',
                   c('Fixed (bulk AML compendium)' = 'fixed', 'Estimated from uploaded samples' = 'estimated'),
                   selected = 'fixed'),
      fileInput('health', 'Healthy reference (optional)', accept = c('.csv', '.tsv', '.txt')),
      actionButton('go', 'Calculate RNHI', class = 'btn-primary w-100'),
      actionButton('back', 'Input guide', class = 'w-100'),
      actionButton('reset_app', 'Reset', class = 'btn-outline-secondary w-100'),
      downloadButton('example_download', 'Download example', class = 'w-100'),
      helpText('Uploaded files are processed only within the current session and are not stored.')),
    uiOutput('error_panel'), uiOutput('content'))),
  nav_panel('Method', method_panel()),
  nav_panel('About', about_panel()))

# Server -----------------------------------------------------------------------
server <- function(input, output, session) {
  result <- reactiveVal(NULL)
  error <- reactiveVal(NULL)
  source_kind <- reactiveVal('none')
  clear <- function() { result(NULL); error(NULL) }
  current_example <- reactive(EXAMPLES[[example_key(input$mode, input$bulk_example)]])

  observeEvent(input$file, { source_kind(if (is.null(input$file)) 'none' else 'upload'); clear() })
  observeEvent(input$use_example, {
    source_kind('example'); clear()
    showNotification(paste('Example loaded:', current_example()$label), type = 'message', duration = 5)
  })
  observeEvent(input$bulk_example, clear(), ignoreInit = TRUE)
  observeEvent(input$mode, {
    clear()
    if (identical(source_kind(), 'example')) source_kind('none')
    updateRadioButtons(session, 'dl_mode', selected = if (identical(input$mode, 'bulk')) 'fixed' else 'estimated')
  }, ignoreInit = TRUE)
  observeEvent(input$dl_mode, clear(), ignoreInit = TRUE)
  observeEvent(input$health, clear(), ignoreInit = TRUE)
  observeEvent(input$back, clear())
  observeEvent(input$reset_app, session$reload())

  output$source_label <- renderText({
    switch(source_kind(),
      none = 'Load an example or upload a file.',
      example = paste('Selected:', current_example()$label),
      upload = paste('Selected:', input$file$name))
  })

  observeEvent(input$go, {
    clear()
    tryCatch(withProgress(message = 'Calculating RNHI', value = .15, {
      if (source_kind() == 'none') stop('Load an example or upload a table first.')
      mode <- input$mode
      dat <- if (source_kind() == 'example') read_input(current_example()$path) else read_input(input$file$datapath, input$file$name)
      raw_dat <- dat
      custom <- if (!is.null(input$health)) read_input(input$health$datapath, input$health$name) else NULL
      incProgress(.25)
      composition_source <- 'composition table'
      decon_qc <- NULL
      if (mode == 'bulk' && is_expression_matrix(dat)) {
        dc <- deconvolve_expression(dat, REFERENCE_BASIS)
        dat <- dc$data; decon_qc <- dc$qc
        composition_source <- sprintf('MuSiC deconvolution, %s of %s reference genes matched', fmt(dc$reference_genes), fmt(dc$total_reference_genes))
      }
      if (mode %in% c('sc', 'spatial')) {
        dat <- prepare_modality_table(dat, mode)
        composition_source <- if (identical(dat, raw_dat)) 'composition table' else 'annotated cells counted per lineage'
      }
      if (!is.null(custom) && mode == 'bulk' && is_expression_matrix(custom)) custom <- deconvolve_expression(custom, REFERENCE_BASIS)$data
      if (!is.null(custom) && mode %in% c('sc', 'spatial')) custom <- prepare_modality_table(custom, mode)
      incProgress(.3)
      r <- calculate_rnhi(dat, DEFAULT_HEALTH, custom, query_dl = if (identical(input$dl_mode, 'fixed')) COMPENDIUM_DL else NULL)
      r$mode <- mode
      r$source <- if (source_kind() == 'example') current_example()$label else input$file$name
      r$composition_source <- composition_source
      r$decon_qc <- decon_qc
      r$health_source <- if (any(r$data$healthy_reference == 1)) 'samples flagged healthy_reference = 1'
                         else if (!is.null(custom)) input$health$name
                         else sprintf('%d BoneMarrowMap donors (default)', N_DEFAULT_HEALTH)
      r$plot_data <- if (mode == 'sc') sc_umap_data(r, raw_dat) else if (mode == 'spatial') spatial_data(r, raw_dat) else NULL
      incProgress(.3); result(r)
    }), error = function(e) { error(conditionMessage(e)); showNotification(conditionMessage(e), type = 'error', duration = 12) })
  })

  output$error_panel <- renderUI({
    if (!is.null(error())) div(class = 'alert alert-danger', role = 'alert', strong('Calculation failed: '), error())
  })

  plot_card <- function(title, id, height = '400px', note = NULL)
    card(card_header(title), if (!is.null(note)) p(class = 'card-note', note), plotOutput(id, height = height),
         div(class = 'download-row', downloadButton(paste0(id, '_pdf'), 'PDF'), downloadButton(paste0(id, '_png'), 'PNG')))

  output$content <- renderUI({
    if (is.null(result())) {
      notice <- if (source_kind() == 'example')
        div(class = 'alert alert-success example-selected', strong('Example loaded: '), current_example()$label,
            p('Click Calculate RNHI to score it.'))
      return(tagList(notice, instructions(if (is.null(input$mode)) 'bulk' else input$mode)))
    }
    r <- result()
    n_single <- sum(r$data$qc_flag != 'ok')
    mode_card <- if (is.null(r$plot_data)) NULL else switch(r$mode,
      sc = plot_card('UMAP', 'umap', '520px', 'Cells colored by the RNHI of their sample.'),
      spatial = plot_card('Spatial RNHI map', 'spatial', '520px', 'Spots colored by RNHI.'),
      NULL)
    n_poor <- if (is.null(r$decon_qc)) 0 else sum(grepl('poor fit', r$decon_qc$status))
    tagList(
      div(id = 'results_ready', class = 'alert alert-success', 'Calculation complete.'),
      card(card_header('RNHI scores'),
           div(class = 'result-summary', uiOutput('summary')),
           if (n_single > 0) div(class = 'alert alert-warning',
             sprintf('%d %s a single detected lineage (see qc_flag).', n_single, if (n_single == 1) 'sample has' else 'samples have')),
           if (n_poor > 0) div(class = 'alert alert-warning',
             sprintf('Poor deconvolution fit (R² < 0) for %d %s; check that genes are symbols and values are on a linear scale.',
                     n_poor, if (n_poor == 1) 'sample' else 'samples')),
           DTOutput('table', fill = FALSE), div(class = 'download-row', downloadButton('download', 'Download CSV'))),
      if (!is.null(r$decon_qc)) card(card_header('Deconvolution'), DTOutput('qc_table', fill = FALSE),
                                     div(class = 'download-row', downloadButton('download_qc', 'Download CSV'))),
      plot_card('RNHI distribution', 'dist', '360px'),
      plot_card('Lineage composition', 'heat', '400px', 'Within-sample lineage proportions; samples ordered by RNHI.'),
      plot_card('Lineage profile versus healthy', 'profile', '400px',
                'CLR difference from the healthy center (red: relatively enriched; blue: relatively depleted).'),
      mode_card)
  })

  output$summary <- renderUI({
    r <- result(); req(r)
    nq <- sum(r$data$healthy_reference == 0); nh <- sum(r$data$healthy_reference == 1)
    unit <- if (r$mode == 'spatial') 'spots' else 'samples'
    tags$ul(
      tags$li(sprintf('Scored: %s %s (%s query, %s healthy reference)', fmt(r$n), unit, fmt(nq), fmt(nh))),
      tags$li(paste0('Input: ', r$source, '; ', r$composition_source)),
      tags$li(paste0('Healthy reference: ', r$health_source)),
      tags$li(paste0('Detection limits: ', if (r$dl_source == 'fixed') 'fixed (bulk AML compendium)' else 'estimated from uploaded query samples')))
  })
  round_cols <- function(d) setdiff(names(d)[vapply(d, is.double, logical(1))], 'lineage_total')
  output$table <- renderDT({
    r <- result(); req(r)
    datatable(r$data, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) |>
      formatRound(round_cols(r$data), digits = 4)
  })
  output$qc_table <- renderDT({
    r <- result(); req(r, r$decon_qc)
    datatable(r$decon_qc, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) |>
      formatRound('r_squared', digits = 4)
  })
  output$download <- downloadHandler(filename = function() 'RNHI_scores.csv',
    content = function(file) { req(result()); utils::write.csv(result()$data, file, row.names = FALSE) })
  output$download_qc <- downloadHandler(filename = function() 'RNHI_deconvolution.csv',
    content = function(file) { req(result()$decon_qc); utils::write.csv(result()$decon_qc, file, row.names = FALSE) })
  output$example_download <- downloadHandler(
    filename = function() basename(current_example()$path),
    content = function(file) file.copy(current_example()$path, file, overwrite = TRUE))

  plots <- list(dist = distribution_plot, heat = composition_plot, profile = profile_plot,
                umap = sc_umap_plot, spatial = spatial_plot)
  sizes <- list(dist = c(10, 4.5), heat = c(10, 4.5), profile = c(10, 4.5), umap = c(8, 6), spatial = c(8, 6))
  for (k in names(plots)) local({
    key <- k
    output[[key]] <- renderPlot({ r <- result(); req(r); if (key %in% c('umap', 'spatial')) req(r$plot_data); plots[[key]](r) }, res = 120)
    for (ext in c('pdf', 'png')) local({
      e <- ext
      output[[paste(key, e, sep = '_')]] <- downloadHandler(
        filename = function() paste0('RNHI_', key, '.', e),
        content = function(file) {
          r <- result(); req(r)
          ggsave(file, plots[[key]](r), device = e, width = sizes[[key]][1], height = sizes[[key]][2], dpi = 300, bg = 'white')
        })
    })
  })
}

shinyApp(ui, server)
