library(shiny)
library(bslib)
library(ggplot2)
library(DT)
source('rnhi_core.R', local=TRUE)
default_comp <- read_input('composition_A4.csv')
default_parsed <- parse_table(default_comp)
default_health <- default_parsed$x[default_parsed$healthy,,drop=FALSE]
options(shiny.maxRequestSize=100*1024^2)

demo_path <- function(mode, bulk_variant='composition') {
  if (mode == 'bulk' && identical(bulk_variant,'expression')) return('data/bulk_expression_demo.csv')
  if (mode == 'bulk' && identical(bulk_variant,'beataml')) return('data/beatAML2_demo.csv')
  switch(mode, bulk='composition_A4.csv', sc='data/Van_scRNA_demo.csv', spatial='data/Visium_spot_demo.csv', 'composition_A4.csv')
}
demo_label <- function(mode, bulk_variant='composition') {
  p <- demo_path(mode,bulk_variant)
  if (mode == 'bulk' && identical(bulk_variant,'expression')) paste(basename(p),'| MuSiC') else basename(p)
}

example_table <- function(mode, bulk_variant='composition') {
  if (mode == 'bulk' && identical(bulk_variant,'expression')) {
    headers <- c('gene','BeatAML2_expr_01','BeatAML2_expr_02','BeatAML2_expr_03')
    caption <- 'Expression format — first column is a gene symbol; each remaining column is one bulk sample.'
    rows <- list(c('A1BG','194311.7','239293.6','203775.9'),c('AAAS','10234.5','12880.2','11021.7'),c('ABL1','5821.3','4900.6','6310.4'))
  } else if (mode == 'bulk') {
    headers <- c('sample_id',LINEAGES)
    if (identical(bulk_variant,'beataml')) {
      values <- c('BeatAML2_01','173.84','281.747','532.106','416.365','355.474','267.441','225.573','166.094','469.511')
      caption <- 'BeatAML2 composition demo — one row per patient; lineage abundances are shown.'
    } else {
      values <- c('AML_001','120','80','450','15','210','95','40','30','20')
      caption <- 'Wide format — one row per patient (values are lineage abundances or fractions).'
    }
    rows <- list(values)
  } else if (mode == 'sc') {
    headers <- c('sample_id','cell_type','abundance')
    caption <- 'Long format — annotated cells are aggregated by sample_id and cell_type.'
    rows <- list(c('Van_01','MyeProg','120'),c('Van_01','Mono','80'),c('Van_01','Gran','450'),
      c('Van_01','Mega','15'),c('Van_01','EarlyEry','210'),c('Van_01','LateEry','95'),
      c('Van_01','B','40'),c('Van_01','NK','30'),c('Van_01','DC','20'))
  } else {
    headers <- c('spot_id','cell_type','abundance','x','y')
    caption <- 'Long spot format — keep each Visium spot separate and include optional x/y coordinates.'
    rows <- list(c('spot_001','MyeProg','25','1','1'),c('spot_001','Mono','15','1','1'),c('spot_001','Gran','45','1','1'),
      c('spot_001','Mega','2','1','1'),c('spot_001','EarlyEry','12','1','1'),c('spot_001','LateEry','8','1','1'),
      c('spot_001','B','4','1','1'),c('spot_001','NK','3','1','1'),c('spot_001','DC','2','1','1'))
  }
  tags$div(class='example-table-wrap',
    tags$div(class='example-caption',caption),
    tags$table(class='table table-sm table-bordered example-table',
      tags$thead(tags$tr(lapply(headers,tags$th))),
      tags$tbody(lapply(rows,function(v) tags$tr(lapply(v,tags$td))))))
}
expression_example_table <- function() {
  headers <- c('gene','BeatAML2_expr_01','BeatAML2_expr_02','BeatAML2_expr_03')
  rows <- list(c('A1BG','194311.7','239293.6','203775.9'),c('AAAS','10234.5','12880.2','11021.7'),c('ABL1','5821.3','4900.6','6310.4'))
  tags$div(class='example-table-wrap music-example',
    tags$div(class='example-caption','MuSiC expression example — gene-by-sample matrix (first column = gene; each additional column = one bulk sample).'),
    tags$table(class='table table-sm table-bordered example-table',
      tags$thead(tags$tr(lapply(headers,tags$th))),
      tags$tbody(lapply(rows,function(v) tags$tr(lapply(v,tags$td))))))
}

mode_result_items <- function(mode) switch(mode,
  bulk=tags$ul(tags$li('RNHI score table'),tags$li('Healthy/query density plot'),tags$li('Ordered nine-lineage composition heatmap'),tags$li('CSV and publication-ready PDF/PNG export')),
  sc=tags$ul(tags$li('RNHI score table per sample'),tags$li('UMAP overlay; supplied UMAP coordinates are retained, otherwise a composition-derived preview is generated'),tags$li('CSV and UMAP PDF/PNG export')),
  spatial=tags$ul(tags$li('RNHI score table per spot'),tags$li('Visium spatial RNHI map with a print-safe continuous colour scale'),tags$li('Composition-space UMAP overview of all spots'),tags$li('CSV and spatial/UMAP PDF/PNG export')))

instructions <- function(mode, bulk_variant='composition') {
  title <- c(bulk='Bulk RNA / pseudobulk',sc='scRNA-seq',spatial='Visium spot')[[mode]]
  tagList(
    card(class='hero-card',card_header(paste(title,'— input guide')),
      div(class='rnhi-hero',div(class='eyebrow','AML · residual normal hematopoiesis'),h1('RNHI'),h3('Residual Normal Hematopoiesis Index'),p('A single composition-based readout of how closely the remaining normal blood-forming structure resembles healthy marrow. RNHI is independent of survival, ELN and mutation labels.')),
      if(mode=='bulk') p('Composition input: one row per patient. Expression input: one row per gene, with one column per bulk sample; MuSiC deconvolution is automatic.'),
      if(mode=='sc') p('Metadata input: one row per annotated cell. The app aggregates cell counts by sample and celltype before scoring.'),
      if(mode=='spatial') p('Metadata input: one row per annotated spot/cell. The app aggregates cell counts by spot and celltype before scoring.'),
      div(class='workflow',div('01 · Prepare',p('Nine normal lineages')),div('02 · Upload',p('CSV or TSV')),div('03 · Calculate',p('Frozen A4 definition')),div('04 · Explore',p('Scores, figures & export'))),
      p('Example input table:'),example_table(mode,bulk_variant),
      if(mode=='bulk') tagList(div(class='music-note',strong('Automatic bulk recognition. '),'A gene-by-sample table (first column = gene, remaining columns = samples) is detected automatically and sent through the bundled locked MuSiC reference before RNHI. A pre-computed nine-lineage table is also accepted.'),if(!identical(bulk_variant,'expression')) expression_example_table()),
      if(mode=='bulk') p('Composition tables can be wide (nine lineage columns) or long (sample_id, cell_type, abundance); repeated long-format keys are summed.'),
      if(mode=='sc') p('A long composition table with sample_id, cell_type and abundance is also accepted; repeated keys are summed.'),
      if(mode=='spatial') p('A long composition table with spot_id/sample_id, cell_type and abundance is also accepted; repeated keys are summed.'),
      p('Use exact lineage names, finite non-negative values and non-empty sample IDs. Every lineage must be positive in at least one sample within each group used for RNHI; otherwise a detection limit cannot be estimated.'),
      p('Optional healthy_reference = 1 identifies healthy rows (0 = patient). Each sample must have one consistent flag. Healthy and patient zero replacements are calculated separately.'),
      if(mode=='sc') p('Upload a Seurat meta.data export with a sample_id (or orig.ident) column and a celltype column using the nine frozen names. The app counts cells by sample_id × celltype to create abundance automatically; RNHI is computed per sample, not per single cell.'),
      if(mode=='spatial') p('Upload a spot-level metadata export with spot_id (or sample_id), a celltype column and optional x/y or imagecol/imagerow coordinates. The app counts annotated cells per spot × celltype to create abundance automatically.'),
      if(mode=='sc') p('The bundled example is a composition-derived UMAP preview when true UMAP coordinates are not supplied; upload umap_1 and umap_2 to display your own embedding.'),
      if(mode=='spatial') p('The spatial figure uses your coordinates when present; otherwise it lays spots out on a compact grid for a runnable preview.')),
    card(card_header('What the results mean'),
      p('RNHI is the negative Aitchison distance to the healthy CLR center: values closer to zero indicate a more similar composition; more negative values indicate greater structural imbalance.'),
      mode_result_items(mode),
      p('The nine-lineage abundance sum is not the full residual-normal fraction when excluded lineages are absent.'),
      if(mode=='bulk') p('Local validation example: composition_A4.csv, 45 healthy donors + 3,310 AML samples. The bundled expression example is intentionally small and uses the locked MuSiC basis.'),
      if(mode!='bulk') p('The bundled ',if(mode=='sc') 'Van-style' else 'Visium-style',' example is a small, fully runnable demonstration; replace it with your own annotated table for analysis.')))
}
ui <- page_navbar(title=span(class='brand-lockup',
  img(src='logo.png',class='brand-logo',alt='Fujian Medical University and affiliated hospital'),
  span(class='brand-copy',span('AML RNHI',class='brand-title'),span('Residual Normal Hematopoiesis Index',class='brand-subtitle'))), fillable=FALSE,
  header=tags$head(includeCSS(file.path('www','style.css'))),
  theme=bs_theme(version=5,bootswatch='flatly',primary='#176b87',secondary='#ee6c4d'),
  nav_panel('Calculate',layout_sidebar(fillable=FALSE,sidebar=sidebar(width=320,
    radioButtons('mode','Input modality',c('Bulk RNA / pseudobulk'='bulk','scRNA-seq'='sc','Visium spots'='spatial')),
    conditionalPanel("input.mode == 'bulk'",selectInput('bulk_demo','Bulk example',choices=c('Local A4 composition table'='composition','BeatAML2 composition demo'='beataml','Bulk expression matrix (MuSiC)'='expression'),selected='composition')),
    actionButton('use_demo','example',class='btn-outline-primary w-100'),
    fileInput('file','Upload your table',accept=c('.csv','.tsv','.txt')),
    textOutput('source_label'),
    fileInput('health','Optional separate healthy reference',accept=c('.csv','.tsv','.txt')),
    helpText('Frozen A4 fraction: 0.65. Embedded healthy_reference rows replace the bundled reference.'),
    actionButton('go','Calculate RNHI',class='btn-primary w-100'),
    actionButton('back','Input guide',class='w-100'),
    actionButton('refresh_app','Refresh / reset',class='btn-outline-secondary w-100'),
    downloadButton('example_download','Download example input',class='w-100'),
    helpText('Uploads are processed in this session. Local example data contain research sample IDs; public redistribution needs a separate data review.')
  ),uiOutput('error_panel'),uiOutput('content'))),
  nav_panel('Method',card(h2('Frozen RNHI A4'),
    p('A4 is the measurable nine-lineage normal-hematopoiesis subcomposition: MyeProg, Mono, Gran, Mega, EarlyEry, LateEry, B, NK and DC. HSC_MPP, LMPP and T are excluded from this frozen definition because they do not meet the project detection criterion.'),
    p('RNHI = -sqrt(sum((CLR(sample) - mean(CLR(healthy donors)))^2))'),
    p('Close the nine lineage abundances to sum 1. For each lineage, take the minimum positive closed proportion within the patient group or healthy group separately. Replace zeros with 0.65 times that limit, and rescale positive entries to retain sum 1.'),
    p('The healthy center is the equal-donor mean CLR. Do not pool patient and healthy detection limits. No survival cutoff is used.'),
    p('Detection limits depend on the sample set. Upload the same complete set and the same healthy reference to reproduce local scores. Subsetting or switching platforms can change RNHI.'),
    p('A single detected lineage produces a degenerate score. The app flags it while retaining the frozen formula.'))),
  nav_panel('About',card(h2('AML RNHI research calculator'),p('Local A4 implementation validation. Research use only. Raw-expression deconvolution and cross-platform validation are separate tasks.')))
)
server <- function(input,output,session) {
  result <- reactiveVal(NULL)
  error <- reactiveVal(NULL)
  source_kind <- reactiveVal('none')
  observeEvent(input$file, { if (is.null(input$file)) source_kind('none') else source_kind('upload'); result(NULL); error(NULL) })
  observeEvent(input$use_demo, {
    source_kind('demo'); result(NULL); error(NULL)
    showNotification(paste('Example selected:', demo_label(input$mode,input$bulk_demo)), type='message', duration=8)
  })
  observeEvent(input$bulk_demo, { if (identical(source_kind(),'demo')) result(NULL); error(NULL) }, ignoreInit=TRUE)
  observeEvent(input$mode, { result(NULL); error(NULL) },ignoreInit=TRUE)
  observeEvent(input$health, { result(NULL); error(NULL) },ignoreNULL=FALSE,ignoreInit=TRUE)
  observeEvent(input$back, { result(NULL); error(NULL) })
  observeEvent(input$refresh_app, { source_kind('none'); result(NULL); error(NULL) })
  output$source_label <- renderText({
    switch(source_kind(),
      none='Choose example or upload a file.',
      demo=paste('Selected example:',demo_label(input$mode,input$bulk_demo)),
      upload=paste('Selected:',input$file$name))
  })
  observeEvent(input$go, {
    error(NULL); result(NULL)
    tryCatch(withProgress(message='Calculating RNHI',value=.15,{
      if(source_kind()=='none') stop('Choose the local A4 example or upload a table first.')
      dat <- if(source_kind()=='demo') read_input(demo_path(input$mode,input$bulk_demo)) else read_input(input$file$datapath,input$file$name)
      raw_dat <- dat
      custom <- if(!is.null(input$health)) read_input(input$health$datapath,input$health$name) else NULL
      incProgress(.3)
      upstream <- 'Composition table'; music_qc <- NULL
      if (input$mode == 'bulk' && is_expression_matrix(dat)) {
        mu <- music_deconvolve(dat)
        dat <- mu$data; music_qc <- mu$qc
        upstream <- paste0(mu$method,' (',mu$common_genes,' common genes; ',sum(mu$qc$status=='success'),' successful samples)')
      }
      if (input$mode %in% c('sc','spatial')) dat <- prepare_modality_table(dat,input$mode)
      if (!is.null(custom) && is_expression_matrix(custom)) custom <- music_deconvolve(custom)$data
      r <- calculate_rnhi(dat,default_health,custom,frac=.65)
      r$mode <- input$mode
      r$source <- if(source_kind()=='demo') demo_label(input$mode,input$bulk_demo) else input$file$name
      r$upstream <- upstream; r$music_qc <- music_qc
      if (input$mode == 'sc') r$plot_data <- make_sc_umap(r, raw_dat)
      if (input$mode == 'spatial') { r$plot_data <- make_spatial_data(r, raw_dat); r$plot_data_umap <- make_spot_umap(r) }
      incProgress(.4); result(r)
    }),error=function(e) { error(conditionMessage(e)); showNotification(conditionMessage(e),type='error',duration=12) })
  })
  output$error_panel <- renderUI({ if(!is.null(error())) div(class='alert alert-danger',role='alert',strong('Calculation failed: '),error()) })
  output$content <- renderUI({
    if(is.null(result())) {
      notice <- if (source_kind() == 'demo') div(class='alert alert-success example-selected',strong('Example selected: '),demo_label(input$mode,input$bulk_demo),p('Click Calculate RNHI to generate the example result.')) else NULL
      return(tagList(notice,instructions(if(is.null(input$mode)) 'bulk' else input$mode,input$bulk_demo)))
    }
    r <- result(); mode <- r$mode
    mode_card <- switch(mode,
      bulk=tagList(
        card(card_header('RNHI distribution'),plotOutput('dist',height='360px'),downloadButton('dist_pdf','PDF'),downloadButton('dist_png','PNG')),
        card(card_header('Ordered composition profile'),plotOutput('heat',height='400px'),downloadButton('heat_pdf','PDF'),downloadButton('heat_png','PNG'))),
      sc=card(card_header('scRNA-seq UMAP — RNHI overlay'),p('Each point is an annotated cell when UMAP coordinates are supplied; otherwise this is a deterministic composition-derived preview. Colour shows sample-level RNHI.'),plotOutput('umap',height='520px'),downloadButton('umap_pdf','PDF'),downloadButton('umap_png','PNG')),
      spatial=tagList(card(card_header('Visium spot spatial RNHI map'),p('Each tile is a spot; colour shows its RNHI.'),plotOutput('spatial',height='520px'),downloadButton('spatial_pdf','PDF'),downloadButton('spatial_png','PNG')),
        card(card_header('Visium spots — UMAP view'),p('A compact composition-space view of the spot-level RNHI values.'),plotOutput('spatial_umap',height='420px'),downloadButton('spatial_umap_pdf','PDF'),downloadButton('spatial_umap_png','PNG'))))
    tagList(div(id='results_ready',class='alert alert-success','Calculation complete — results below.'),
      card(card_header('RNHI results'),textOutput('message'),DTOutput('table'),downloadButton('download','Download scores CSV')),mode_card)
  })
  output$table <- renderDT({ req(result()); datatable(result()$data,rownames=FALSE,options=list(pageLength=15,scrollX=TRUE)) })
  output$message <- renderText({ r<-result(); req(r); paste('Scored',r$n,'records;',sum(r$data$healthy_reference==0),'query samples;',r$nhealth,'healthy reference donors. Input:',r$source,'Upstream:',r$upstream) })
  output$dist <- renderPlot({ req(result()); distribution_plot(result()) },res=120)
  output$heat <- renderPlot({ req(result()); composition_plot(result()) },res=120)
  output$umap <- renderPlot({ req(result(),result()$mode == 'sc'); sc_umap_plot(result()) },res=120)
  output$spatial <- renderPlot({ req(result(),result()$mode == 'spatial'); spatial_plot(result()) },res=120)
  output$spatial_umap <- renderPlot({ req(result(),result()$mode == 'spatial'); spot_umap_plot(result()) },res=120)
  output$download <- downloadHandler(filename=function() 'RNHI_scores.csv',content=function(file) { req(result()); readr::write_csv(result()$data,file) })
  output$example_download <- downloadHandler(
    filename=function() basename(demo_path(input$mode,input$bulk_demo)),
    content=function(file) file.copy(demo_path(input$mode,input$bulk_demo),file,overwrite=TRUE))
  for (kind in c('dist','heat')) for (ext in c('pdf','png')) local({
    k<-kind; e<-ext
    output[[paste(k,e,sep='_')]] <- downloadHandler(filename=function() paste0('RNHI_',k,'.',e),content=function(file) {
      req(result()); p<-if(k=='dist') distribution_plot(result()) else composition_plot(result())
      ggsave(file,p,device=e,width=10,height=4.5,dpi=300,bg='white')
    })
  })
  for (kind in c('umap','spatial')) for (ext in c('pdf','png')) local({
    k<-kind; e<-ext
    output[[paste(k,e,sep='_')]] <- downloadHandler(filename=function() paste0('RNHI_',k,'.',e),content=function(file) {
      req(result())
      p<-if(k=='umap') sc_umap_plot(result()) else spatial_plot(result())
      ggsave(file,p,device=e,width=8,height=6,dpi=300,bg='white')
    })
  })
  for (ext in c('pdf','png')) local({
    e<-ext
    output[[paste('spatial_umap',e,sep='_')]] <- downloadHandler(filename=function() paste0('RNHI_spatial_umap.',e),content=function(file) {
      req(result(),result()$mode == 'spatial'); ggsave(file,spot_umap_plot(result()),device=e,width=8,height=6,dpi=300,bg='white')
    })
  })
}
shinyApp(ui,server)
