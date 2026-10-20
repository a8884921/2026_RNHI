# Browser E2E check for direct Seurat/Visium metadata uploads.
# Start the app on port 3840 first.
library(chromote)

lineages <- c('MyeProg','Mono','Gran','Mega','EarlyEry','LateEry','B','NK','DC')
sc_meta <- data.frame(orig.ident=rep(c('Van_meta_01','Van_meta_02'), each=length(lineages)),
                      celltype=rep(lineages, 2),
                      umap_1=rep(seq_along(lineages), 2), umap_2=rep(rev(seq_along(lineages)), 2))
sp_meta <- data.frame(spot_id=rep(sprintf('spot_%02d',1:4), each=length(lineages)),
                      celltype=rep(lineages, 4),
                      x=rep(rep(1:2, each=2), each=length(lineages)),
                      y=rep(rep(1:2, 2), each=length(lineages)))
write.csv(sc_meta,'/tmp/rnhi_sc_meta.csv',row.names=FALSE,quote=FALSE)
write.csv(sp_meta,'/tmp/rnhi_sp_meta.csv',row.names=FALSE,quote=FALSE)

run_upload <- function(mode, path, plot_selector, expected) {
  b <- ChromoteSession$new(width=1440,height=1000)
  on.exit(b$close(), add=TRUE)
  b$Page$navigate('http://127.0.0.1:3840/')
  b$Page$loadEventFired(); Sys.sleep(1)
  js_mode <- sprintf("var x=document.querySelector('#mode input[value=\\\"%s\\\"]'); x.click(); x.dispatchEvent(new Event('change',{bubbles:true}));", mode)
  b$Runtime$evaluate(js_mode); Sys.sleep(.7)
  doc <- b$DOM$getDocument(); node <- b$DOM$querySelector(nodeId=doc$root$nodeId,selector='#file')
  b$DOM$setFileInputFiles(files=list(path),nodeId=node$nodeId); Sys.sleep(1)
  b$Runtime$evaluate('document.getElementById("go").click()'); Sys.sleep(8)
  js <- sprintf("({ready:!!document.getElementById('results_ready'),msg:document.getElementById('message')?.innerText||'',plot:document.querySelector('%s img')?.naturalWidth||0,rows:document.querySelectorAll('#table tbody tr').length})",plot_selector)
  state <- b$Runtime$evaluate(js,returnByValue=TRUE)$result$value
  stopifnot(isTRUE(state$ready), state$plot > 0, state$rows > 0, grepl(expected,state$msg,fixed=TRUE))
  cat(mode,'metadata upload PASS —',state$msg,'\n')
}
run_upload('sc','/tmp/rnhi_sc_meta.csv','#umap','Scored 2 records')
run_upload('spatial','/tmp/rnhi_sp_meta.csv','#spatial','Scored 4 records')
cat('BROWSER METADATA INPUT PASS\n')
