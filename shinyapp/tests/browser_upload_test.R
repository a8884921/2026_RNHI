# End-to-end browser check. Start the app first:
#   Rscript -e "shiny::runApp('.', host='127.0.0.1', port=3840, launch.browser=FALSE)"
library(chromote)

file.copy('/Volumes/mac2/Project2/PP_AML空间免疫2/09_shinyapp网站/composition_A4.csv', '/tmp/composition_A4.csv', overwrite=TRUE)

b <- ChromoteSession$new(width=1440, height=1200)
on.exit(b$close(), add=TRUE)
b$Page$navigate('http://127.0.0.1:3840/')
b$Page$loadEventFired()
Sys.sleep(1)
doc <- b$DOM$getDocument()
file_node <- b$DOM$querySelector(nodeId=doc$root$nodeId, selector='#file')
stopifnot(!is.null(file_node$nodeId))
b$DOM$setFileInputFiles(files=list('/tmp/composition_A4.csv'), nodeId=file_node$nodeId)
Sys.sleep(1)
b$Runtime$evaluate('document.getElementById("go").click()')
Sys.sleep(10)
state <- b$Runtime$evaluate(
  '({ready:!!document.getElementById("results_ready"),msg:document.getElementById("message")?.innerText||"",rows:document.querySelectorAll("#table tbody tr").length,imgs:Array.from(document.querySelectorAll("#dist img, #heat img")).map(x=>x.naturalWidth)})',
  returnByValue=TRUE)$result$value
stopifnot(isTRUE(state$ready), grepl('Scored 3355 records; 3310 query samples; 45 healthy reference donors', state$msg, fixed=TRUE))
stopifnot(state$rows > 0, length(state$imgs) == 2, all(state$imgs > 0))
b$screenshot('tests/browser_real_data_results.png', selector='body', delay=1)
cat('BROWSER E2E PASS\n', state$msg, '\n', 'Visible table rows:', state$rows, '\n', 'Plot widths:', paste(state$imgs, collapse=', '), '\n')
