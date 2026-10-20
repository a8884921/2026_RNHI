# Browser smoke test for the three mode-specific example workflows.
# Start the app first on port 3840, then run: Rscript tests/browser_modalities_test.R
library(chromote)

run_mode <- function(mode, selector, expected_text, extra_selector=NULL) {
  b <- ChromoteSession$new(width=1440,height=1000)
  on.exit(b$close(), add=TRUE)
  b$Page$navigate('http://127.0.0.1:3840/')
  b$Page$loadEventFired()
  Sys.sleep(1)
  js_mode <- sprintf("(function(){var x=document.querySelector('#mode input[value=\\\"%s\\\"]'); if(!x) return false; x.click(); x.dispatchEvent(new Event('change',{bubbles:true})); return true;})()", mode)
  stopifnot(isTRUE(b$Runtime$evaluate(js_mode, returnByValue=TRUE)$result$value))
  Sys.sleep(.8)
  b$Runtime$evaluate("document.getElementById('use_demo').click()")
  Sys.sleep(.8)
  b$Runtime$evaluate("document.getElementById('go').click()")
  Sys.sleep(8)
  extra <- if (is.null(extra_selector)) 'true' else sprintf("!!document.querySelector('%s img') && (document.querySelector('%s img')?.naturalWidth||0)>0", extra_selector, extra_selector)
  js <- sprintf("({ready:!!document.getElementById('results_ready'),plot:!!document.querySelector('%s img'),width:document.querySelector('%s img')?.naturalWidth||0,extra:%s,msg:document.getElementById('message')?.innerText||''})", selector, selector, extra)
  state <- b$Runtime$evaluate(js, returnByValue=TRUE)$result$value
  stopifnot(isTRUE(state$ready), isTRUE(state$plot), state$width > 0, isTRUE(state$extra), grepl(expected_text, state$msg, fixed=TRUE))
  cat(mode, 'PASS —', state$msg, '\n')
}

run_mode('bulk', '#dist', 'Scored 3355 records')
run_mode('sc', '#umap', 'Scored 6 records')
run_mode('spatial', '#spatial', 'Scored 144 records', '#spatial_umap')
cat('BROWSER MODALITY EXAMPLE PASS\n')
