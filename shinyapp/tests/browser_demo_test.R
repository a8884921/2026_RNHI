# End-to-end check for the no-upload demo button. Start the app on port 3840.
library(chromote)
b <- ChromoteSession$new(width=1000,height=800)
on.exit(b$close(),add=TRUE)
b$Page$navigate('http://127.0.0.1:3840/')
b$Page$loadEventFired(); Sys.sleep(1)
header <- b$Runtime$evaluate('({subtitle:document.querySelector(".brand-subtitle")?.innerText||"",logo:document.querySelector(".brand-logo")?.getBoundingClientRect().width||0,tables:document.querySelectorAll(".example-table").length,pre:document.querySelectorAll("pre").length})',returnByValue=TRUE)$result$value
stopifnot(identical(header$subtitle,'Residual Normal Hematopoiesis Index'), header$logo >= 200, header$tables >= 2, header$pre == 0)
bulk_options <- b$Runtime$evaluate('Object.keys(document.getElementById("bulk_demo")?.selectize?.options||{})',returnByValue=TRUE)$result$value
stopifnot(all(c('composition','beataml','expression') %in% bulk_options))
b$Runtime$evaluate('document.getElementById("use_demo").click()')
Sys.sleep(1)
b$Runtime$evaluate('document.getElementById("go").click()')
Sys.sleep(10)
message <- b$Runtime$evaluate('document.getElementById("message")?.innerText||""')$result$value
stopifnot(grepl('Scored 3355 records; 3310 query samples; 45 healthy reference donors',message,fixed=TRUE))
cat('DEMO BUTTON PASS\n',message,'\n')

# The expression selector must execute the bundled locked MuSiC/fallback path.
b$Runtime$evaluate('document.getElementById("refresh_app").click()')
Sys.sleep(.5)
b$Runtime$evaluate('document.getElementById("bulk_demo-selectized").click()')
Sys.sleep(.2)
b$Runtime$evaluate('document.querySelector("#bulk_demo + div .option[data-value=expression]")?.click()')
Sys.sleep(.5)
b$Runtime$evaluate('document.getElementById("use_demo").click()')
b$Runtime$evaluate('document.getElementById("go").click()')
Sys.sleep(15)
music_message <- b$Runtime$evaluate('document.getElementById("message")?.innerText||""')$result$value
music_alert <- b$Runtime$evaluate('document.querySelector(".alert-danger")?.innerText||""')$result$value
stopifnot(grepl('Scored 8 records',music_message,fixed=TRUE), !nzchar(music_alert),
          isTRUE(b$Runtime$evaluate('!!document.querySelector("#dist img") && (document.querySelector("#dist img")?.naturalWidth||0)>0',returnByValue=TRUE)$result$value))
cat('BULK EXPRESSION EXAMPLE PASS\n',music_message,'\n')

# BeatAML2 composition is a separate selectable bulk example.
b$Runtime$evaluate('document.getElementById("refresh_app").click()')
Sys.sleep(.5)
b$Runtime$evaluate('document.getElementById("bulk_demo-selectized").click()')
Sys.sleep(.2)
b$Runtime$evaluate('document.querySelector("#bulk_demo + div .option[data-value=beataml]")?.click()')
Sys.sleep(.5)
b$Runtime$evaluate('document.getElementById("use_demo").click()')
b$Runtime$evaluate('document.getElementById("go").click()')
Sys.sleep(8)
beataml_message <- b$Runtime$evaluate('document.getElementById("message")?.innerText||""')$result$value
beataml_alert <- b$Runtime$evaluate('document.querySelector(".alert-danger")?.innerText||""')$result$value
stopifnot(grepl('Scored 8 records',beataml_message,fixed=TRUE), grepl('beatAML2_demo.csv',beataml_message,fixed=TRUE), !nzchar(beataml_alert))
cat('BEATAML2 COMPOSITION EXAMPLE PASS\n',beataml_message,'\n')
