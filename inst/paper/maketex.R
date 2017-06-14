
#render to markdown
rmarkdown::render("~/src/rnabioco/valr/inst/paper/figs.Rmd",
                  output_format = "html_document",
                  clean = F)

# render to latex
rmarkdown::pandoc_convert("~/src/rnabioco/valr/inst/paper/figs.knit.md", to = "latex", output = "figs.tex")
