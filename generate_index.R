
setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples")

library(GetoptLong)
files = list.files(pattern = "html$")
files = setdiff(files, c("index.html", "readme.html", "README.html"))

qqcat("- [@{files}](@{files})\n", file = "readme.md")

library(knitr)
knit2html("readme.md")

servr::httd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples")
