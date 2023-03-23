# Create scripts necessary for translocation analyses

knitr::purl("code/translocation_survival_analysis.Rmd")
file.rename("translocation_survival_analysis.R", "code/translocation_survival_analysis.R")

knitr::purl("code/translocation_survival_analysis_supp.Rmd")
file.rename("translocation_survival_analysis_supp.R", "code/translocation_survival_analysis_supp.R")

knitr::purl("code/maps.qmd")
file.rename("maps.R", "code/maps.R")