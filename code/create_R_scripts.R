# Create scripts necessary for translocation analyses

knitr::purl("code/translocation_survival_analysis.Rmd")
file.rename("translocation_survival_analysis.R", "code/translocation_survival_analysis.R")