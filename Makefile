# Makefile for translocation analyses
# 
# Instructions
# ------------
# 1. Run `make install_packages` to load the necessary R packages into your environment
# 2. Run `make`

# .PHONY: pngs
# pngs: translocation_survival_bysiteid.png mcmc_areas_m1d.png cond_effects_plot.png mcmc_areas_m2b.png

doc/manuscript/figures/translocation_survival_bysiteid.png doc/manuscript/figures/mcmc_areas_m1d.png doc/manuscript/figures/cond_effects_plot.png doc/manuscript/figures/mcmc_areas_m2b.png: code/translocation_survival_analysis.Rmd code/classification_summary.R data/clean/frog_translocation_final.csv
	Rscript code/create_R_scripts.R; \
	Rscript code/translocation_survival_analysis.R; \
	rm -f code/translocation_survival_analysis.R 

# Note: `rm -f code/translocation_survival_analysis.R is not deleted. See GitHub issue for details. 


# Run make install_packages to ensure R environment has the necessary packages
install_packages:
	cd code; \
	Rscript install_R_packages.R; \
	cd ..

