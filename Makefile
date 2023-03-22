# Makefile for translocation analyses
# 
# Instructions
# ------------
# 1. Run `make install_packages` to load the necessary R packages into your environment
# 2. Run `make create_rscripts` to create the necessary analysis R scripts
# 3. Run `make`

# .PHONY: pngs
# pngs: translocation_survival_bysiteid.png mcmc_areas_m1d.png cond_effects_plot.png mcmc_areas_m2b.png

translocation_survival_bysiteid.png mcmc_areas_m1d.png cond_effects_plot.png mcmc_areas_m2b.png &: code/translocation_survival_analysis.R code/classification_summary.R data/clean/frog_translocation_final.csv
	Rscript code/translocation_survival_analysis.R

.PHONY: clean
clean:
	rm -f code/translocation_survival_analysis.R 

# Run make install_packages to ensure R environment has the necessary packages
install_packages:
	cd code; \
	Rscript install_R_packages.R; \
	cd ..

# Run `make create_rscripts` to ensure R environment has the necessary .R files
create_rscripts:
	Rscript code/create_R_scripts.R