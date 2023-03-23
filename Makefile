# Makefile for translocation analyses
# 
# Instructions
# ------------
# 1. Run `make install_packages` to load the necessary R packages into your environment
# 2. Run `make`

#.PHONY: pngs
#pngs: doc/manuscript/figures/translocation_survival_bysiteid.png doc/manuscript/figures/mcmc_areas_m1d.png doc/manuscript/figures/cond_effects_plot.png doc/manuscript/figures/mcmc_areas_m2b.png doc/manuscript/figures/bdload_beforeafter.png doc/manuscript/figures/map_translocation_points.png

doc/manuscript/figures/translocation_survival_bysiteid.png doc/manuscript/figures/mcmc_areas_m1d.png doc/manuscript/figures/cond_effects_plot.png doc/manuscript/figures/mcmc_areas_m2b.png: code/translocation_survival_analysis.Rmd code/classification_summary.R data/clean/frog_translocation_final.csv
	Rscript code/create_R_scripts.R; \
	Rscript code/translocation_survival_analysis.R; \
	rm -f code/translocation_survival_analysis.R 

# Note: `rm -f code/translocation_survival_analysis.R is not deleted. See GitHub issue for details. 

doc/manuscript/figures/bdload_beforeafter.png: code/translocation_survival_analysis_supp.Rmd data/clean/bd_beforeafter_translocation.csv
	Rscript code/create_R_scripts.R; \ # This script only needs to run once to create all necessary .R files, currently is run in each rule. 
	Rscript code/translocation_survival_analysis_supp.R; \
	rm -f code/translocation_survival_analysis_supp.R 

doc/manuscript/figures/map_translocation_points.png: code/maps.qmd data/maps/translocation_points.csv data/maps/nps-boundaries/yose_boundary.shp data/maps/hillshade-90m/california-color-hillshade-90m.tif data/maps/hillshade-300m/california-color-hillshade-300m.tif
	Rscript code/create_R_scripts.R; \
	Rscript code/maps.R; \
	rm -f code/maps.R


# Run make install_packages to ensure R environment has the necessary packages
install_packages:
	cd code; \
	Rscript install_R_packages.R; \
	cd ..

