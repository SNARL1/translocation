# Makefile for translocation analyses
# 
# Instructions
# ------------
# 1. Run `make install_packages` to load the necessary R packages into your environment
# 2. Run `make`
# 3. Run `make clean` to remove temporary .R scripts not removed by other rules

# Outputs
OUTPUT1 = doc/manuscript/figures/translocation_survival_bysiteid.png
OUTPUT2 = doc/manuscript/figures/mcmc_areas_m1d.png
OUTPUT3 = doc/manuscript/figures/cond_effects_plot.png
OUTPUT4 = doc/manuscript/figures/mcmc_areas_m2b.png
OUTPUT5 = doc/manuscript/figures/bdload_beforeafter.png
OUTPUT6 = doc/manuscript/figures/map_translocation_points.png
OUTPUTS_ALL = $(OUTPUT1) $(OUTPUT2) $(OUTPUT3) $(OUTPUT4) $(OUTPUT5) $(OUTPUT6)
OUTPUTS_MAIN = $(OUTPUT1) $(OUTPUT2) $(OUTPUT3) $(OUTPUT4)

# Dependencies
DEPEND1 = code/translocation_survival_analysis.Rmd
DEPEND2 = code/classification_summary.R
DEPEND3 = data/clean/frog_translocation_final.csv
DEPEND4 = code/translocation_survival_analysis_supp.Rmd
DEPEND5 = data/clean/bd_beforeafter_translocation.csv
DEPEND6 = code/maps.qmd
DEPEND7 = data/maps/translocation_points.csv
DEPEND8 = data/maps/nps-boundaries/yose_boundary.shp
DEPEND9 = data/maps/hillshade-90m/california-color-hillshade-90m.tif
DEPEND10 = data/maps/hillshade-300m/california-color-hillshade-300m.tif

# .R scripts
SCRIPT1 = code/create_R_scripts.R
SCRIPT2 = code/translocation_survival_analysis.R
SCRIPT3 = code/translocation_survival_analysis_supp.R
SCRIPT4 = code/maps.R

.PHONY: pngs
pngs: $(OUTPUTS_ALL)

$(OUTPUTS_MAIN) &: $(DEPEND1) $(DEPEND2) $(DEPEND3)
	Rscript $(SCRIPT1); \
	Rscript $(SCRIPT2); \
	rm -f $(SCRIPT2)

# Used `&:` to separate targets and dependencies to prevent `translocation_survival_analysis.R` from running for each target (instead of once for all targets).
# Removal of SCRIPT2 is not occurring, perhaps due to problem in terminating parallel processing (described in GitHub issue). 

$(OUTPUT5): $(DEPEND4) $(DEPEND5)
	Rscript $(SCRIPT1); \
	Rscript $(SCRIPT3); \
	rm -f $(SCRIPT3)
	
# The script `create_R_scripts.R` only needs to run once to create all .R files for this Makefile, currently is run in each rule. Can this be improved?
# Removal of SCRIPT3 is not occurring, perhaps due to problem in terminating parallel processing (described in GitHub issue). 


$(OUTPUT6): $(DEPEND6) $(DEPEND7) $(DEPEND8) $(DEPEND9) $(DEPEND10)
	Rscript $(SCRIPT1); \
	Rscript $(SCRIPT4); \
	rm -f $(SCRIPT4)

# SCRIPT4 is removed as expected.

# Run `make install_packages` to ensure R environment has the necessary packages
install_packages:
	cd code; \
	Rscript install_R_packages.R; \
	cd ..

# Run `make clean` to remove temporary files that are not removed by the above rules
clean:
	rm -f $(SCRIPT2) $(SCRIPT3)
