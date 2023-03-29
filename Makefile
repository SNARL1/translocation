# Makefile for translocation analysis and population analysis for the manuscript
# 
#
# This makefile has been successfully run on a MacBook pro with Big Sur and a Linux
# machine running Ubuntu.
# 
# 
# Instructions
# ------------
# 1. First, run `make install_packages` to ensure that you have the necessary packages
#    installed in your R environment.  Note that on a Linux machine, you will likely
#	 need to install additional command-line packages to successfully install the
#	 devtools R package.  The error message you get will probably tell you what packages
#	 you need to install. You can also open up R and install each of the packages listed
#	 in your R environment if you are getting errors from install_packages.
# 2. To run all of the analyses, type `make all`
# 2. To run only the translocation analysis execute `make trans_analysis`
# 3. To run only the viability analysis execute `make viability_analysis`
# 4. To clean you directory of extraneous files execute `make clean`

############## RUN ALL ANALYSIS ###################################

.PHONY: all
all: $(OUTPUTS_ALL) $(VIABILITY_OUTPUT)

###################################################################


############## TRANSLOCATION ANALYSIS #############################

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

.PHONY: trans_analysis
trans_analysis: $(OUTPUTS_ALL)

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

# Run `make clean` to remove temporary files that are not removed by the above rules
clean:
	rm -f $(SCRIPT2) $(SCRIPT3)

############## END TRANSLOCATION ANALYSIS #########################


############## START VIABILITY ANALYSIS ############################

VIABILITY_OUTPUT1 = out/lambda_estimates.pdf
VIABILITY_OUTPUT2 = out/compare_surv_probs.csv
VIABILITY_OUTPUT3 = out/compare_surv_probs.jpg
VIABILITY_OUTPUT4 = out/pop_viability_figures_for_manuscript.jpg
VIABILITY_OUTPUT5 = out/pop_viability_figures_for_supp.jpg
VIABILITY_OUTPUT = $(VIABILITY_OUTPUT1) $(VIABILITY_OUTPUT2) $(VIABILITY_OUTPUT3) $(VIABILITY_OUTPUT4) $(VIABILITY_OUTPUT5)
VIABILITY_SCRIPT = code/extinction_and_growth_rate_analysis.R
POP_FXNS = code/population_viability_functions.R

# CMR output
CMR_SURVEYS = $(wildcard data/raw/cmr-analysis/survey/*_survey.csv) 
CMR_MODELS = $(wildcard data/raw/cmr-analysis/model/*_model.rds)
CMR_TABLES = $(wildcard data/raw/cmr-analysis/survival/*survival_cohort.csv)
CMR_CAPTURE = $(wildcard data/raw/cmr-analysis/capture/*_capture.csv) 
CMR_TRANS = $(wildcard data/raw/cmr-analysis/translocation/*_translocation.csv) 

# Survival script output
SURV_OUTPUT1 = out/yearly_survival_estimates_all_individuals.csv
SURV_OUTPUT2 = out/yearly_survival_estimates_only_recruited_individuals.csv
SURV_OUTPUT = $(SURV_OUTPUT1) $(SURV_OUTPUT2)
SURV_SCRIPT = code/extract_survival_estimates.R
SURV_FXNS = code/survival_table_nontranslocated.R

# Translocation script output
TRANS_OUTPUT1 = data/translocation_recruitment_values.csv
TRANS_OUTPUT2 = out/adult_survival_probability_model.rds
TRANS_OUTPUT3 = out/lakes_to_use_in_analysis.rds
TRANS_OUTPUT = $(TRANS_OUTPUT1) $(TRANS_OUTPUT2) $(TRANS_OUTPUT3)
TRANS_SCRIPT = code/extract_translocation_data.R

# Abundance script output. 70449 is just a placeholder to stimulate building
ABUND_OUTPUT1 = data/clean/abundance_and_recruitment/70449_abundance.csv
ABUND_OUTPUT2 = data/clean/abundance_and_recruitment/70449_recruitment.csv
ABUND_OUTPUT = $(ABUND_OUTPUT1) $(ABUND_OUTPUT2)
ABUND_SCRIPT = code/extract_abundance_trajectories.R

.PHONY: viability_analysis
viability_analysis : $(VIABILITY_OUTPUT)

# Population viability analysis
$(VIABILITY_OUTPUT) : $(TRANS_OUTPUT) $(TRANS_SCRIPT) $(SURV_OUTPUT) $(SURV_SCRIPT) $(ABUND_OUTPUT) $(ABUND_SCRIPT) $(VIABILITY_SCRIPT) $(POP_FXNS)
	cd code; \
	Rscript $(notdir $(VIABILITY_SCRIPT)); \
	cd ..

# Translocation output
$(TRANS_OUTPUT) : $(CMR_TABLES) $(CMR_MODELS) $(CMR_CAPTURE) $(CMR_SURVEYS) $(CMR_TRANS) $(TRANS_SCRIPT)
	cd code; \
	Rscript $(notdir $(TRANS_SCRIPT)); \
	cd ..

# Abundance trajectories
$(ABUND_OUTPUT) : $(CMR_SURVEYS) $(CMR_MODELS) $(ABUND_SCRIPT) $(POP_FXNS)
	cd code; \
	Rscript $(notdir $(ABUND_SCRIPT)); \
	cd ..

# Survival estimates
$(SURV_OUTPUT) : $(CMR_SURVEYS) $(CMR_MODELS) $(SURV_SCRIPT) $(SURV_FXNS)
	cd code; \
	Rscript $(notdir $(SURV_SCRIPT)); \
	cd ..

################## END VIABILITY ANALYSIS #################################

.PHONY: install_packages

# Run make install_packages to ensure R environment has the necessary packages
install_packages:
	cd code; \
	Rscript install_R_packages.R; \
	cd ..

