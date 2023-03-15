# Makefile for population viability analysis
# 
# 
# Instructions
# ------------
# 1. First run `make install_packages` to load the necessary R packages into your environment
# 2. Then run `make`
# 
# Note that the make file assumes that the `cmr-analysis` repository is cloned and is living in the
# same directory as the `translocation` repository


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

.PHONY: install_packages

# Run make install_packages to ensure R environment has the necessary packages
install_packages:
	cd code; \
	Rscript install_R_packages.R; \
	cd ..

