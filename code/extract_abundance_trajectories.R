## Extract abundance and recruitment trajectories for different lakes

library(cmdstanr)
library(data.table)
library(magrittr)
source("population_viability_functions.R")

cmr_path = file.path("..", "data", "raw", "cmr-analysis")

# Lakes to include
lake_ids = c(70414, 70370, 70134, 70505, 70550, 70641, 70619, 70449, 70413,
			 70628, 70556, 74976)
thin_values = rep(10, length(lake_ids))

for(l in 1:length(lake_ids)){

	lake_number = lake_ids[l]
	cat("Working on lake", lake_number, "\n")
	modelpath = file.path(cmr_path, "model", paste0(lake_number, "_model.rds"))
	datapath = file.path(cmr_path, "survey", paste0(lake_number, "_survey.csv"))

	model = readRDS(modelpath)
	abundance_table = abundances_from_model(model, "abundance")
	recruitment_table = abundances_from_model(model, "recruitment")

	# Save tables to disk
	dpath = file.path("..", "data", "clean", "abundance_and_recruitment")
	if(!dir.exists(dpath)){
		dir.create(dpath)
	}
	fwrite(abundance_table, file.path("..", "data", "clean", "abundance_and_recruitment", paste0(lake_number, "_abundance.csv")))
	fwrite(recruitment_table, file.path("..", "data", "clean", "abundance_and_recruitment", paste0(lake_number, "_recruitment.csv")))
}






