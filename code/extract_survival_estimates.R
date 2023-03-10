## Script to extract all individuals using the mrmr output
## Here we extract the time varying survival estimates,
## calculate the time period over which these phi estimates were
## estimated, and then calculate a weighted mean of the survival rate.
## This survival estimate includes the survival probability for all individuals
## including translocated individuals and survival estiamtes for non-translocated 
## individuals

library(cmdstanr)
library(data.table)
library(magrittr)
source("survival_table_nontranslocated.R")

# Lakes to include
lake_ids = c(70414, 70370, 70134, 70505, 70550, 70641, 70619, 70449, 70413,
			 70628, 70556, 74976)

thin_values = rep(10, length(lake_ids))

all_surv_ests = list()
only_recruited_surv_ests = list()

count = 1
for(l in 1:length(lake_ids)){

	lake_number = lake_ids[l]
	cat("Working on lake", lake_number, "\n")
	modelpath = file.path("..", "cmr-analysis", "out", "model", paste0(lake_number, "_model.rds"))
	datapath = file.path("..", "cmr-analysis", "data", "clean", paste0(lake_number, "_survey.csv"))

	model = readRDS(modelpath)
	survey_dat = fread(datapath)

	### Extract survival of only non-translocated individuals ###
	res = survival_table(model, thin=thin_values[l])

	# Only execute if there are non-translocated pittag animals
	if(res$num_pittags > 0) { 

		surv_probs = res$surv_probs

		# Find the first stretch of primary periods where survival is 1 and remove them
		# Doing this because individuals were not actually observed during this period
		# and it is always conditional on recruitment.  Thus, survival has to be 1 
		# as we know they were alive later.

		nm_append = c("lower", "med", "upper")
		rle_res = rle(surv_probs[["med"]] == 1)
		cut = ifelse(rle_res$values[1] == TRUE, rle_res$lengths[1], 0)
		surv_probs = surv_probs[(cut + 1):nrow(surv_probs)]
		surv_list = list()
		for(nm in nm_append){

			# Get the weighted survival probabilities
			total_days = sum(surv_probs$days)
			weights = (surv_probs$days / total_days)
			mean_rate = sum(weights * surv_probs[[paste0('death_rate_day_', nm)]])

			yearly_surv_prob = exp(-mean_rate*365)
			surv_list[[nm]] = yearly_surv_prob
		}
	} else{
		surv_list = list(lower=NA, med=NA, upper=NA)
	}

	# Store survival probabilities
	only_recruited_surv_ests[[count]] = data.frame(surv_prob_med=surv_list[['med']],
												   surv_prob_lower=surv_list[['lower']],
												   surv_prob_upper=surv_list[['upper']],
											   	   lake_id=lake_number)

	### Extract survival probabilities of non-translocated and translocated individuals ###

	introduced = model$data$stan_d$introduced
	T = model$data$stan_d$T # Primary periods
	Y = model$data$stan_d$Y

	beta_phi = as.data.table(model$m_fit$draws("beta_phi", format = "draws_df"))[, c(".draw", ".chain", ".iteration"):=NULL]
	eps_phi = as.data.table(model$m_fit$draws("eps_phi", format = "draws_df"))[, c(".draw", ".chain", ".iteration"):=NULL]
	sigma_phi = as.data.table(model$m_fit$draws("sigma_phi", format = "draws_df"))[, c(".draw", ".chain", ".iteration"):=NULL]

	# Generate time-varying survival probabilities
	phi = array(NA, dim=c(nrow(beta_phi), T))
	for(t in 1:T){
		phi[, t] = 1 / (1 + exp(-(beta_phi[[1]] + eps_phi[[t]])))
	}

	# Look at the unique primary periods
	unq_primary = survey_dat[!duplicated(primary_period)]
	dates = as.POSIXlt(unq_primary$survey_date)
	unq_primary$deltat_days = c(0, diff(dates))

	# Get the overall survival probabilities
	lower = apply(phi, 2, quantile, 0.25)
	upper = apply(phi, 2, quantile, 0.75)
	med = apply(phi, 2, quantile, 0.5)
	unq_primary$surv_prob = med[2:length(med)]

	# Put the probabilities in terms of rates per day
	rates = t(t(-log(phi[, 3:length(med)])) / unq_primary$deltat_days[2:nrow(unq_primary)])
	med_rates = apply(rates, 2, quantile, 0.5)
	lower_rates = apply(rates, 2, quantile, 0.25)
	upper_rates = apply(rates, 2, quantile, 0.75)

	# Get the average death rate per day over the time scale of interest
	all_rates = list(lower=lower_rates, med=med_rates, upper=upper_rates)
	probs = list()
	for(j in 1:3){

		ar = all_rates[[j]]
		unq_primary$death_rate = c(0, ar)
		total_days = sum(unq_primary$deltat_days)
		weights = (unq_primary$deltat_days / total_days)
		mean_rate = sum(weights * unq_primary$death_rate)
		probs[[names(all_rates)[j]]] = exp(-mean_rate*365)

	}

	# Yes, the upper and lower are purposly mis-matched because
	# lower death rate is upper survival probability
	all_surv_ests[[count]] = data.frame(surv_med=probs[['med']],
										surv_upper=probs[['lower']],
										surv_lower=probs[['upper']],
										lake_id=lake_number)
	count = count + 1
}

all_surv_ests_dt = data.table(do.call(rbind, all_surv_ests))
fwrite(all_surv_ests_dt, file.path("..", "out", "yearly_survival_estimates_all_individuals.csv"))

only_recruited_surv_ests_dt = data.table(do.call(rbind, only_recruited_surv_ests))
fwrite(only_recruited_surv_ests_dt, file.path("..", "out", "yearly_survival_estimates_only_recruited_individuals.csv"))



