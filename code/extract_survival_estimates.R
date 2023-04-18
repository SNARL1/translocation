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


cmr_path = file.path("..", "data", "raw", "cmr-analysis")

# Lakes to include
lake_ids = c(70414, 70370, 70134, 70505, 70550, 70641, 70619, 70449, 70413,
			 70628, 70556, 74976)

thin_values = rep(10, length(lake_ids))

all_surv_ests = list()
only_recruited_surv_ests = list()
only_translocated_surv_ests = list()

count = 1
for(l in 1:length(lake_ids)){

	lake_number = lake_ids[l]
	cat("Working on lake", lake_number, "\n")
	modelpath = file.path(cmr_path, "model", paste0(lake_number, "_model.rds"))
	datapath = file.path(cmr_path, "survey", paste0(lake_number, "_survey.csv"))

	model = readRDS(modelpath)
	survey_dat = fread(datapath)

	### Extract survival of only non-translocated individuals ###

	# Calling function survival_table from survival_table_nontranslocated.R
	res = survival_table(model, thin=thin_values[l], lowerupper=c(0.025, 0.975),
						 translocated=FALSE)

	# Get averaged survival probabilities
	surv_list = build_surv_list(res, drop_ones=TRUE)

	# Store survival probabilities
	only_recruited_surv_ests[[count]] = data.frame(surv_prob_med=surv_list[['med']],
												   surv_prob_lower=surv_list[['lower']],
												   surv_prob_upper=surv_list[['upper']],
											   	   lake_id=lake_number)

	### Extract survival probabilities of translocated individuals ###

	# Calling function survival_table from survival_table_nontranslocated.R

	res_trans = survival_table(model, thin=thin_values[l], lowerupper=c(0.025, 0.975),
						 	   translocated=TRUE)

	# When calculating the mean death rate, only use values that were estimated
	# from more than 10 individuals.  Otherwise, we can over weight survival
	# estimates that come from just one or two individuals. 
	res_trans$surv_probs = res_trans$surv_probs[num_obs > 10]
	surv_list_trans = build_surv_list(res_trans, drop_ones=FALSE)

	# Store survival probabilities
	only_translocated_surv_ests[[count]] = data.frame(surv_prob_med=surv_list_trans[['med']],
												   	  surv_prob_lower=surv_list_trans[['lower']],
												      surv_prob_upper=surv_list_trans[['upper']],
											   	      lake_id=lake_number)

	### Extract survival probabilities of non-translocated and translocated individuals ###

	# Note that this analysis explicitly extracts the survival transitions probabilities from the MRMR model,
	# which are time-varying.  The transition probabilities are accounting for transitions
	# of translocated and non-translocated individuals.

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
	lower = apply(phi, 2, quantile, 0.025)
	upper = apply(phi, 2, quantile, 0.975)
	med = apply(phi, 2, quantile, 0.5)
	unq_primary$surv_prob = med[2:length(med)]

	# Put the probabilities in terms of rates per day
	rates = t(t(-log(phi[, 3:length(med)])) / unq_primary$deltat_days[2:nrow(unq_primary)])
	med_rates = apply(rates, 2, quantile, 0.5)
	lower_rates = apply(rates, 2, quantile, 0.025)
	upper_rates = apply(rates, 2, quantile, 0.975)

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

# Save results
all_surv_ests_dt = data.table(do.call(rbind, all_surv_ests))
fwrite(all_surv_ests_dt, file.path("..", "out", "yearly_survival_estimates_all_individuals.csv"))

only_recruited_surv_ests_dt = data.table(do.call(rbind, only_recruited_surv_ests))
fwrite(only_recruited_surv_ests_dt, file.path("..", "out", "yearly_survival_estimates_only_recruited_individuals.csv"))

only_translocated_surv_ests_dt = data.table(do.call(rbind, only_translocated_surv_ests))
fwrite(only_recruited_surv_ests_dt, file.path("..", "out", "yearly_survival_estimates_only_translocated_individuals.csv"))

# Make a plot comparing the survival estimates for three lakes with naturally recruited adults
# and translocated adults

# These three lakes have both naturally recruited adults and translocated adults
three_lakes = c(70550, 70413, 70449)

or_dt = only_recruited_surv_ests_dt[lake_id %in% three_lakes]
or_dt$type = "Naturally recruited"

t_dt = only_translocated_surv_ests_dt[lake_id %in% three_lakes]
t_dt$type = "Translocated" 

c_dt = all_surv_ests_dt[lake_id %in% three_lakes][, .(surv_med, surv_lower, surv_upper, lake_id)]
colnames(c_dt) = c("surv_prob_med", "surv_prob_lower", "surv_prob_upper", "lake_id")
c_dt$type = "Combined"

all_dt = rbind(or_dt, t_dt)[lake_id %in% three_lakes]
all_dt$lake_id = as.factor(all_dt$lake_id)

ggplot(all_dt) + 
				geom_point(aes(x=lake_id, y=surv_prob_med, color=type), position=position_dodge(width=0.55)) +
				geom_errorbar(aes(x=lake_id, 
													ymin=surv_prob_lower, 
													ymax=surv_prob_upper, color=type), width=0.25, position=position_dodge(width=0.55)) +
				scale_color_manual(values=c("#1b9e77", "#d95f02", "black")) +
				ylab("Yearly survival probability") + xlab("Lake ID") + 
				theme_classic() + theme(legend.title=element_blank()) + ylim(c(0, 1))

ggsave(file.path("..", "out", "compare_surv_probs.jpg"), width=6, height=4)





