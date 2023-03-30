library(data.table)
library(ggplot2)
library(rstan)
library(rethinking)
library(mrmr)
library(lubridate)

## Extract and summarize relevant survival analysis data for MYL Frogs
## Fit a Stan survival model to the analyses
## Extract empirical patterns of recruitment

# The base path of the CMR analysis.  Assumes the working directory is code
# and that CMR analyses are stored two directories up.
cmr_path = file.path("..", "data", "raw", "cmr-analysis")

#####################################
#### Extract survival estimates #####
#####################################

# Extract survival estimates for all translocated lakes
table_files = Sys.glob(file.path(cmr_path, "survival", "*survival_cohort.csv"))

# Loop through each file and load table
yearly_surv_probs = list()
count = 1
for(file in table_files){

	lake_id = strsplit(basename(file), "_")[[1]][1]
	surv_table = fread(file)
	surv_table[, release_date:=as.character(release_date)]

	release_dates = unique(surv_table$release_date)
	for(rl in release_dates){

		# Focus on a specific release date
		ttab = surv_table[release_date == rl & years_since_introduction > 0][order(years_since_introduction)]

		if(nrow(ttab) > 0){
			med_surv = ttab$median_survival
			med_surv = c(1, med_surv)
			yearly_probs = med_surv[2:length(med_surv)] / med_surv[1:(length(med_surv) - 1)]
			new_dat = data.frame(surv_prob=yearly_probs, 
								 lake_id=lake_id, 
								 release_date=rl, 
								 total_surv=med_surv[1:(length(med_surv) - 1)])
			yearly_surv_probs[[count]] = new_dat
			count = count + 1
		}
	}
}


yearly_surv_probs_dt = data.table(do.call(rbind, yearly_surv_probs))[!is.na(surv_prob)]

# Look at yearly survival probabilities only where TOTAL survival probabilities where
# > 0.4.  Did this to remove high variability estimates as cohort population size dwindled
ggplot(yearly_surv_probs_dt[total_surv > .4]) + geom_boxplot(aes(x=lake_id, y=surv_prob))

# Only including the six lakes with "suitable" habitat based on discussion with Roland
keep_lakes = yearly_surv_probs_dt[total_surv > .4][, 
								 .(med_surv=median(surv_prob)), by=.(lake_id)]$lake_id
saveRDS(keep_lakes, file.path("..", "out", "lakes_to_use_in_analysis.rds"))
surv_probs_for_stan = yearly_surv_probs_dt[total_surv > .4 & lake_id %in% keep_lakes]
fwrite(surv_probs_for_stan, file.path("..", "data", "surv_probs_for_stan.csv"))

###################################################
### Stan Model for survival analysis ##############
###################################################

# Now we want to fit a heirarchical beta regression to generate a distribution
# from which I can draw survival values.  The idea is that there is within and between
# site variability in survival.  I would like a distribution from which I can draw
# survival values for any given site as the population proceeds through time. Ignoring
# any temporal correlation in survival probabilities

beta_mod = "
data {

	int<lower=1> N; // Data points
	int<lower=1> S; // Number of unique sites
	int site_ids[N]; // Site random effects
	real<lower=0, upper=1> surv_probs[N]; // Yearly survival probabilities

} parameters {

	real beta0; // Overall intercept
	real alpha_z[S]; // Random effects
	real<lower=0> sigma_alpha; // Among site variability

	real<lower=0> logphi; // Beta distribution dispersion

} transformed parameters{

	real alpha[S];
	real mu[N];

	for(s in 1:S){
		alpha[s] = sigma_alpha*alpha_z[s];
	}

	for(i in 1:N){
		mu[i] = inv_logit(beta0 + alpha[site_ids[i]]);
	}

} model {

	// Natural prior on transformed scale
	alpha_z ~ normal(0, 1);

	// Dispersion parameters
	sigma_alpha ~ normal(0, 2);
	logphi ~ cauchy(0, 2);

	for(i in 1:N){
		surv_probs[i] ~ beta_proportion(mu[i], exp(logphi)); 
	}

}
"

# Build stan data
surv_probs = surv_probs_for_stan$surv_prob 
surv_probs[surv_probs == 1] = 0.999 # No such thing as exactly one in the beta distribution, so set to really close to 1
standata = list(N=nrow(surv_probs_for_stan),
				S=length(unique(surv_probs_for_stan$lake_id)),
				surv_probs=surv_probs,
				site_ids=rethinking::coerce_index(as.factor(surv_probs_for_stan$lake_id)))

# Compile and fit stan model
mod = stan_model(model_code=beta_mod)
fit = sampling(mod, data=standata, cores=4, iter=2000, chains=4)
saveRDS(fit, file.path("..", "out", "adult_survival_probability_model.rds"))
fit

# Reparameterize to define beta distribution in R
mu = 1 / (1 + exp(-extract(fit, "beta0")$beta0)) # Mean of beta
phi = exp(extract(fit, "logphi")$logphi) # disperson parameter of beta

# Shape 1 of beta
a = mu*phi

# Shape 2 of beta
b = phi - a

# Save the parameter estimates for use in the simulation model
saveRDS(list(a=a, b=b), file.path("..", "out", "survival_distribution_parameters.rds"))


############################################
##### Extract patterns of recruitment ######
############################################

# Want to look at patterns of recruitment in Conness pond and translocated lakes
# to qualitatively compare with patterns of recruitment predicted by the model.
keep_lakes = readRDS(file.path("..", "out", "lakes_to_use_in_analysis.rds"))

# Look at Conness data
# con_dat = fread("../data/conness_frog_counts.csv")
# ggplot(con_dat) + geom_line(aes(x=date, y=count, color=lifestage))

# Extract recruitment patterns from cmr fits
all_recruit = list()
for(lake in keep_lakes){

	file = file.path(cmr_path, "model", paste0(lake, "_model.rds"))
	fit = readRDS(file)
	recruitment = data.table(fit$m_fit$summary(variables="B"))
	pop_size = data.table(fit$m_fit$summary(variables="N"))

	# Build the stan data so we know when translocations occurred
	captures = readr::read_csv(file.path(cmr_path, "capture", paste0(lake, "_capture.csv")))
	surveys = readr::read_csv(file.path(cmr_path, "survey", paste0(lake, "_survey.csv")))
	translocations = readr::read_csv(file.path(cmr_path, "translocation", paste0(lake, "_translocation.csv")))
	mrmr_dat = clean_data(captures = captures, surveys = surveys, translocations = translocations)

	# Extract year of primary period
	surveys_dt = data.table(surveys)[, year:=year(survey_date)][, 
									  .(year=unique(year)[1]), by=.(primary_period)][, 
									  setnames(.SD, "primary_period", "time")]

	# Set up introduction times
	introductions = table(mrmr_dat$stan_d$t_intro) 
	introductions = introductions[2:length(introductions)]
	names(introductions) = as.integer(names(introductions)) - 1 # Looks like we need to substract 1
	intro_dt = data.table(time=names(introductions), num_intro=as.vector(unname(introductions)))
	intro_dt$time = as.integer(intro_dt$time)

	dt = data.table(recruitment=recruitment$median,
					pop_size=pop_size$median,
					lake_id=lake,
					time=1:nrow(recruitment))
	dt_merge = merge(dt, intro_dt, key="time", all.x=TRUE)
	dt_merge$num_intro[is.na(dt_merge$num_intro)] = 0
	dt_merge = dt_merge[, corrected_recruitment:=recruitment - num_intro] %>% merge(surveys_dt, key="time", all.x=TRUE)

	all_recruit[[lake]] = dt_merge

}

all_recruit_dt = do.call(rbind, all_recruit)

# Group by year and site
recruit_by_year = all_recruit_dt[, .(total_recruit=sum(corrected_recruitment), pop_size=max(pop_size)), by=.(lake_id, year)]
fwrite(recruit_by_year, file.path("..", "data", "translocation_recruitment_values.csv"))


ggplot(recruit_by_year) + geom_line(aes(x=year, y=total_recruit)) + facet_wrap(~lake_id, scale="free")
ggplot(recruit_by_year) + geom_point(aes(x=log10(pop_size), y=log10(total_recruit + 1), color=lake_id))


