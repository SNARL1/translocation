library(tidyverse)
library(data.table)

assign_parameters = function(adult_survR, adult_survT, new_params=list()){
	# Assign the model parameters based on a value for yearly
	# adult survival for naturally recruited adults and translocated adults
	#
	# Parameters
	# ----------
	# adult_survR : float
	# 	Yearly probability of adult survival for recruits, between 0 and 1
	# 
	# adult_survT : float
	#	  Yearly probability of adult survival for translocated adults, between 0 and 1
	#
	# new_params : list
	#		A list of any new parameters that you want to update
	#
	# Returns
	# -------
	# : list
	#	List of parameters

	params = list(sigma_L1=0.7, 
			  sigma_L2=0.7, 
			  sigma_L3=0.7, 
			  sigma_J1=0.25, 
			  sigma_J2=0.5,
			  sigma_AR=adult_survR,
			  sigma_AT=adult_survT,
			  p_L1=1, p_L2=0.25, p_J1=0.25,
			  p_F=0.5, F=100, 
			  reproduction="poisson")

	# Set any new parameters
	for(p in names(new_params)){
		params[[p]] = new_params[[p]]
	}

	return(params)
}

build_transition_matrix = function(p){
	# Build the density-independent transition matrix
	#
	# Parameters
	# ----------
	# p : list
	# 	List of parameters needed in the transition matrix
	#
	# Returns
	# -------
	# : the transition matrix

	# Columns of matrix
	L1 = c(0, p$sigma_L1*p$p_L1, 0, p$sigma_L1*(1 - p$p_L1), 0, 0, 0)
	L2 = c(0, 0, p$sigma_L2*p$p_L2, p$sigma_L2*(1 - p$p_L2), 0, 0, 0)
	L3 = c(0, 0, 0, p$sigma_L3, 0, 0, 0)
	J1 = c(0, 0, 0, 0, p$sigma_J1*p$p_J1, p$sigma_J1*(1 - p$p_J1), 0)
	J2 = c(0, 0, 0, 0, 0, p$sigma_J2, 0)
	AR = c(0, 0, 0, 0, 0, p$sigma_AR, 0)
	AT = c(0, 0, 0, 0, 0, 0, p$sigma_AT)

	T = cbind(L1, L2, L3, J1, J2, AR, AT)
	return(T)

}

transition_matrix_w_death = function(p){
	# Augment the transition matrix with a death class

	T = build_transition_matrix(p)

	# Add death
	death_probs = 1 - colSums(T)
	Tnew = rbind(T, death_probs)
	Tnew = cbind(Tnew, c(rep(0, 7), 1))
	return(Tnew)

}

environmental_simulation = function(steps, initial_values, muR, muT, phi, new_params=list()){
	# Iterate deterministic projections with environmentally
	# variable transition matrices arbitrarily far foward in time
	# There is no demographic stochasticity in this simulation.
	# 
	# Parameters
	# ----------
	# steps : int
	# 	yearly time steps to simulate over
	# initial_values : vector
	#		Vector of length 6 specifying 
	# phi : float
	#		Dispersion parameter of a beta distribution
	# muR : float
	#		Average recruited adult survival probability
	# muT : float
	#		Average translocated adult survival probability
	#
	# Returns
	# -------
	# : Projection of initial values arbitrarily far forward in time

	aR = muR*phi
	aT = muT*phi

	bR = phi - aR
	bT = phi - aT

	for(i in 1:(steps)){

		adult_survR = rbeta(1, aR, bR)
		adult_survT = rbeta(1, aT, bT)

		params = assign_parameters(adult_survR, adult_survT, new_params=new_params)

		T = build_transition_matrix(params)
		F = build_fecundity_matrix(params)
		A = T + F

		if(i == 1){
			Anext = A
		} else{
			Anext = A %*% Anext
		}
	}

	return(Anext %*% initial_values)

}

stochastic_step = function(state_vars, p){
	# Iterate the matrix forward one time step using a stochastic update
	#
	# Parameters
	# ----------
	# state_vars : array
	#		Values of the state variables at the current time step
	# p : list
	# 	List of parameter values
	#
	# Returns
	# -------
	# : array
	#	Updated state variables after stochastic step

	# Stochastic transitions
	Tnew = transition_matrix_w_death(p)
	D = sapply(1:length(state_vars), function(i) { 
										if(state_vars[i] > 0){
											if(state_vars[i] < 10000){
												x = rmultinom(1, state_vars[i], Tnew[, i]) 
											}
											else{
												# Deterministic calculation
												x = round(state_vars[i]*Tnew[, i])
											}
										}
										else {
											x = rep(0, length(state_vars) + 1)
										}
										return(x)})
	new_states = rowSums(D)[1:length(state_vars)]
	recruited = D[6, 5] + D[6, 4] # Number recruited from Juvenile 1, 2 to adult R

	# How many are reproducing?
	n = length(state_vars)

	# Account for overflow
	if(state_vars[(n - 1)] < 10000){
		reproducing_adultsR = rbinom(1, state_vars[(n - 1)], p$sigma_AR*p$p_F)
		reproducing_adultsT = rbinom(1, state_vars[n], p$sigma_AT*p$p_F)

		# Only females reproduce so divide by 2
		if(p$reproduction == "nbd"){
			new_larvaeR = rnbinom(1, mu=(reproducing_adultsR / 2)*p$F, size=1)
			new_larvaeT = rnbinom(1, mu=(reproducing_adultsT / 2)*p$F, size=1)
		} else{
			new_larvaeR = rpois(1, (reproducing_adultsR) / 2 * p$F)
			new_larvaeT = rpois(1, (reproducing_adultsT) / 2 * p$F)
		}
	} else{

		# Just calculate deterministically if population is large...
		new_larvaeR = state_vars[(n - 1)]*p$sigma_AR*p$p_F*(p$F / 2)
		new_larvaeT = state_vars[n]*p$sigma_AR*p$p_F*(p$F / 2)
	}


	new_states[1] = new_states[1] + new_larvaeR + new_larvaeT
	return(list(new_states=new_states, num_recruited=recruited))
}

build_fecundity_matrix = function(p){
	# Build the density-independent fecundity matrix
	#
	# Parameters
	# ----------
	# p : list
	# 	List of parameters needed in the fecundity matrix
	#
	# Returns
	# -------
	# : the fecundity matrix

	# Columns of matrix
	L1 = rep(0, 7)
	L2 = rep(0, 7)
	L3 = rep(0, 7)
	J1 = rep(0, 7)
	J2 = rep(0, 7)
	AR = c(p$sigma_AR*p$p_F*p$F, 0, 0, 0, 0, 0, 0)
	AT = c(p$sigma_AT*p$p_F*p$F, 0, 0, 0, 0, 0, 0)

	F = cbind(L1, L2, L3, J1, J2, AR, AT)
	return(F)

}

stochastic_simulation = function(steps, initial_values, muR, muT, phi, phi_sigma_J1, 
								 new_params=list(), sigma_J1_traj=NA){
	# Stochastic simulation of projection matrix
	#
	# Parameters
	# ----------
	# steps : int
	# 	yearly time steps to simulate over
	# initial_values : vector
	#		Vector of length 6 specifying 
	# phi : float
	#		Dispersion parameter of a beta distribution
	# muR : float
	#		Average recruited adult survival probability
	# muT : float
	#		Average translocated adult survival probability
	# mu_sigma_J1 : float
	# 	Mean successful recruitment probability
	# phi_sigma_J1 : float
	#		Dispersion of successful recruitment probability
	# new_params : list
	#		Any new parameter values to assign to params
	# sigma_J1_traj: array or NA
	# 	If array, you can specifiy a fixed trajectory.
	#
	# Returns
	# -------
	# : list
	#	results: projected state variables
	#	num_recruited: number adults recruited in a time step
	#	sigma_J1_probs: Probability of recruitment in a time step

	results = array(NA, dim=c(length(initial_values), steps + 1))
	results[, 1] = initial_values
	num_recruited = array(NA, dim=steps)
	sigma_J1_probs = array(NA, dim=steps)
	# adult_surv_values = array(NA, dim=steps)

	aR = muR*phi
	aT = muT*phi

	bR = phi - aR
	bT = phi - aT

	for(t in 2:(steps + 1)){

		values_now = results[, t - 1]

		# Stochastically update adult survival based on empirical data
		adult_survR = rbeta(1, aR, bR)
		adult_survT = rbeta(1, aT, bT)

		# adult_surv_values[t - 1] = adult_surv
		params = assign_parameters(adult_survR, adult_survT, new_params=new_params)

		# Update year 1 juvenile survival stochastically 
		if(all(is.na(sigma_J1_traj))) {
			a_sigma_J1 = params$sigma_J1 * phi_sigma_J1
			b_sigma_J1 = phi_sigma_J1 - a_sigma_J1
			params$sigma_J1 = rbeta(1, a_sigma_J1, b_sigma_J1)
		} else{
			params$sigma_J1 = sigma_J1_traj[t - 1]
		}

		# Stochastically update population
		stoch_step = stochastic_step(values_now, params)
		results[, t] = stoch_step$new_states

		num_recruited[t - 1] = stoch_step$num_recruited
		sigma_J1_probs[t - 1] = params$sigma_J1
	}

	return(list(results=results, num_recruited=num_recruited, sigma_J1_probs=sigma_J1_probs))
}

build_projection_matrix = function(params){
	# Build the projection matrix
	#
	# Parameters
	# ----------
	# params : list
	# 	List of model parameters
	# 
	# Return
	# ------
	# Full projection matrix

	T = build_transition_matrix(params)
	F = build_fecundity_matrix(params)
	P = T + F
	return(P)
}

# We also need the derivatives with respect to parameters
get_sensitivity = function(params){
	# Calculates sensitivites of lambda to sigma_AR, F, sigma_J1, sigma_J2
	#
	# Parameters
	# ----------
	# params : list
	#
	# Returns
	# -------
	# : list
	#	sens: length 5 array with sensitivities of lambda to sigmaAR, omega, F, sigma_J1, sigma_J2
	# 	elas: Length 5 array with elasticities of lambda to sigmaAR, omega, F, sigma_J1, sigma_J2	

	A = build_projection_matrix(params)

	# Get max eigenval eigenvect
	# Get left and right eigen vect
	right = eigen(A)
	left = eigen(t(A))
	lam = max(Re(right$values))

	# Get left and right eigenvectors corresponding to
	# maximum eigenvalue. Vectors already normalized
	w = Re(right$vectors[, which.max(Re(right$values)), drop=F])
	v = Re(left$vectors[, which.max(Re(right$values)), drop=F])

	# Parameters: [sigma_AR, omega, F]. From Caswell 2019, 3.42
	part1 = kronecker(t(w), t(v)) / (t(v) %*% w)[1]
	
	# Derivatives with respect to sigma_AR, omega, F*pF, sigma_J1, sigma_J2
	dAdsar = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, params$F*params$p_F, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
	dAdF = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, params$p_F*params$sigma_AR, 0, 0, 0, 0, 0, 0, params$p_F*params$sigma_AT, 0, 0, 0, 0, 0, 0) 
	dAdsj1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, params$p_J1, (1 - params$p_J1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
	dAdsj2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
	dmatA = t(rbind(dAdsar, dAdF, dAdsj1, dAdsj2))

	# From Caswell 2019
	D_sens_params = diag(c(params$sigma_AR, params$F, params$sigma_J1, params$sigma_J2))
	sens = part1 %*% dmatA
	elas = (sens %*% D_sens_params) / lam

	return(list(sens=sens, elas=elas))

}

get_lambda = function(params){
	# From params, compute deterministic growth rate
	#
	# Return
	# ------
	# : lambda

	P = build_projection_matrix(params)
	lam = max(abs(eigen(P)$values))
	return(lam)

}


abundances_from_model <- function(model, what) {
	# Plot a mark-recapture model's output.  Adapted from plot_model in the mrmr package
	#
	# Parameters
	# ----------
	# model : Model object 
	#	A model objectgenerated by \code{mrmr:fit_model()}.
	# what : str 
	#	What to plot. Must be a character string of either "abundance" or "recruitment"
	# 
	# Return
	# ------
	# : data.table
	#	A data.table with estimates of abundance or recruitment for each primary period

	valid_plots <- c("abundance", "recruitment")
	stopifnot(what %in% valid_plots)
	any_translocations <- 'data.frame' %in% class(model$data$translocations)

	survey_prim_periods <- model$data$surveys %>%
	  group_by(.data$primary_period) %>%
	  filter(.data$secondary_period == min(.data$secondary_period))

	survey_prim_periods <- ungroup(survey_prim_periods)

	if(what == "abundance"){
		p <- model$m_fit$draws("N", format = "draws_df") %>%
			tidyr::pivot_longer(tidyselect::starts_with("N")) %>%
			suppressWarnings() %>%
			mutate(
		  		get_numeric_indices(.data$name),
		  		primary_period = .data$index_1 + 1
				  ) %>%
			group_by(.data$primary_period) %>%
			summarize(lo = quantile(.data$value, .025),
		          	  med = median(.data$value),
		          	  hi = quantile(.data$value, .975),
		              .groups = "drop") %>%
			left_join(survey_prim_periods, by = "primary_period")
	}

	if(what == "recruitment"){
	  p <- model$m_fit$draws("B", format = "draws_df") %>%
	    tidyr::pivot_longer(tidyselect::starts_with("B")) %>%
	    suppressWarnings() %>%
	    mutate(
	      get_numeric_indices(.data$name),
	      primary_period = .data$index_1 + 1
	    ) %>%
	    group_by(.data$primary_period) %>%
	    summarize(lo = quantile(.data$value, .025),
	              med = median(.data$value),
	              hi = quantile(.data$value, .975),
	              .groups = "drop") %>%
	    left_join(survey_prim_periods, by = "primary_period")
	}

	return(as.data.table(p))
}


# Code from MRMR package with slight modifications

get_numeric_indices <- function(string) {
  idx_mat <- stringr::str_extract_all(string, "[0-9]+", simplify = TRUE)
  colnames(idx_mat) <- paste0("index_", seq_len(ncol(idx_mat)))
  tibble::as_tibble(idx_mat) %>%
    dplyr::mutate_all(readr::parse_integer)
}

mean_na = function(x){

  if(all(is.na(x))) {
    return(NA)
  }else{
    return(mean(x, na.rm=TRUE))
  }
}

get_surv_vector = function(primary_period, value){

  val1 = value[primary_period < max(primary_period)]
  val2 = value[primary_period > min(primary_period)]

  # Survived 
  died = (val1 == 2) & (val2 == 3) # Died
  survived =  (val1 == 2) & (val2 == 2) # Survived
  died_vect = array(NA, length(died))
  died_vect[died] = 0
  died_vect[survived] = 1
  return(as.list(died_vect))

}

convert_to_death_rates = function(surv_probs, deltat){
  # Convert survival probabilities to death rates given
  # deltat.  Assuming geometric survival proability
  #
  # Parameters
  # ----------
  # surv_probs : array-like
  #   Survival probabilities
  # deltat : array-like
  #   Time period over which survival probabiliites are calculated
  #
  # Returns
  # -------
  # : death ratess

  rates = -log(surv_probs) / deltat
  return(rates)

}

survival_table <- function(model, thin=1, lower_upper=c(0.05, 0.95)) {

  any_translocations <- 'data.frame' %in% class(model$data$translocations)

  if (!any_translocations) {
    stop(paste("No translocation data are present, so a cohort survival",
               "table cannot be created."))
  }

   primary_period_dates <- model$data$surveys %>%
    group_by(.data$primary_period) %>%
    summarize(date = min(.data$survey_date),
              year = min(.data$year)) %>%
    ungroup

  survival_summary <- model$m_fit$draws("s", format = "draws_df") %>%
    posterior::thin_draws(thin) %>%
    tidyr::pivot_longer(tidyselect::starts_with("s")) %>%
    suppressWarnings() %>%
    mutate(
      get_numeric_indices(.data$name),
      primary_period = .data$index_2,
      pit_tag_id = dimnames(model$data$stan_d$Y)[[1]][.data$index_1]
    ) %>%
    filter(!(.data$pit_tag_id %in% as.character(model$data$translocations$pit_tag_id))) %>%
    left_join(primary_period_dates) %>%
    dplyr::transmute(
      .data$.draw, .data$value, .data$primary_period, .data$pit_tag_id,
      .data$date, .data$year
    ) 

  # Drop all augmented data and just extract survival probabiliites from
  # pit tagged, non-translocated individuals.
  ss_dt = data.table(survival_summary)[!(pit_tag_id %like% "aug")]
  num_pittags = length(ss_dt$pit_tag_id %>% unique())

  if(num_pittags > 0){

    ss_dt_id = ss_dt[order(pit_tag_id, .draw, primary_period)]
    prob_dt = ss_dt_id[, .(surv=get_surv_vector(primary_period, value),
                           primary_period=primary_period[1:(length(primary_period) - 1)]), by=.(.draw, pit_tag_id)]

    # Compute the survival probability conditional on being alive for
    # each primary period.
    pp_dates = as.data.table(primary_period_dates)
    diff_dates = pp_dates[, .(days=date[2:length(date)] - date[1:(length(date) - 1)],
                              primary_period=primary_period[1:(length(date) - 1)])]

    # Upper lower 
    nontrans_surv = prob_dt[, .(surv_prob=as.numeric(mean_na(as.integer(surv)))), by=.(.draw, primary_period)][, 
                            .(lower=quantile(surv_prob, lower_upper[1], na.rm=T),
                              med=quantile(surv_prob, 0.5, na.rm=T),
                              upper=quantile(surv_prob, lower_upper[2], na.rm=T)), by=.(primary_period)] %>% 
                    merge(diff_dates, by="primary_period")

    # nontrans_surv = prob_dt[, .(surv_prob=as.numeric(mean_na(as.integer(surv)))), by=.(primary_period)] %>% 
    #                 merge(diff_dates, by="primary_period")
    
    # Convert surv probs to death rates
    nontrans_surv$death_rate_day_med = convert_to_death_rates(nontrans_surv$med, 
                                                              as.integer(nontrans_surv$days))
    nontrans_surv$death_rate_day_lower = convert_to_death_rates(nontrans_surv$lower, 
                                                              as.integer(nontrans_surv$days))
    nontrans_surv$death_rate_day_upper = convert_to_death_rates(nontrans_surv$upper, 
                                                              as.integer(nontrans_surv$days))
    nontrans_surv$days = as.integer(nontrans_surv$days)

    res = list(num_pittags=num_pittags, surv_probs=nontrans_surv[!is.na(med)])
  } else {
    res = list(num_pittags=num_pittags, surv_probs=NA)
  }

  return(res)

}


