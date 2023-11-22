library(data.table)
source("population_viability_functions.R")
library(ggplot2)
library(patchwork)
library(parallel)

#### Data loading and preliminary analysis ######

# Load Stan fits and lake ids.  The stan fits contain a model that has analyzed
# The MRMR results of yearly adult survival probability.
fit = readRDS(file.path("..", "out", "adult_survival_probability_model.rds"))
keep_lakes = readRDS(file.path("..", "out", "lakes_to_use_in_analysis.rds"))

colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f',
					 '#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

# Color blind-friendly palette
# colors = c("#000000","#004949","#009292","#ff6db6","#ffb6db",
#  "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
#  "#920000","#924900","#db6d00","#24ff24","#ffff6d")

# Extract lake-specific probabilities of translocated survival
beta0 = rstan::extract(fit, "beta0")$beta0
alpha = rstan::extract(fit, "alpha")$alpha
mus = apply(alpha, 2, function(x) 1 / (1 + exp(-(beta0 + x))))
colnames(mus) = keep_lakes
mu_dist = melt(data.table(mus), variable.name = "lake_id", value.name="surv_est")
mean_mu_dist = mu_dist[, .(mean_est=mean(surv_est)), by=.(lake_id)]

phi = exp(rstan::extract(fit, "logphi")$logphi)
phi = 10000 # Set to be large to ignore inter-annual variability in survival

# Load in survival estimates from all individuals
all_surv = fread(file.path("..",
						   "out", 
						   "yearly_survival_estimates_all_individuals.csv")
				)[order(lake_id)]

# Load survival probabilities estimated for recruited individuals
recruited_surv = fread(file.path("..",
						   "out", 
						   "yearly_survival_estimates_only_recruited_individuals.csv")
				)[order(lake_id)]

#############################################################
###########  Analysis 1: Lambda calculation #################
#############################################################

# Calculate deterministic growth rates for site

lake_growth_rates = list()
nm_append = c("med", "lower", "upper")
quant_vals = c(0.5, 0.025, 0.975)
med_surv_vals = array(NA, dim=nrow(all_surv))

for(l in 1:nrow(all_surv)){

	lake = all_surv$lake_id[l]

	lambdas_mlu = array(NA, dim=3)
	for(j in 1:length(nm_append)){

		nm = nm_append[j]

		muT = quantile(1 / (1 + exp(-(beta0 + alpha[, l]))), quant_vals[j])

		if(is.na(recruited_surv[lake_id == lake]$surv_prob_med) | lake %in% c(70641, 70556)){
			# For these populations, there is no estimate of naturally recruited survival
			# or the estimate is highly variable (i.e., between 0 and 1)
			muR = muT
		} else{
			muR = recruited_surv[lake_id == lake][[paste0('surv_prob_', nm)]]
		}

		if(nm == "med")
			med_surv_vals[l] = muR

		params = assign_parameters(muR, muT)
		tlam = get_lambda(params)
		lambdas_mlu[j] = tlam
	}

	lake_growth_rates[[paste0(lake)]] = lambdas_mlu
}

lake_growth_rates_dt = data.table(do.call(rbind, lake_growth_rates))
colnames(lake_growth_rates_dt) = nm_append
lake_growth_rates_dt$lake_id = names(lake_growth_rates)
lake_growth_rates_dt$surv_med = med_surv_vals

p1 = ggplot(data=lake_growth_rates_dt) + geom_point(aes(x=lake_id, y=med, color=lake_id)) +
		   						    geom_errorbar(aes(x=lake_id, ymin=lower, ymax=upper, color=lake_id), width=0.25) +
		   						    scale_color_manual(values=colors) +
		   						    geom_hline(aes(yintercept=1), linetype="dashed") +
		   						    theme_classic() + ylab("Long-run growth rate, lambda") +
		   						    xlab("Site ID")

ggsave(file.path("..", "out", "lambda_estimates.pdf"), width=8, height=5)

#############################################################
###########  Analysis 1.5: Lambda sensitivity ###############
#############################################################


muR_range = seq(0.01, 0.9, len=50)
juv_surv_range = seq(0.01, .25, len=50)
lambdas = array(NA, dim=c(length(muR_range), length(juv_surv_range)))

for(i in 1:length(muR_range)){
	for(j in 1:length(juv_surv_range)){

		muR = muR_range[i]
		juv_surv = juv_surv_range[j]
		params = assign_parameters(muR, muR, new_params=list(sigma_J1=juv_surv))
		lambdas[i, j] = get_lambda(params)

	}
}

dat = expand.grid(muR=muR_range, juv_surv=juv_surv_range)
dat$lambda = as.vector(lambdas)
ptile = ggplot(dat) + geom_tile(aes(x=muR, y=juv_surv, fill=lambda)) + 
							geom_hline(aes(yintercept=0.08), linetype="dashed") +
							geom_contour(aes(x=muR, y=juv_surv, z=lambda), breaks=c(1), color="red", size=1) +
							geom_point(data=lake_growth_rates_dt, aes(x=surv_med, 
																												y=rep(0.08, length(surv_med)), 
																													color=lake_id), size=2) +
		   				scale_color_manual(values=colors) +
							# scale_fill_gradient("\u03bb", low = "white", high = "black") +
							scale_fill_gradient2(midpoint=1) + #guide=guide_legend(title.position="bottom")) +
							theme_classic() + xlab(bquote("Adult survival probability (\u03c3"~.[AR]~")")) +
							ylab(bquote("Juvenile (Y1) survival probability (\u03c3"~.[J1]~")")) + 
							annotate("text", x=0.25, y=0.05, label="Population declines", size=3) +
							annotate("text", x=0.7, y=0.2, label="Population grows", size=3) +
							guides(color="none", fill=guide_colorbar(title="\u03bb")) + theme(legend.position="bottom", 
																					 legend.text=element_text(size=7), 
																					 legend.title=element_text(size=7))


#############################################################
###########  Analysis 2: Lambda sensitivity #################
#############################################################

# How sensitive is lambda to a change in key parameters?

lake_sens = list()
lake_elas = list()
nm_append = c("med")
quant_vals = c(0.5, 0.025, 0.975)

for(l in 1:nrow(all_surv)){

	lake = all_surv$lake_id[l]

	muT = quantile(1 / (1 + exp(-(beta0 + alpha[, l]))), 0.5)

	if(is.na(recruited_surv[lake_id == lake]$surv_prob_med) | lake %in% c(70641, 70556)){
		# For these populations, there is no estimate of naturally recruited survival
		muR = muT
	} else{
		muR = recruited_surv[lake_id == lake][[paste0('surv_prob_', 'med')]]
	}

	params = assign_parameters(muR, muT)
	sens_elas = get_sensitivity(params)

	lake_sens[[paste0(lake)]] = sens_elas$sens
	lake_elas[[paste0(lake)]] = sens_elas$elas
}

lake_elas_dt = data.table(do.call(rbind, lake_elas))
colnames(lake_elas_dt) = c("\U03C3_AR", "F", "\U03C3_J1", "\U03C3_J2") #c("sigma_AR", "F", "sigma_J1", "sigma_J2")
lake_elas_dt$lake = names(lake_elas)

lake_elas_dt_melt = melt(lake_elas_dt, id.vars="lake")

p2 = ggplot(lake_elas_dt_melt) + geom_bar(aes(x=variable, y=value, fill=lake), stat="identity", position="dodge") +
														scale_fill_manual(values=colors) +
														xlab("Parameter") + ylab("Elasticity of \U03BB") + theme_classic() +
														guides(fill=guide_legend(title="Site"))
														# theme(legend.position = c(0.7, 0.7))

#############################################################
###########  Analysis 3: Extinction curves ##################
#############################################################

#### Run extinction curve simulation ####

sims = 1000 # number of replicates
steps = 50 # years
beta0 = rstan::extract(fit, "beta0")$beta0
alpha = rstan::extract(fit, "alpha")$alpha
keep_lakes_new = c(keep_lakes)#, "average")
phi_omega = 2 # Inverse dispersion of juvenile dispersion and recruitment probability
omegas = seq(0, 1, by=0.05) 
extinction_at_fifty = list()

# Specify whether we want to allow recruited adult survival probability to differ
# from translocated adult survival probability
run_type = "Recruit survival != Translocated survival"

extinction_curve_for_omega = function(omega, surv_juv){

	cat("Working on omega =", omega, "\n")

	compare_muT_mR = list()
  extinction_curves = list()
  pop_trajectories = list()
  recruitment_trajectories = list()
  omega_trajectories = list()
  
  for(l in 1:length(keep_lakes_new)){
    
    lake = keep_lakes_new[l]

    if(lake %in% c(keep_lakes_new)){

	    cat("Working on", lake, "\n")
	    
	    extinction_array = array(NA, dim=c(sims, steps + 1))
	    recruitment_array = array(NA, dim=c(sims, steps))
	    omega_array = array(NA, dim=c(sims, steps))
	    total_pop_array = array(NA, dim=c(sims, steps + 1))
	    
	    # Mu and phi parameters on beta distribution of adult survival
	    if(lake != "average"){

	      muT =  median(1 / (1 + exp(-(beta0 + alpha[, l]))))
	      lower_muT = quantile(1 / (1 + exp(-(beta0 + alpha[, l]))), c(0.05))
	      upper_muT = quantile(1 / (1 + exp(-(beta0 + alpha[, l]))), c(0.95))
	      
	      if(run_type == "Recruit survival == Translocated survival"){
	        muR = muT
	      }
	      else{

	      	if(is.na(recruited_surv[lake_id == lake]$surv_prob_med) | lake %in% c(70641, 70556)){

	      		# Lakes don't really have recruited individuals
	      		# So set to translocated survival
	        	# muR = all_surv[lake_id == lake]$surv_med
	        	muR = muT
	       		lower_muR = muT
	       		upper_muR = muT

	      	} else{
	      		muR = recruited_surv[lake_id == lake]$surv_prob_med
	      		lower_muR = recruited_surv[lake_id == lake]$surv_prob_lower
	      		upper_muR = recruited_surv[lake_id == lake]$surv_prob_upper
	      	}
	      }
	      
	    } else{
	      muT =  median(1 / (1 + exp(-(beta0))))
	      
	      if(run_type == "Recruit survival == Translocated survival"){
	        muR = muT
        	lower_muR = muT
       		upper_muR = muT
	      }
	      else{
	        muR = mean(all_surv$surv_med)
        	lower_muR = muT
       		upper_muR = muT
	      }
	    }
	    phi = 10000 #median(phi)
	   	
	   	if(run_type != "Recruit survival == Translocated survival"){ 
	    	compare_muT_mR[[as.character(lake)]] = c(muR=muR, 
	    																					 lower_muR=lower_muR,
	    																					 upper_muR=upper_muR,
	    																					 muT=muT, 
	    																					 lower_muT=lower_muT, 
	    																					 upper_muT=upper_muT)
	    }

	    # Set initial conditions
			recruitment = fread(file.path("..", "data", "clean", "abundance_and_recruitment", paste0(lake, "_recruitment.csv")))
			med_recruit = recruitment[, .(recruit=sum(med)), by=.(year)]
			med_recruit[, index:=.(1:nrow(med_recruit))]
			initial_values = c(0, 0, 0, 0, 0, 0, 40) #med_recruit$recruit[1]) # Initial introduction size
			omega_traj = NA # No prespecified omega trajectory

	    # Run simulations
	    for(k in 1:sims){

	      sim_res = stochastic_simulation(steps, initial_values, muR, muT, phi, phi_omega, 
	      																new_params=list(omega=omega, reproduction="nbd", sigma_J1=surv_juv),
	      																omega_traj=omega_traj)
	      results = sim_res$results
	      num_recruited = sim_res$num_recruited
	      extinction_array[k, ] = as.integer(apply(results, 2, function(x) all(x == 0)))
	      recruitment_array[k, ] = num_recruited
	      omega_array[k, ] = sim_res$omega_probs
	      total_pop_array[k, ] = apply(results[6:7, ], 2, sum) # Just look at adults
	    }
	    
	    extinction_curves[[lake]] = data.frame(ext_prob=colMeans(extinction_array), lake_id=lake, time=0:steps)
	    pop_trajectories[[lake]] = total_pop_array
	    recruitment_trajectories[[lake]] = recruitment_array
	    omega_trajectories[[lake]] = omega_array

	  } ## End lake loop
    
  }
  
  extinction_curves_dt = as.data.table(do.call(rbind, extinction_curves))

  extdt = extinction_curves_dt[time == 50]
  extdt[, omega:=omega]
  return(list(extinction=extdt, compare=compare_muT_mR))
	  

} # End function


# Run for different year 1 juvenile survivals
baseline_surv = c(0.09, 0.25)
for(tsurv in baseline_surv){

	# Parallelize calculation of extinction curves
	parallel_res = mclapply(omegas, function(x) extinction_curve_for_omega(x, surv_juv=tsurv), mc.cores=8)
	extinction_at_fifty = lapply(parallel_res, function(x) x[['extinction']])
	compare_muT_mR = lapply(parallel_res, function(x) x[['compare']])[[1]]

	extinction_at_fifty_dt = do.call(rbind, extinction_at_fifty) %>% 
													 merge(lake_growth_rates_dt[, .(lake_id, surv_med)], key="lake_id")
	extinction_at_fifty_dt[, lake_id_surv:=paste0(lake_id, ", ", round(surv_med, 2))]


	# Set up colors for plotting
	cdat = data.table(extinction_at_fifty_dt)[, .(lake_id_surv=lake_id_surv[1],
																								surv_med=surv_med[1]), 
																						by=.(lake_id)][order(lake_id)]
	cdat[, col:=colors]
	cdat = cdat[order(surv_med)]
	extinction_at_fifty_dt$lake_id_surv = factor(extinction_at_fifty_dt$lake_id_surv, levels=cdat$lake_id_surv)


	pext = ggplot(extinction_at_fifty_dt[order(surv_med, lake_id)], aes(x=100*(1 - omega), y=ext_prob)) + 
								#geom_vline(aes(xintercept=1 - 0.5), linetype="dashed") +
								geom_line(aes(color=lake_id_surv)) + 
								geom_point(aes(color=lake_id_surv)) + 
								scale_color_manual(values=cdat$col) +
								theme_classic() + xlab("Average % reduction in juvenile survival") + ylab("Extinction prob. in 50 years") +
								guides(color=guide_legend(title=bquote("Site, \u03c3"~.[AR]), title.position="top", ncol=3)) + 
								theme(legend.position="bottom", legend.text=element_text(size=7), legend.title=element_text(size=8))

	if(tsurv == 0.09){
		# Only save the extinction curves for different baseline survival
		ggsave(file.path("..", "out", "extinction_surves_surv=0.09.jpg"), plot=pext, width=5, height=5, dpi=300)
	}

}

# Save comparison between translocated and survival probabilities from all individuals

compare_muT_mR_mat = do.call(rbind, compare_muT_mR)
compare_muT_mR_dt = data.table(compare_muT_mR_mat)
compare_muT_mR_dt$lake_id = row.names(compare_muT_mR_mat)
colnames(compare_muT_mR_dt) = c("yearly_surv_prob_recruit", 
																"yearly_surv_prob_recruit_lower",
																"yearly_surv_prob_recruit_upper",
																"yearly_surv_prob_translocated", 
																"yearly_surv_prob_translocated_lower", 
																"yearly_surv_prob_translocated_upper", 
																"lake_id")
fwrite(compare_muT_mR_dt, file.path("..", "out", "compare_surv_probs.csv"))

#############################################################
###########  Analysis 4: Model projections ##################
#############################################################

# Compare the projected model trajectories to the observed trajectories
# Just focus on 70550

sims = 50000 # number of replicates
steps = 20 # years
beta0 = rstan::extract(fit, "beta0")$beta0
alpha = rstan::extract(fit, "alpha")$alpha
keep_lakes_new = c(keep_lakes)#, "average")
plots = list()
phi_omega = 2 # Inverse dispersion of recruitment probability
omega = 0.5

extinction_curves = list()
pop_trajectories = list()
recruitment_trajectories = list()
omega_trajectories = list()
sigma_J1_results = list()

for(l in 1:length(keep_lakes_new)){
  
  lake = keep_lakes_new[l]

  if(lake %in% c("70550")){

	  cat("Working on", lake, "\n")
	  
	  extinction_array = array(NA, dim=c(sims, steps + 1))
	  recruitment_array = array(NA, dim=c(sims, steps))
	  omega_array = array(NA, dim=c(sims, steps))
	  total_pop_array = array(NA, dim=c(sims, steps + 1))
	  sigma_J1_array = array(NA, dim=c(sims))
	  
	  muT =  median(1 / (1 + exp(-(beta0 + alpha[, l]))))
	  muR = recruited_surv[lake_id == lake]$surv_prob_med
	  phi = 10000  # Set as large so there is no yearly variation in survival probability

	  # Set initial conditions
		recruitment = fread(file.path("..", "data", "clean", "abundance_and_recruitment", paste0(lake, "_recruitment.csv")))
		med_recruit = recruitment[, .(recruit=sum(med)), by=.(year)]
		med_recruit[, index:=.(1:nrow(med_recruit))]
		initial_values = c(0, 0, 0, 0, 0, 0, med_recruit$recruit[1]) # Initial introduction size
		omega_traj = NA # No prespecified omega trajectory

	  # Run simulations
	  for(k in 1:sims){

	  	tsigma_J1 = runif(1, 0, 0.25) # Prior draw for sigma_J1
	  	sigma_J1_array[k] = tsigma_J1
	    sim_res = stochastic_simulation(steps, initial_values, muR, muT, phi, phi_omega, 
	    																new_params=list(omega=omega, reproduction="nbd", sigma_J1=tsigma_J1),
	    																omega_traj=omega_traj)
	    results = sim_res$results
	    num_recruited = sim_res$num_recruited
	    extinction_array[k, ] = as.integer(apply(results, 2, function(x) all(x == 0)))
	    recruitment_array[k, ] = num_recruited
	    omega_array[k, ] = sim_res$omega_probs
	    total_pop_array[k, ] = apply(results[6:7, ], 2, sum) # Just look at adults
	  }
	  
	  extinction_curves[[lake]] = data.frame(ext_prob=colMeans(extinction_array), lake_id=lake, time=0:steps)
	  pop_trajectories[[lake]] = total_pop_array
	  recruitment_trajectories[[lake]] = recruitment_array
	  omega_trajectories[[lake]] = omega_array
	  sigma_J1_results[[lake]] = sigma_J1_array


	} ## End lake loop
  
}
  
# Look at predictions for lake 70550
lake = '70550'

# Extract obs abundance and recruitment data
abundance = fread(file.path("..", "data", "clean", "abundance_and_recruitment", paste0(lake, "_abundance.csv")))
recruitment = fread(file.path("..", "data", "clean", "abundance_and_recruitment", paste0(lake, "_recruitment.csv")))
med_abund = abundance[, .(abund=mean(med)), by=.(year)]
med_abund[, index:=.(1:nrow(med_abund))]
med_recruit = recruitment[, .(recruit=sum(med)), by=.(year)]
med_recruit[, index:=.(1:nrow(med_recruit))]

# Compare observed and predicted abundance using sum of squares
obs = med_abund$abund
pred_ss = apply(pop_trajectories[['70550']][1:sims, 1:16], 1, function(x) sum((obs - x)^2))
pred_ss_dt = data.table(ss=pred_ss, index=1:length(pred_ss))[order(ss)]
ind = which.min(pred_ss)

# Plot posterior distribution for omega
best_omegas = reshape2::melt(omega_trajectories[['70550']][pred_ss_dt$index[1:250], 1:16])
ggplot(best_omegas) + geom_histogram(aes(x=value), bins=15) + facet_wrap(~Var2)
colMeans(omega_trajectories[['70550']][pred_ss_dt$index[1:50], 1:16])

# Plot posterior distribution for sigma_J1
best_sigmas = sigma_J1_results[['70550']][pred_ss_dt$index[1:1000]]

# Plot and compare observed and predicted abundance
pred_traj = data.table(reshape2::melt(pop_trajectories[['70550']][pred_ss_dt$index[1:100], 1:16]))
colnames(pred_traj) = c("sim", "time", "abund")
map = data.table(year=2006:2021, time=1:16)
pred_traj = merge(pred_traj, map, by="time")
ptraj = ggplot() + geom_line(data=pred_traj, aes(x=year, y=abund, group=sim, color="Predicted"), alpha=0.25) +
										geom_line(data=NULL, aes(x=2006:2021, y=obs, color="Observed, 70550"), size=1) + 
										scale_color_manual(values=c(colors[7], "black")) +
										theme_classic() + xlab("Time (year)") + ylab("Adult abundance") +
										scale_x_continuous(breaks=seq(2006, 2021, 2)) + 
										theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
											   legend.position="bottom", legend.text=element_text(size=8), legend.title=element_blank())

# Plot and compare observed and predicted recruitment
# matplot(t(recruitment_trajectories[['70550']][pred_ss_dt$index[1:100], 1:15]), type="l", ylab="Recruitment", xlab="Time (years)")
# lines(1:15, med_recruit$recruit[2:16], col="red", lwd=5)

#########################################
##### Join all of the plots together ####
#########################################

myplot = (ptile + pext + ptraj) + plot_annotation(tag_levels="A", tag_suffix="")
ggsave(file.path("..", "out", "pop_viability_figures_for_manuscript.jpg"), width=11, height=5, dpi=300)

myplot = (p2) + plot_annotation(tag_levels="A", tag_suffix="")
ggsave(file.path("..", "out", "pop_viability_figures_for_supp.jpg"), width=5, height=4, dpi=300)

