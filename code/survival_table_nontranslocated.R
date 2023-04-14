library(tidyverse)

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

survival_table <- function(model, thin=1, lowerupper=c(0.25, 0.75), translocated=FALSE) {
  # Take in MRMR model and return survival probabilities for pit-tagged 
  # individuals that are either translocated or naturall recruited
  #
  # Parameters
  # ----------
  # model : MRMR fitted model
  # thin : int, thin the bayesian chains
  # lowerupper : vector, the lower and upper bounds for survival estimates
  # translocated : bool, if TRUE only focus on translocated individuals, 
  #               else focus on naturally recruited individuals
  #
  # Returns
  # -------
  # : list 
  #   number of pit tags in population, survival probabilities of non-translocated

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
    )


  if (translocated) {

    # Get individuals that were translocated
    survival_summary <- survival_summary %>% filter(.data$pit_tag_id %in% as.character(model$data$translocations$pit_tag_id)) %>%
      left_join(primary_period_dates) %>%
      dplyr::transmute(
        .data$.draw, .data$value, .data$primary_period, .data$pit_tag_id,
        .data$date, .data$year
      ) 

  } else {

    # Get individuals that weren't translocated and were naturally recruited
    survival_summary <- survival_summary %>% filter(!(.data$pit_tag_id %in% as.character(model$data$translocations$pit_tag_id))) %>%
      left_join(primary_period_dates) %>%
      dplyr::transmute(
        .data$.draw, .data$value, .data$primary_period, .data$pit_tag_id,
        .data$date, .data$year
      ) 

  }

  # Drop all augmented data
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
    est_surv = prob_dt[, .(surv_prob=as.numeric(mean_na(as.integer(surv))), num_obs=sum(!is.na(surv))), by=.(.draw, primary_period)][, 
                            .(lower=quantile(surv_prob, lowerupper[1], na.rm=T),
                              med=quantile(surv_prob, 0.5, na.rm=T),
                              upper=quantile(surv_prob, lowerupper[2], na.rm=T),
                              num_obs=round(mean(num_obs))), by=.(primary_period)] %>% 
                    merge(diff_dates, by="primary_period")

    # Convert surv probs to death rates
    est_surv$death_rate_day_med = convert_to_death_rates(est_surv$med, 
                                                              as.integer(est_surv$days))
    est_surv$death_rate_day_lower = convert_to_death_rates(est_surv$lower, 
                                                              as.integer(est_surv$days))
    est_surv$death_rate_day_upper = convert_to_death_rates(est_surv$upper, 
                                                              as.integer(est_surv$days))
    est_surv$days = as.integer(est_surv$days)

    res = list(num_pittags=num_pittags, surv_probs=est_surv[!is.na(med)])
  } else {
    res = list(num_pittags=num_pittags, surv_probs=NA)
  }

  return(res)

}

build_surv_list = function(res, drop_ones){
  # Build the time averaged survival probabilities for each lake
  #
  # Parameters
  # -----------
  # res : result from survival_table
  # drop_ones : bool, if TRUE, drop the initial string of ones. 
  #             If FALSE, keep the initial string of ones.


  if(res$num_pittags > 0) { 

    surv_probs = res$surv_probs

    # Find the first stretch of primary periods where survival is 1 and remove them
    # Doing this because individuals were not actually observed during this period
    # and it is always conditional on recruitment.  Thus, survival has to be 1 
    # as we know they were alive later.

    nm_append = c("lower", "med", "upper")

    if(drop_ones){
      rle_res = rle(surv_probs[["med"]] == 1)
      cut = ifelse(rle_res$values[1] == TRUE, rle_res$lengths[1], 0)
      surv_probs = surv_probs[(cut + 1):nrow(surv_probs)]
    }


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

  return(surv_list)
  
}
