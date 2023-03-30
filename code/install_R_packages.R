# Install the packages necessary for R analyses

# List of CRAN packages needed for the analysis.  Uncomment and run to load packages in your 
# native R environment. Otherwise, packages are assumed to be loaded through conda.
# list.of.packages <- c("data.table", "ggplot2", "patchwork", "cmdstanr", "magrittr",
# 					  "rstan", "lubridate", "tidyverse", 
# 					  "coda","mvtnorm","devtools","loo","dagitty", "tmap",
# 					  "RPostgreSQL", "DBI", "here", "brms", "rstanarm", 
# 					  "tidybayes", "broom.mixed", "janitor", "arm", "stars")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

# Load Github packages
github_packages = c("SNARL1/mrmr", "rmcelreath/rethinking")

for(package in github_packages){
	devtools::install_github(package)
}
