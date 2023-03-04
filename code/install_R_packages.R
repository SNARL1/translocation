# Install the packages necessary for R analyses

list.of.packages <- c("data.table", "ggplot2", "patchwork", "cmdstanr", "magrittr",
					  "rstan", "lubridate", "tidyverse", 
					  "coda","mvtnorm","devtools","loo","dagitty")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load Github packages
github_packages = c("rethinking", "mrmr")

for(package in github_packages){

	if(!(package %in% installed.packages()[,"Package"])){
		if(package == "rethinking"){
			devtools::install_github("rmcelreath/rethinking")
		} else if(package == "mrmr"){
			devtools::install_github("SNARL1/mrmr")
		} 
	}
}
