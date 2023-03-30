# Translocations to reestablish populations of the mountain yellow-legged frog

## Authors of this repository

Roland A. Knapp (roland.knapp(at)ucsb.edu) [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--1954--2745-green.svg)](https://orcid.org/0000-0002-1954-2745)

Mark Q. Wilber (mqwilber(at)gmail.com) [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--8274--8025-green.svg)](https://orcid.org/0000-0002-8274-8025) 

Thomas C. Smith (tcsmith(at)ucsb.edu [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7908--438X-green.svg)](https://orcid.org/0000-0001-7908-438X)

## Overview of contents

This repository is organized as a reproducible research compendium, and describes frog translocations to reestablish populations of the endangered [mountain yellow-legged frog](https://www.fws.gov/sites/default/files/documents/Mountain-Yellow-Legged-Frog-Conservation-Strategy.pdf). 
It contains data on frog survival following translocations conducted during the period 2006-2020, code to identify predictors of survival using [R](https://www.r-project.org/), a notebook to document project progress to date, and a manuscript (in preparation). 
Throughout this repository, frog populations are referenced only by 5-digit unique site identifiers.
No site names or x-y coordinates are provided to protect these sensitive populations to the maximum extent possible.

This repository contains the following directories and files:

- `code/` directory: `Rmd` files that contain code to create and analyze the project datasets. 
- `data/` directory: Raw data and cleaned data.
- `doc/` directory: Manuscript and notebook files.
- `out/` directory: Output files.

## Notebooks

Notebooks are available in the [doc/](https://github.com/SNARL1/translocation/tree/main/doc/notebook#readme) and [out/notebooks_code/](https://github.com/SNARL1/translocation/tree/main/out/notebooks_code#readme) directories. Links are to README files that describe how to view notebooks directly from GitHub. The `doc/` notebook describes issues of interest related to dataset creation and analysis. The `out/notebooks_code/` notebooks are rendered from the `Rmd` files in the `code/` directory, and describe dataset creation/data analysis steps for each population. 

## Reproducing this analysis


### Preliminary instructions

1. Install the package manager anaconda or miniconda (https://docs.continuum.io/anaconda/install/).
2. Build the conda environment specified in the environmental.yml to install most of the R packages needed to run the analysis.  This can be done on the command line using the command
	- `conda env create -f environment.yml`.  
	- See https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file.
3. Activate the environment in the terminal
	- `conda activate r_env_translocation`
4. NOTE: You don't have to use the conda environment and can use your native R build. Just ensure you have the packages installed that are listed in the environmental.yml file.
 
### Running the Makefile

1. First, run `make install_packages` to ensure that you have the necessary packages installed in your R environment.
2. To run all of the analyses, type `make all` in the command line.
2. To run only the translocation analysis execute `make trans_analysis` in the command line.
3. To run only the viability analysis execute `make viability_analysis` in the command line.
4. To clean you directory of extraneous files execute `make clean` in the command line.

## Contact

Roland Knapp, Research Biologist, University of California Sierra Nevada Aquatic Research Laboratory, Mammoth Lakes, CA 93546 USA; rolandknapp(at)ucsb.edu, <https://mountainlakesresearch.com/roland-knapp/>
