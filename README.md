# Translocations to reestablish populations of the mountain yellow-legged frog

## Authors of this repository

Roland A. Knapp (roland.knapp(at)ucsb.edu) [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--1954--2745-green.svg)](https://orcid.org/0000-0002-1954-2745)

Mark Q. Wilber (mqwilber(at)gmail.com) [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--8274--8025-green.svg)](https://orcid.org/0000-0002-8274-8025) 

## Overview of contents

This repository is organized as a reproducible research compendium, and describes frog translocations to reestablish populations of the endangered [mountain yellow-legged frog](https://www.fws.gov/sites/default/files/documents/Mountain-Yellow-Legged-Frog-Conservation-Strategy.pdf). 
It contains data on frog survival following translocations conducted during the period 2006-2020, [R](https://www.r-project.org/) code to identify predictors of survival and estimate population viability, a notebook to document translocation-related analyses, and a manuscript (submitted). 
Throughout this repository, frog populations are referenced only by 5-digit unique site identifiers.
No site names or x-y coordinates are provided to protect these sensitive populations to the maximum extent possible.

This repository contains the following directories and files:

- `code/` directory: `Rmd` and `qmd` files to create and analyze the project datasets, and create figures. 
- `data/` directory: Raw and cleaned data, and map data/layers.
- `doc/` directory: Manuscript and notebook files.
- `out/` directory: Output files, including files related to the viability analyses.

## Notebooks

A notebook describing dataset creation and analysis is available in the [doc/notebook](https://github.com/SNARL1/translocation/tree/main/doc/notebook#readme) directory. Link is to a README file that describes how to view the notebook directly from GitHub.

## Reproducing this analysis

### Instructions

1. Install the package manager anaconda or miniconda (https://docs.continuum.io/anaconda/install/).
2. Build the conda environment specified in the environmental.yml to install most of the R packages needed to run the analysis.  This can be done on the command line using the command
	- `conda env create -f environment.yml`.  
	- See [link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) for details about conda enviroments.
3. Activate the environment in the terminal
	- `conda activate r_env_translocation`
4. NOTE: You don't have to use the conda environment and can instead use your native R build. Just ensure you have the packages installed that are listed in the `environment.yml` file.
	- For Linux, we have found that to install the R package `devtools` you are often prompted to install additional software outside of R.  Linux typically tells you what additional packages you need to successfully install `devtools` in any error message you might get when installing this R package.
 
### Running the Makefile

1. First, run `make install_packages` to ensure that you have the necessary packages installed in your R environment.
2. To run all of the analyses, type `make all` in the command line.
3. To run only the translocation analysis, execute `make trans_analysis` in the command line.
4. To run only the viability analysis, execute `make viability_analysis` in the command line.
5. To clean your directory of extraneous files create during these runs, execute `make clean` in the command line.

## Creating journal-formatted PDF documents

1. Render `translocation.qmd` to latex format. 
2. From the doc/manuscript/ subdirectory, run `python convert_qmd_to_pnas_latex.py` in Terminal. 
3. To compile the journal-formatted PDFs, run the following code in Terminal:
 `pdflatex translocation_pnas.tex`
 `bibtex translocation_pnas`
 `pdflatex translocation_pnas.tex`
 `pdflatex translocation_pnas.tex`

 `pdflatex translocation_pnas_SI.tex`
 `bibtex translocation_pnas_SI`
 `pdflatex translocation_pnas_SI.tex`
 `pdflatex translocation_pnas_SI.tex`
4.  If the files fail to compile due to missing latex packages, [install packages](https://en.wikibooks.org/wiki/LaTeX/Installing_Extra_Packages) and rerun the code. In Linux/Ubuntu, packages can be installed using `tlmgr install <package_name>`. 
 
## Contact

Roland Knapp, Research Biologist, University of California Sierra Nevada Aquatic Research Laboratory, Mammoth Lakes, CA 93546 USA; rolandknapp(at)ucsb.edu, <https://mountainlakesresearch.com/roland-knapp/>
