import subprocess
import sys
import re

"""
Convert Quarto to PNAS format for translocation manucript

Generates two files

1. transoclation_pnas.tex : Main text in PNAS format
2. transoclation_pnas_SI.tex : SI in PNAS format
"""

# Files for converting
args = ['translocation.qmd', 'translocation.tex']

### Build the main manuscript ####

# Convert from qmd to latex with pandoc
subprocess.call(['quarto', 'render', args[0], '--to', "latex"])

# Read in the latex manuscript
with open(args[1], 'r') as fin:
    latex_manuscript = fin.readlines()
latex_manuscript = "".join(latex_manuscript)

# Remove all the references to width and height
latex_manuscript = re.sub("includegraphics\[.*\]", 
						  "includegraphics", 
						  latex_manuscript)

with open("PNAS-template-empty.tex") as fin:
	pnas_template = fin.readlines()
pnas_template = "".join(pnas_template)

# Fill in each section of manuscript

## Title ## 

title_start = latex_manuscript.find("\\subsection{Title")
title_end = latex_manuscript[title_start:].find("}") + title_start
title = latex_manuscript[title_start:title_end].replace("\n", " ").split("\\subsection{Title: ")[1]

# Insert title
pnas_ms = pnas_template.replace("\\title{}", "\\title{" + title + "}")

## Abstract ##

abkey1 = "\\hypertarget{abstract}{%\n\\subsubsection{Abstract}\\label{abstract}}"
abkey2 = "\\hypertarget{significance}"
abstract_start = latex_manuscript.find(abkey1)
abstract_end = latex_manuscript.find(abkey2)
abstract = latex_manuscript[abstract_start:abstract_end].split(abkey1)[1].strip("\n")

pnas_ms = pnas_ms.replace("\\begin{abstract}\n\\end{abstract}", 
						  "\\begin{abstract}\n" + abstract + "\n\\end{abstract}")


## Significance statement ##

sigkey = "\\hypertarget{significance}{%\n\\subsubsection{Significance}\\label{significance}}"
intkey1 = "\\hypertarget{introduction}{%\n\\subsubsection{Introduction}\\label{introduction}}"
sig_start = latex_manuscript.find(sigkey)
sig_end = latex_manuscript.find(intkey1)
sig = latex_manuscript[sig_start:sig_end].split(sigkey)[1].strip("\n")

pnas_ms = pnas_ms.replace("\\significancestatement{}", 
						  "\\significancestatement{" + sig + "}")

## Introduction ##

reskey1 = "\\hypertarget{results}{%\n\\subsubsection{Results}\\label{results}}"
intro_start = latex_manuscript.find(intkey1)
intro_end = latex_manuscript.find(reskey1)
intro = latex_manuscript[intro_start:intro_end].split(intkey1)[1].strip("\n")

# Add drop cap
intro = "\\dropcap{" + intro[0] + "}" + intro[1:]

pnas_ms = pnas_ms.replace("% Introduction\n", intro)

## Results ##

diskey1 = "\\hypertarget{discussion}{%\n\\subsubsection{Discussion}\\label{discussion}}"
res_start = latex_manuscript.find(reskey1)
res_end = latex_manuscript.find(diskey1)
results = latex_manuscript[res_start:res_end].split(reskey1)[1].strip("\n")

pnas_ms = pnas_ms.replace("% Results\n",  "\n\\section*{Results}\n\n" + results)

# Specific cleaning of results
targ1 = "\\hypertarget{frog-evolution-in-response-to-bd}{%\n\\paragraph{Frog evolution in response to\nBd}\\label{frog-evolution-in-response-to-bd}}"
pnas_ms = pnas_ms.replace(targ1, "\\subsection*{Frog evolution in response to Bd}")

targ2 = "\\hypertarget{frog-population-recovery}{%\n\\paragraph{Frog population recovery}\\label{frog-population-recovery}}"
pnas_ms = pnas_ms.replace(targ2, "\\subsection*{Frog population recovery}")

targ3 = "\\hypertarget{long-term-population-viability}{%\n\\paragraph{Long-term population\nviability}\\label{long-term-population-viability}}"
pnas_ms = pnas_ms.replace(targ3, "\\subsection*{Long-term population viability}")

## Discussion ##

matmethkey1 = "\\hypertarget{materials-and-methods}{%\n\\subsubsection{Materials and Methods}\\label{materials-and-methods}}"
dis_start = latex_manuscript.find(diskey1)
dis_end = latex_manuscript.find(matmethkey1)
discussion = latex_manuscript[dis_start:dis_end].split(diskey1)[1].strip("\n")

pnas_ms = pnas_ms.replace("% Discussion\n", "\n\\section*{Discussion}\n\n" + discussion)

## Methods ##

ackkey1 = "\\hypertarget{acknowledgments}{%\n\\subsubsection{Acknowledgments}\\label{acknowledgments}}"
methods_start = latex_manuscript.find(matmethkey1)
methods_end = latex_manuscript.find(ackkey1)
methods = latex_manuscript[methods_start:methods_end].split(matmethkey1)[1].strip("\n")

# Specific fixes in the methods

# Update equation to span two columns
methods = methods.replace("\\[\n\\begin{bmatrix}", "\\begin{figure}\\[\n\\begin{bmatrix}")
methods = methods.replace("\\end{bmatrix}(t) \n\\]", "\\end{bmatrix}(t)\\numberthis \\label{eqn:matrix} \n\\]\\end{figure}")

methods = methods.replace("\\paragraph{", "\\section*{")
methods = methods.replace("\\subparagraph{", "\\subsection*{")
methods = methods.replace("Incorporating yearly variability in vital rates", "\\subsubsection*{Incorporating yearly variability in vital rates}")
methods = methods.replace("Estimating model parameters", "\\subsubsection*{Estimating model parameters}")
methods = methods.replace("Model analysis and simulation", "\\subsubsection*{Model analysis and simulation}")

pnas_ms = pnas_ms.replace("\\matmethods{}", "\n\n\\matmethods{\n" + methods + "}\n")

## Acknowledgements ##

figkey = "\\hypertarget{figures}{%\n\\subsubsection{Figures}\\label{figures}}"
ack_start = latex_manuscript.find(ackkey1)
ack_end = latex_manuscript.find(figkey)
acknowledge = latex_manuscript[ack_start:ack_end].split(ackkey1)[1].strip("\n")

pnas_ms = pnas_ms.replace("\\acknow{}", "\\acknow{" + acknowledge + "}")

## Figures ##

suppmatkey = "\\hypertarget{supporting-information}{%\n\\subsection{Supporting Information}\\label{supporting-information}}"
fig_start = latex_manuscript.find(figkey)
fig_end = latex_manuscript.find(suppmatkey)
figures = latex_manuscript[fig_start:fig_end].split(figkey)[1].strip("\n")#.strip("\\newpage")#.strip("hfill\\break")
pnas_ms = pnas_ms.replace("% Figures\n", "\\clearpage\n" + figures)

pnas_ms = pnas_ms.replace("\\begin{figure", "\\begin{figure*")
pnas_ms = pnas_ms.replace("\\end{figure", "\\end{figure*")
pnas_ms = pnas_ms.replace("\\includegraphics", "\\includegraphics[width=0.8\\textwidth]")


# Clean up references to SI material

supp_map = {"Figure~\\ref{fig-selectionresults} A: SI": "Figure S1A",
			"Figure~\\ref{fig-selectionresults} B, C: SI": "Figure S1B,C",
			"Figure~\\ref{fig-selectionresults} B, C: SI": "Figure S1B,C",
			"Figure~\\ref{fig-synteny-plot} SI": "Figure S2",
			"Figure~\\ref{fig-yosemap} SI": "Figure S3",
			"Figure~\\ref{fig-yosemap}\nSI": "Figure S3",
			"Figure~\\ref{fig-transsurvival-postdens} SI": "Figure S4",
			"Figure~\\ref{fig-bdload-beforeafter} SI": "Figure S5",
			"Figure~\\ref{fig-bdload-beforeafter}\n SI": "Figure S5",
			"Figure~\\ref{fig-bdload-beforeafter}\nSI": "Figure S5",
			"Figure~\\ref{fig-bdload-beforeafter}": "Figure S5",
			"Figure~\\ref{fig-survival-postdens} SI": "Figure S6",
			"Figure~\\ref{fig-compare_surv_probs} : SI": "Figure S7",
			"Figure~\\ref{fig-compare_surv_probs} SI": "Figure S7",
			"Table~\\ref{tbl-survival-earlylate} : SI": "Table S1",
			"Table~\\ref{tbl-survival-earlylate} SI": "Table S1",
			"Table~\\ref{tbl-survival-earlylate}": "Table S1",
			"Table~\\ref{tbl-param_values} : SI": "Table S2",
			"Table~\\ref{tbl-param_values}\n: SI": "Table S2",
			"Table~\\ref{tbl-param_values} SI": "Table S2"
			}

for key, value in supp_map.items():
	pnas_ms = pnas_ms.replace(key, value)

# Save result
with open("translocation_pnas.tex", "w") as fout:
	fout.writelines(pnas_ms)

#### Build the Supporting Information #####

# Load in supportin template
with open("PNAS-template-supp-info-empty.tex", 'r') as fin:
	pnas_supp = fin.readlines()
pnas_supp = "".join(pnas_supp)

## Main text of SI ##

figkey2 = "\\hypertarget{figures-1}{%\n\\subsubsection{Figures}\\label{figures-1}}"
supp_start = latex_manuscript.find(suppmatkey)
supp_end = latex_manuscript.find(figkey2)
supp_info = latex_manuscript[supp_start:supp_end].split(suppmatkey)[1].strip("\n")

supp_info = supp_info.replace("\\paragraph", "\\subsection")
supp_info = supp_info.replace("\\subsubsection", "\\section")

pnas_supp = pnas_supp.replace("\\SItext", "\\SItext\n" + supp_info)

## SI figures ##

datakey = "\\hypertarget{datasets}{%\n\\subsubsection{Datasets}\\label{datasets}}"
fig_start = latex_manuscript.find(figkey2)
fig_end = latex_manuscript.find(datakey)
figures_supp = latex_manuscript[fig_start:fig_end].split(figkey2)[1].strip("\n")

figures_supp = figures_supp.replace("\\end{figure}", "\\end{figure}\\clearpage")
figures_supp = figures_supp.replace("\\includegraphics", "\\includegraphics[width=0.8\\textwidth]")

pnas_supp = pnas_supp.replace("% Figures", figures_supp)

# Account for main text reference in SI

main_map = {
			"Figure~\\ref{fig-spline-manhattan}": "Fig. 2",
			"Figure~\\ref{fig-allelefrequencies}": "Fig. 1",
			"Figure~\\ref{fig-cond-effects}": "Fig. 4",
			"Figure~\\ref{fig-translocation-survival}": "Fig. 3",
			} 

for key, value in main_map.items():
	pnas_supp = pnas_supp.replace(key, value)


# Remove verbatim
verb_start = latex_manuscript.find("\\begin{verbatim}")
verb_end = latex_manuscript.find("\\end{verbatim}")
verb = latex_manuscript[verb_start:(verb_end + len("\\end{verbatim}"))]

pnas_supp = pnas_supp.replace(verb, "")

pnas_supp = pnas_supp.replace("\\title{}", "\\title{" + title + "}")

with open("translocation_pnas_SI.tex", "w") as fout:
	fout.writelines(pnas_supp)

# Remove transcloation.tex
# subprocess.call(['rm', args[1]])


