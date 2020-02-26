# Tidy Transcriptomics [SOME OTHER TITLE?]

Dr. Maria Doyle (EMAIL) and Dr. Stefano Mangiola (EMAIL)

# Workshop Description

This workshop will present how to perform analysis of RNA sequencing data following the tidiness paradigm. The tidy paradigm follows the representation of the data where each variable is a column and each observation is a row, and the manipulation of data using verbs with an easy-to-understand vocabulary. Most importantly the data structure can be forwarded across manipulation and analysis function keeping its properties (i.e., endomorphic). 

This can be achieved for RNA sequencing data with tidybulk, tidyHeatmap and tidyverse packages. The package tidybulk provides a modular framework and a tidy data structure for bulk transcriptional analyses. tidyHeatmap provides a tidy implementation of the function ComplexHeatmaps::Heatmap. These packages are part of the tidytranscriptomics suite that introduces semantics and properties of a tidy approach to RNA sequencing data.

The topics presented in this workshop will be

- Data exploration
- Data preprocessing and quality control
- Data dimensionality reduction
- Differential transcript abundance analyses 
- Data visualisation

## Pre-requisites

* Basic knowledge of RStudio
* Familiarity with tidyverse syntax

Recommended Background Reading 
https://mblue9.github.io/r-intro-biologists/intro_r_biologists.html

## Workshop Participation

Sudents will be expected to participate in the workshop in a hands-on way, following along with the code provided and performing  exercises.

## _R_ / _Bioconductor_ packages used

* tidyverse
* tidybulk
* tidyHeatmap
* edgeR
* devtools

## Time outline

| Activity                               | Time |
|----------------------------------------|------|
| Data import and quality control        | 30m  |
| Data preprocessing and visualisation   | 30m  |
| Differential transcript abundance      | 30m  |
| Data visualisation                     | 30m  |

# Workshop goals and objectives

In exploring and analysing RNA sequencing data, there are recurrent concept such as filtering, scaling, dimensionality reduction, hypothesis test, clustering and visualisation. These concepts can be intuitively explained to new users, however (i) the heterogeneity of the vocabulary across methodologies/algorithms/packages, (ii) the complexity of data wrangling, and (iii) the coding burden play against an effective learning experience of the statistics and biology underlying an informed data analysis. 

The tidyverse (tidytranscriptomics) approach to (RNA sequencing) data analysis abstracts out the coding related complexity and presents tools that use an intuitive and jargon-free vocabulary, so allowing to focus to the biological and the statistical challenges. Furthermore, it allows the compilation of modular and pipe oriented workflows that allow new user to compare alternative analysis algorithms and workflows.

## Learning goals

* To understand of the recurrent concepts and steps of bulk RNA sequencing data analysis
* To approach data representation and analysis though a tidy paradigm, integrating tidyverse with tidybulk.

## Learning objectives

* To remember the recurrent concept and vocabulary of RNA sequencing data analysis
* To understand the statistical concepts of the most common RNA sequencing data analysis procedures
* To apply the learned concepts to experimental, publicly available data
* To compare for each analysis step different methodologies (e.g., dimensionality reduction with PCA or MDS)
* To evaluate alternative workflows and their consequences/biases on the information that can be learned
* To create plots that summarise the informtion concent of the data and analysis results
