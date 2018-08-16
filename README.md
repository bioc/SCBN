# SCBN
This package take into account the available knowledge of conserved orthologous 
genes and the hypothesis testing framework to get scaling factor and detect differential
expression orthologous genes. The methods on this package are described 
in the article 'A statistical normalization method and differential expression
analysis for RNA-seq data between different species'[1]. 

# Introduction
High-throughput techniques bring novel tools and also statistical challenges to genomic
research. Identifying genes with differential expression between different species is an
effective way to discover evolutionarily conserved transcriptional responses. To remove
systematic variation between different species for a fair comparison, normalization
serves as a crucial pre-processing step that adjusts for the varying sample sequencing
depths and other confounding technical effects.
In this paper, we propose a scale based normalization (SCBN) method by taking
into account the available knowledge of conserved orthologous genes and the hypothesis
testing framework. Considering the different gene lengths and unmapped genes between
species, we formulate the problem from the perspective of hypothesis testing and search
for an optimal scaling factor that minimizes the deviation between the empirical and
nominal type I errors. 

# User's Guide
Please refer to the "SCBN.pdf" vignetee for detailed function 
instructions.

# Reference
1. Zhou Y, Zhu JD, Zhao MT, Zhang BX, Jiang CF, and Yang XY. "Methylation-level 
Inferences and Detection of Differential Methylation with Medip-seq Data". 
Plos one. Accepted.

[![Travis-CI Build Status](https://travis-ci.org/FocusPaka/SIMD.svg?branch=master)](https://travis-ci.org/FocusPaka/SIMD)
