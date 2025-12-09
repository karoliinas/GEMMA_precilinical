# GEMMA_preclinical

This repository contains scripts for preclinical mouse data analysis in GEMMA project.

Initially, Spearman correlations for each variable of interest were calculated, and the correlations were compared between autism and control groups using Fishers R-to-Z transformation (correlation_analysis.R).

Different measurements / data types were integrated using sPLS-DA in mixOmics R-package (mixomics1.R) separately for the mouse strains. Data inputs are normalized matrices with identical samples (in rows) and variables in columns.
