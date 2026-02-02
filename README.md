# CRC_project
This repository contains the four files that were used for performing data analysis in the project "Discrete Transcriptional States Define Biphasic Immune Response and  Dynamic CMS Transitions in Colorectal Cancer"

The file QC_plots_script.R reads in the mouse data, cleans/normalizes it, and produces the mouse figures (figure 3, figure 4, figure 5). Additionally, this file trains the cross-species classifier and produces the plots with the AK-like, AKP-like, and AKPS-like groups (figure 7, figure 8A,8B,8C,8D,8F). Finally, this file contains the code for producing the CMS-alluvial plot (figure 1D). 

The file CIBERSORT_HUMAN.R performs cibersort on the human data and produces the alluvial plots (figure 1C, figure 8E).

The file tcga_gbm_processing.R processes the TCGA-COAD data. Additionally, this code makes the plots with the AK-alterations, AKP-alterations, and AKPS-alterations labels (figure 1A, 1B, 1D).

The file survivalAnalysis.Rmd performs the survival analysis between the AK-like, AKP-like, and AKPS-like groups (figure 8G, 8H). 
