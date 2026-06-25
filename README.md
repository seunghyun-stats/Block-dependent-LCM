# Block-dependent-LCM

This repository contains the R codes for the manuscript [Lee, S. and Gu, Y. (2026), Learning High-Dimensional Block-Dependent Latent
Class Models](https://arxiv.org/abs/).

### For simulations:
The folder `Simulations-Main` contains the codes for the proposed three-step method (see Section 3 in the manuscript). `simulation_main.R` contains both the code for the method and main simulations. `simulation_CV.R` contains the code that is used to select the graphical lasso tuning parameter. The `Comparion` folder contains additional codes to implement and compare the performance with the joint graphical lasso method.

### For real data analysis:
To replicate data analysis, go to the folder `Real data`. Each R script corresponds to the main code for the analysis in Section 5 of the manuscript. Guidelines for downloading the data are provided as comments.
