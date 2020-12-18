# Code of Huisman, Scire et al.

Repository containing the code pertaining to the manuscript ["Estimation and worldwide monitoring of the effective reproductive number of SARS-CoV-2"](https://doi.org/10.1101/2020.11.26.20239368) by Jana S. Huisman,  Jérémie Scire, Daniel C. Angst,  Richard A. Neher, Sebastian Bonhoeffer, Tanja Stadler.

The scripts in this repository assume that the code of the estimation pipeline is available in a folder "covid-19-re-shiny-app" at the same level as this repository. The git-repo for the pipeline can be found [here](https://github.com/covid-19-Re/shiny-dailyRe).  

## Simulation code
- generateSimulations.R contains functions to generate simulations, and estimate Re from simulated observations
- compareSimulations.R contains functions to calculate performance metrics on simulated data
- interactiveSimulations.R is used to generate Fig. 1A
- paperSimulationGrid.R is used to generate Fig. 1B and supplementary Figs.  (since the simulations take some time to compute, the summarised files needed for plotting can be found in the simulations folder)

## Empirical analysis code
- empiricalPaper.R is used to generate Table 1, S1; and Fig S10, S11
- interventionImpact.R is used to generate Fig 3b, 4, S12, and tables S3, S4