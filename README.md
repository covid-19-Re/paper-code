# Code of Huisman, Scire et al.

Repository containing the code pertaining to the manuscript ["Estimation and worldwide monitoring of the effective reproductive number of SARS-CoV-2"](https://doi.org/10.1101/2020.11.26.20239368) by Jana S. Huisman,  Jérémie Scire, Daniel C. Angst,  Jinzhou Li, Richard A. Neher, Marloes H.  Maathuis, Sebastian Bonhoeffer, Tanja Stadler.

The scripts in this repository assume that the code of the estimation pipeline is available in a folder "covid-19-re-shiny-app" at the same level as this repository. The git-repo for the pipeline can be found [here](https://github.com/covid-19-Re/shiny-dailyRe).  

## Simulation code
- generateSimulations.R contains functions to generate simulations, and estimate Re from simulated observations.
- paperSimulationGrid.R is used to generate Fig. 1 and supplementary Figs. S2-8 (since the simulations take some time to compute, the pre-computed results can be found in the simulations folder; this folder is ~380 Mb).
- interactiveSimulations.R can be used to investigate the outcome of individual simulation scenarios. 

## Estimate stability
- bootstrap_feathers.R is used to generate Fig. 2, it uses a data-file called `all_bootstrap_estimates.csv' which contains the outcome of our Re estimation pipeline on 7 months of daily data files for Switzerland.
- aggregate_bootstrap_estimates.R was included for transparency as it is the file we use to generate `all_bootstrap_estimates.csv' from the individual estimate files.
- ReCountryCluster.R was similarly included for transparency as it is used to estimate Re from each individual incidence data file. This is analogue to the ReCountry.R script used in the [general pipeline](https://github.com/covid-19-Re/shiny-dailyRe), but adapted to accept different incidence data files.

## Dashboard
- Fig. 3 was taken from [our dashboard](https://ibz-shiny.ethz.ch/covid-19-re-international/).

## Empirical analysis code
- empiricalPaper.R is used to generate Table 1 and S2.
- Fig4.R is used to generate Figs. 4 and S12.

## Currently not included
- Status 24.6: The code to generate figures S9-S11 is currently not included since it requires access to Swiss linelist data to run. We are working to provide a processed form of the data and the corresponding scripts for transparency.

