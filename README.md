# Optimal timing and effectiveness on control strategies for COVID19 outbreak in China: a modelling study

Anthony Zhenhuan Zhang (zhan4490@umn.edu)   Eva A. Enns (eenns@umn.edu)


University of Minnesota

## Introduction
This Github repository contains Rscripts to simulation COVID-19 outbreak in major Chinese cities: Wuhan, Chongqing, Beijing, and Shanghai. Below, we first introduce the functionality of each script, and present the model calibration. Lastly, we demonstrate how to run simulation to generate results.

## Scripts
1. functional scripts which simulates SARS-CoV-2 dynamics in Wuhan and other Chinese cities: wuhan_simulation_policy_by_age.R and other_city_simulaiton_policy_by_age.R
2. model inputs generation: model_inputs.R
3. Incremental Mixture Importance Sampling model calibration: model_calibration_ver_3.R
4. generate model outcomes under different timing and duration of control policies: all scripts under  "simulation_by_policy" folder.

## Model Calibration
We calibrated our model using Incremental Mixture Importance Sampling methods, see model_calibration_ver_3.R

## Status quo, and disease burden estimation 
Run decompose_economy_loss_probablistic.R twice. In the first round, we estimate the mean disease burden, in the second round, we generate both the epidemiological and economic outcome


## Results Generation
To generate model results, run scripts under  "simulation_by_policy" folder.

