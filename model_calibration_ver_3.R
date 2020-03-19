# setwd("~/Desktop/Spring 2020/coronavirus/codes")
setwd("~/Dropbox/codes")
rm(list = ls())

# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)
library(IMIS)
library(matrixStats) # package used for sumamry statistics
library(lubridate)
require("sfsmisc")
require("deSolve")
# visualization


##### target: infected individuals/prevalence on day 56
##### 030920: change target to confirmed cases by 9th March

########### Model inputs
# initial population N
# Wuhan: 
N_wh = 11e+06 + 8e+06
# Chongqing: 
N_cq = 30.48e+06   
# Beijing: 
N_bj = 21.54e+06
# Shanghai: 
N_sh = 24.24e+06
# calculate age distribution for other cities, 
# obtained from https://www.citypopulation.de/php/china-hubei-admin_c.php?adm1id=4201
f_cq = c(3106599 + 4042103, 3715039 + 4308372 + 4927387 + 3722276, 2910105 + 1540808 + 573481)
f_cq = f_cq/sum(f_cq)
f_bj = c(1187149 + 1571366, 5012042 + 3544940 + 3285159 + 2551604, 1267620 + 890379 + 302109)
f_bj = f_bj/sum(f_bj)
f_sh = c(1426078 + 1677953, 5191551 + 4049413 + 3678414 + 3526132, 1803695 + 1078156 + 587804)
f_sh = f_sh/sum(f_sh)


exportion_perc_cq = c(0.0115, 0.0127, 0.0044) # 1 normal, 2 chunyun 3.after quarantine
# exportion_perc_bj = c(0.0175, 0.0086)
exportion_perc_bj = c(0.0086, 0.0086, 0.0024)
# exportion_perc_sh = c(0.0125, 0.0066)
exportion_perc_sh = c(0.0066, 0.0066, 0.0029)
###### input for models: 
# time = 365 # simulate up to a year
### to save time, simulate up to April 9th (130 times)
time = 130

# R0_0 = 2.68 (mean)
z_0 = 43 * 2 # zoonotic force of infection, 86 cases per day
D_E = 6 # average 6 days of incubation period from Lancet
D_I = 8.4  - D_E # serial interval minus the mean latent period
E_I_beta_ratio = 0.005
# mu = c(1e-06, 1.25e-05, 2.6e-04) # 020920: mortality rate: need to calibrate with China CDC data

##### non-pharmaceutical interventions  ####


simulation_start = as.Date("12/01/19", "%m/%d/%y")
simulation_end = simulation_start + time 
chunyun_start = as.Date("01/10/20", "%m/%d/%y")
chunyun_end = chunyun_start + 40 - 1

quarantine_start = as.Date("01/23/20", "%m/%d/%y")
# quarantine_end = as.Date("03/01/20", "%m/%d/%y")
quarantine_end = simulation_end
social_distancing_start = quarantine_start
# social_distancing_end = quarantine_end
nage = length(f_cq)
SD_end_dates = structure(rep(NA_real_, length(nage)), class="Date")
#### for other cities
ths_start = quarantine_start ###ths: travel history screening
ths_end = quarantine_end ### ths: travel history screening
SD_end_dates[1] = as.Date("03/31/20", "%m/%d/%y")
SD_end_dates[2] = as.Date("03/31/20", "%m/%d/%y")
SD_end_dates[3] = as.Date("03/31/20", "%m/%d/%y")

# social_distancing_end = delta_t + base_line_social_distancing_end
time_knot_vec = list('simulation_start' = simulation_start,'simulation_end' = simulation_end, 
                     'chunyun_start' = chunyun_start, 'chunyun_end' = chunyun_end, 
                     'quarantine_start' = quarantine_start, 'quarantine_end' = quarantine_end,  
                     'social_distancing_start' = social_distancing_start, 'SD_end_dates' = SD_end_dates)
######### 022720: re-calibrate
reduction = 0.5
contact_mat = rbind(c(21.09312, 20.44581, 15.02198),# normal
                    c(21.09312, 20.44581, 15.02198)* 1.2, # cny
                    # c(3, 3, 3))  #  during outbreak
                    c(21.09312, 20.44581, 15.02198)*(1-reduction))  #  during outbreak

contact_ratio = rbind(c(0.5928983 - (1-0.5928983)/3, 0.8353008 - (1- 0.8353008)/3, 0.29090466 - (1 - 0.29090466)/3), # normal
                      c(0, 0, 0)) # during cny

Wuhan_cases = 49948 # on Mar 9
Wuhan_death = 2404 # on Mar 9
#### 0 - 19: 1/1023*2404   20 - 60: 193/1023*2404 60+: 829/1023*2404   ###
confirmed_perc = c(0.021, 0.667, 0.312) 
#### 0 - 19: 1/1023*2349   20 - 60: 193/1023*2349 60+: 829/1023*2349   ###
death_perc = c(1/1023, 193/1023, 829/1023)
# symptomatic_ratio = c(2/6, (301-2)/(301-2+318-4), (301-2)/(301-2+318-4)) # from Diamond Princess paper
### maybe compute symptomatic ratio of child on our own
symptomatic_ratio = c(2/60, (301-2 - 200)/(301-2+318-4 -200- 265), (200)/(465))
days_symptom_to_death = c(20, 20, 11.3) # from https://www.ncbi.nlm.nih.gov/pubmed/31994742
## We need to match death cases 20 days ago for Youth, Adults, and ???
## to complicated, time-consuming...
# elderly_diagnose_prob = 1



###### source wuhan model and other city model

source('wuhan_simulation_policy_by_age.R')






####################################################################
######  Specify calibration parameters  ######
####################################################################

# names and number of input parameters to be calibrated
v_param_names <- c("S_A_ratio_youth", "S_A_ratio_adult", "S_A_ratio_elderly",
                   "mortality_youth", "mortality_adult", "mortality_elderly", "beta_IS")
n_param <- length(v_param_names)
##### changable values!

# symptomatic_ratio = c(2/60, (301-2 - 200)/(301-2+318-4 -200- 265), (200)/(465))
# mu_1 = c(6.36e-05, 5.56e-03, 5.50e-02)## calibrated assuming social distancing with 50% reduction
# mu_2 = mu_1/7
param_1 <- c(S_A_ratio_youth = 2, S_A_ratio_adult = 99,     S_A_ratio_elderly = 200,
             mortality_youth = 0, mortality_adult = 0.00076, mortality_elderly = 0.038, beta_IS = 0.05)
# param_2 <- c(S_A_ratio_youth = 4, S_A_ratio_adult = 49,     S_A_ratio_elderly = 265,
#              mortality_youth = 0.0014, mortality_adult = 0.034, mortality_elderly = 0.1, beta_IS = 0.07)
###### debug
param_2 <- c(S_A_ratio_youth = 4, S_A_ratio_adult = 49,     S_A_ratio_elderly = 265,
             mortality_youth = 0.00014, mortality_adult = 0.011, mortality_elderly = 0.082, beta_IS = 0.07)
##### S_A_ratio prior from Diamond Princess ship study
##### mortality rate range from https://www.medrxiv.org/content/10.1101/2020.03.04.20031104v1.full.pdf
# number of calibration targets
# v_target_names <- c("Prevalence on 25 Jan", "Exported to CQ", "Exported to BJ", "Exported to SH")
# n_target <- length(v_target_names)
# 
# v_target_names <- c("Prevalence on 25 Jan")
# n_target <- length(v_target_names)
# target_vec = c(75815, 37304, 130330) # mean, min, max
# CQ_vec = c(461, 227, 805)
# BJ_vec = c(98, 49, 168)
# SH_vec = c(80, 40, 139)
v_target_names <- c("Elderly confirmed cases", "Youth death", "Adult death", "Elderly death") # all on Mar 9 2020
n_target <- length(v_target_names)



### Calibration Functions

# Write function to sample from prior
sample_prior <- function(n_samp){
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  
  for (i in 1:3){
    # m_param_samp[, i] <- qunif(m_lhs_unit[,i],
    #                            min = lb[i],
    #                            max = ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
                               shape1 = param_1[i],
                               shape2 = param_2[i])
  }
  for (i in 4:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = param_1[i],
                               max = param_2[i])
   
  }
  # return(as.numeric(m_param_samp))
  return((m_param_samp))
}



###  PRIOR  ### 
# Write functions to evaluate log-prior and prior

# function that calculates the log-prior
calc_log_prior <- function(v_params){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  
  for (i in 1:3){
    lprior <- lprior + dbeta(v_params[, i],
                             shape1 = param_1[i],
                             shape2 = param_2[i],
                             log = T)
  }
  for (i in 4: n_param){
  lprior <- lprior + dunif(v_params[, i],
                           min = param_1[i],
                           max = param_2[i], 
                           log = T)
  }
  
  return(lprior)
}

# function that calculates the (non-log) prior
calc_prior <- function(v_params) { 
  exp(calc_log_prior(v_params)) 
}
calc_log_prior(v_params = sample_prior(10))
calc_prior(v_params =  sample_prior(10))
###  LIKELIHOOD  ###
# Write functions to evaluate log-likelihood and likelihood

# function to calculate the log-likelihood
calc_log_lik <- function(v_params){
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_samp <- nrow(v_params)
  # v_llik <-  rep(0, n_samp, n) # only 1 target
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parametr set "v_params" ###
      # model_res <- run_sick_sicker_markov(v_params[j, ])
      symptomatic_ratio = as.numeric(v_params[j, 1:3])
      mu_1 = as.numeric(v_params[j, 4:6])
      beta = as.numeric(v_params[j, 7])
      wuhan_output <- wuhan_simulation(N_wh, time, beta, z_0, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu_1, symptomatic_ratio)
      # sum(wuhan_output$I[56,]) # simulation start  + 55 days in Jan 25th 2020
      ###  Calculate log-likelihood of model outputs to targets  ###
      v_llik[j, 1] <- sum(dpois(x = (round(wuhan_output$IS[100, 3] + wuhan_output$RS[100, 3])),
                                lambda =  Wuhan_cases * confirmed_perc[3], # confirmed elderly cases,
                                log = T))
      #### calibrate mortality
      # mu_1 = c(6.36e-05, 5.56e-03, 5.50e-02) ### calibrated assuming social distancing with 50% reduction
      # mu_2 = mu_1/7 ## from China CDC case fatality ratio  
       #  wuhan_output$D[80,1], wuhan_output$D[80,2], wuhan_output$D[89,3])
      v_llik[j, 2] <- sum(dpois(x = (round(wuhan_output$D[80,1])),
                                lambda =  Wuhan_death * death_perc[1], # Youth death,
                                log = T))
      v_llik[j, 3] <- sum(dpois(x = (round(wuhan_output$D[80,2])),
                                lambda =  Wuhan_death * death_perc[2], # Adult death,
                                log = T))
      v_llik[j, 4] <- sum(dpois(x = (round(wuhan_output$D[89,3])),
                                lambda =  Wuhan_death * death_perc[3], # Elderly death,
                                log = T))
      
      
      # OVERALL 
      # llik_overall[j] <- sum(v_llik[j])
      llik_overall[j] <- sum(v_llik[j,])
    }, error = function(e) NA) 
    if(is.na(jj)) { llik_overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik_overall)
}

# function to calculate the (non-log) likelihood
calc_likelihood <- function(v_params){ 
  exp(calc_log_lik(v_params)) 
}
calc_log_lik(v_params = sample_prior(10))
calc_likelihood(v_params = sample_prior(10))


###  POSTERIOR  ###
# Write functions to evaluate log-posterior and posterior

# function that calculates the log-posterior
calc_log_post <- function(v_params) { 
  lpost <- calc_log_prior(v_params) + calc_log_lik(v_params)
  return(lpost) 
}

calc_log_post(v_params = sample_prior(10))

# function that calculates the (non-log) posterior
calc_post <- function(v_params) { 
  exp(calc_log_post(v_params)) 
}
# calc_post(v_params = v_params_test)
calc_post(v_params = sample_prior(10))


####################################################################
######  Calibrate!  ######
####################################################################
# record start time of calibration
t_init <- Sys.time()
###  Bayesian calibration using IMIS  ###
# define three functions needed by IMIS: prior(x), likelihood(x), sample.prior(n)
prior <- calc_prior
likelihood <- calc_likelihood
sample.prior <- sample_prior

# run IMIS
# Specify seed (for reproducible sequence of random numbers)
set.seed(02298)

# number of random samples
n_resamp = 1000
fit_imis <- IMIS(B = 2000, # the incremental sample size at each iteration of IMIS
                 B.re = n_resamp, # the desired posterior sample size
                 number_k = 30, # the maximum number of iterations in IMIS
                 D = 0) 

# obtain draws from posterior
m_calib_res <- fit_imis$resample

# Calculate log-likelihood (overall fit) and posterior probability of each sample
m_calib_res <- cbind(m_calib_res,
                      "Overall_fit" = calc_log_lik(m_calib_res),
                     "Posterior_prob" = calc_post(m_calib_res))

# normalize posterior probability
m_calib_res[,"Posterior_prob"] <- m_calib_res[,"Posterior_prob"]/sum(m_calib_res[,"Posterior_prob"])
m_calib_res<- cbind(m_calib_res,
                    "c_posterior_prob" = cumsum(m_calib_res[,"Posterior_prob"]))

# Calculate computation time
comp_time <- Sys.time() - t_init

save(m_calib_res, file = 'm_calib_res.rda')














####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

load('m_calib_res.rda')
v_param_names <- c("S_A_ratio_youth", "S_A_ratio_adult", "S_A_ratio_elderly",
                   "mortality_youth", "mortality_adult", "mortality_elderly", "beta_IS")
n_param <- length(v_param_names)
# Compute posterior mean
v_calib_post_mean <- colMeans(m_calib_res[,v_param_names])
v_calib_post_mean
# Compute posterior median and 95% credible interval
m_calib_res_95cr <- colQuantiles(m_calib_res[,v_param_names], probs = c(0.025, 0.5, 0.975))
m_calib_res_95cr
m_calib_res_CI_and_mean <- m_calib_res_95cr
m_calib_res_CI_and_mean[, 2] <- v_calib_post_mean
round(m_calib_res_CI_and_mean*100,2)


# Compute maximum-a-posteriori (MAP) parameter set
# v_calib_map <- m_calib_res[which.max(m_calib_res[,"Posterior_prob"]),]
# 
# 
# ### Plot model-predicted output at best set vs targets ###
# v_out_best <- run_sick_sicker_markov(v_calib_map[v_param_names])




########################################################################
########## Now, calibrate ratio between symptomatic and clinical confirmation etc.
########################################################################

# symptomatic_ratio = as.numeric(v_calib_post_mean[1:3])
# mu_1 = as.numeric(v_calib_post_mean[4:6])
# beta = as.numeric(v_calib_post_mean[7])
# wuhan_output <- wuhan_simulation(N_wh, time, beta, z_0, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu_1, symptomatic_ratio)
# 
# ratio = (Wuhan_cases * confirmed_perc)/(wuhan_output$IS[100, ] + wuhan_output$RS[100, ])
n_sample = 1000
ratio_mat = matrix(0, n_sample, 3)
symptomatic_mat = matrix(0, n_sample, 3)
D_mat = matrix(0, n_sample, 3)
for (ind in 1:n_sample){ 
  # param_ind = which.min(abs(runif(1) - m_calib_res[,"c_posterior_prob"]))
  param_ind = ind
  symptomatic_ratio = as.numeric(m_calib_res[param_ind, 1:3])
  mu_1 = as.numeric(m_calib_res[param_ind, 4:6])
  beta = as.numeric(m_calib_res[param_ind, 7])
  wuhan_output <- wuhan_simulation(N_wh, time, beta, z_0, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu_1, symptomatic_ratio)
  
  ratio_mat[ind, ] = (Wuhan_cases * confirmed_perc)/(wuhan_output$IS[100, ] + wuhan_output$RS[100, ])
  symptomatic_mat[ind, ] = (wuhan_output$IS[100, ] + wuhan_output$RS[100, ])
  D_mat[ind, ] = c(wuhan_output$D[80,1], wuhan_output$D[80,2], wuhan_output$D[89,3])
}

symptomatic_mat_95cr <- colQuantiles(symptomatic_mat, probs = c(0.025, 0.975))
ratio_mat_95cr <- colQuantiles(ratio_mat, probs = c(0.025, 0.975))
D_mat_95cr  <- colQuantiles(D_mat, probs = c(0.025, 0.975))
symptomatic_mat_mean <- colMeans(symptomatic_mat)
ratio_mat_mean <- colMeans(ratio_mat)
D_mat_mean <- colMeans(D_mat)

round(D_mat_mean/sum(D_mat_mean) * 100, 2)
round(D_mat_95cr[,1]/sum(D_mat_95cr[,1])* 100, 2)
round(D_mat_95cr[,2]/sum(D_mat_95cr[,2])* 100, 2)



