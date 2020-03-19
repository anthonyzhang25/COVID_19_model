##### we evaluate health outcome from various potential intervetions in this script
# setwd("~/Desktop/Spring 2020/coronavirus/codes")
setwd("~/Dropbox/codes")
rm(list = ls())
library(lubridate)
require("sfsmisc")
require("deSolve")
library(matrixStats)



###### source wuhan model and other city model
source('wuhan_simulation_policy_by_age.R')
source('other_city_simulation_policy_by_age.R')
load('m_calib_res.rda')
load('model_inputs.rda')




############### Sensitivity Analysis on different combinations of control policies ###########
################################################################################################
################################      Change in parameters ################################
################################################################################################


################ sens_param_1: change the E_I_beta_ratio
################ sens_param_3: contact rate during chunyun/CNY
# contact_mat = rbind(c(20, 20, 20),# normal
#                     c(20, 20, 20) * 1.2, # cny
#                     # c(3, 3, 3))  #  during outbreak
#                     c(10, 10, 10))  #  during outbreak
id = 4
E_I_beta_ratio = 0.005
milestone = as.Date("03/31/20", "%m/%d/%y")
# base_line_social_distancing_end = as.Date("02/29/20", "%m/%d/%y")
start = 0.25
incre = 0.05
end = 0.75
param_range = seq(start, end, by = incre)

### debug 
# param_range = 25:30
nage = length(f_cq)
CI_and_mean = 3
date = structure(rep(NA_real_, length(param_range)), class="Date")
cum_exportion = array(0, c(length(param_range), CI_and_mean))
cum_infected = array(0, c(length(param_range), 4, CI_and_mean)) # 4 cities: Wuhan, CQ, BJ, SH
cum_death = array(0, c(length(param_range), 4, CI_and_mean)) # 4 cities: Wuhan, CQ, BJ, SH
disease_burden = array(0, c(length(param_range), 4, CI_and_mean)) # 4 cities: Wuhan, CQ, BJ, SH
economy_loss = array(0, c(length(param_range), 4, CI_and_mean)) # 4 cities: Wuhan, CQ, BJ, SH
# THS_effectiveness = array(0, c(length(param_range), 3, CI_and_mean)) # only for CQ, BJ, SH
screening_cost = array(0, c(length(param_range), 3, CI_and_mean))
SD_end_dates = structure(rep(NA_real_, length(nage)), class="Date")
n_sample = 1000
temp_disease_burden = matrix(0, 4, n_sample)
temp_cum_infected = matrix(0, 4, n_sample)
temp_cum_death = matrix(0, 4, n_sample)
temp_economy_loss = matrix(0, 4, n_sample)
temp_screening_cost = matrix(0, 3, n_sample)
# sample_param_indexes = which.min(abs(runif(n_sample) - m_calib_res[,"c_posterior_prob"]))
t_init = Sys.time()
for (reduction in param_range){
  # social_distancing_end = social_distancing_end_vec[ind]
  reduction_ind = as.integer(as.character((reduction - start)/incre)) + 1
  cat("reduction is ",reduction, "reduction_ind is ", reduction_ind, "\n")
  contact_mat = rbind(c(21.09312, 20.44581, 15.02198),# normal
                      c(21.09312, 20.44581, 15.02198) * 1.2, # cny
                      # c(3, 3, 3))  #  during outbreak
                      c(21.09312, 20.44581, 15.02198) * reduction)  #  during outbreak
  base_line_quarantine_startdate = as.Date("01/23/20", "%m/%d/%y")
  for (ind in 1:n_sample){
    # if (ind %% 50 == 0){
    #   cat('index is ', ind, "\n")
    # }
    # param_ind = which.min(abs(runif(1) - m_calib_res[,"c_posterior_prob"]))
    param_ind = ind
    symptomatic_ratio = as.numeric(m_calib_res[param_ind, 1:3])
    mu_1 = as.numeric(m_calib_res[param_ind, 4:6])
    beta = as.numeric(m_calib_res[param_ind, 7])
    # ratio = ratio_range[beta_ind]
    
    
    quarantine_start = base_line_quarantine_startdate
    quarantine_end = simulation_end
    ### we assume social_distancing start as the same time as 
    social_distancing_start = quarantine_start
    social_distancing_end = quarantine_end
    SD_end_dates[1] = as.Date("03/31/20", "%m/%d/%y")
    SD_end_dates[2] = as.Date("03/31/20", "%m/%d/%y")
    SD_end_dates[3] = as.Date("03/31/20", "%m/%d/%y")
    
    # social_distancing_end = delta_t + base_line_social_distancing_end
    time_knot_vec = list('simulation_start' = simulation_start,'simulation_end' = simulation_end, 
                         'chunyun_start' = chunyun_start, 'chunyun_end' = chunyun_end, 
                         'quarantine_start' = quarantine_start, 'quarantine_end' = quarantine_end,  
                         'social_distancing_start' = social_distancing_start, 'SD_end_dates' = SD_end_dates)
    
    mu = mu_1
    WH_output <- wuhan_simulation(N_wh, time, beta, z_0, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu, symptomatic_ratio)
    # ratio = ratio_range[2]
    #### for other cities
    exported = WH_output$exported
    ths_efficiency = 1
    ths_window = 14
    ths_start = quarantine_start ###ths: travel history screening
    ths_end = quarantine_end ### ths: travel history screening
    # social_distancing_start = as.Date("01/23/20", "%m/%d/%y")
    social_distancing_start = ths_start
    SD_end_dates[1] = as.Date("02/29/20", "%m/%d/%y")
    SD_end_dates[2] = as.Date("02/29/20", "%m/%d/%y")
    SD_end_dates[3] = as.Date("03/31/20", "%m/%d/%y")
    time_knot_vec = list('simulation_start' = simulation_start,'simulation_end' = simulation_end, 
                         'chunyun_start' = chunyun_start, 'chunyun_end' = chunyun_end, 
                         'ths_start' = ths_start, 'ths_end' = ths_end,  
                         'social_distancing_start' = social_distancing_start, 'SD_end_dates' = SD_end_dates) 
    
    mu = mu_1/7
    CQ_output <- other_city_simulation(N_cq, f_cq, time, beta, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu, exported, exportion_perc_cq, ths_efficiency, ths_window, symptomatic_ratio)
    BJ_output <- other_city_simulation(N_bj, f_bj, time, beta, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu, exported, exportion_perc_bj, ths_efficiency, ths_window, symptomatic_ratio)
    SH_output <- other_city_simulation(N_sh, f_sh, time, beta, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu, exported, exportion_perc_sh, ths_efficiency, ths_window, symptomatic_ratio)
    diff_t = as.integer(milestone - simulation_start)
    temp_disease_burden[1, ind] = sum((WH_output$IS[diff_t + 1, ] + WH_output$RS[diff_t + 1, ]) * ratio * (cost_treatment + cost_survived) +
                                        (WH_output$D[diff_t + 1, ]) * cost_death) /1e+06
    temp_disease_burden[2, ind] = sum((CQ_output$IS[diff_t + 1, ] + CQ_output$RS[diff_t + 1, ]) * ratio * (cost_treatment + cost_survived) + 
                                        (CQ_output$D[diff_t + 1, ]) * cost_death) /1e+06
    temp_disease_burden[3, ind] = sum((BJ_output$IS[diff_t + 1, ] + BJ_output$RS[diff_t + 1, ]) * ratio * (cost_treatment + cost_survived) + 
                                        (BJ_output$D[diff_t + 1, ]) * cost_death) /1e+06
    temp_disease_burden[4, ind] = sum((SH_output$IS[diff_t + 1, ] + SH_output$RS[diff_t + 1, ]) * ratio * (cost_treatment + cost_survived) + 
                                        (SH_output$D[diff_t + 1, ]) * cost_death) /1e+06 
    temp_cum_infected[1, ind] = sum(WH_output$incidence[1:diff_t + 1])/WH_output$N_t[diff_t + 1] * 10000
    temp_cum_infected[2, ind] = sum(CQ_output$incidence[1:diff_t + 1])/CQ_output$N_t[diff_t + 1] * 10000
    temp_cum_infected[3, ind] = sum(BJ_output$incidence[1:diff_t + 1])/BJ_output$N_t[diff_t + 1] * 10000
    temp_cum_infected[4, ind] = sum(SH_output$incidence[1:diff_t + 1])/SH_output$N_t[diff_t + 1] * 10000
    
    temp_cum_death[1, ind] = sum(WH_output$D[diff_t + 1,]) 
    temp_cum_death[2, ind] = sum(CQ_output$D[diff_t + 1,]) 
    temp_cum_death[3, ind] = sum(BJ_output$D[diff_t + 1,]) 
    temp_cum_death[4, ind] = sum(SH_output$D[diff_t + 1,]) 
    
    temp_screening_cost[1, ind] = CQ_output$screening_counts * 16.438/1e+06
    temp_screening_cost[2, ind] = BJ_output$screening_counts * 16.438/1e+06
    temp_screening_cost[3, ind] = SH_output$screening_counts * 16.438/1e+06
    
    
    temp_economy_loss[, ind] = temp_economy_loss[ ,ind] + 
      GDP_loss_matrix[1,] * (milestone - quarantine_start + 1) + 
      GDP_loss_matrix[2, ] * (SD_end_dates[2] - social_distancing_start + 1) + 
      temp_disease_burden[ ,ind] + c(0, temp_screening_cost[ ,ind])
    temp_economy_loss[ , ind] = round(temp_economy_loss[ ,ind]/1000,2)
    
  }
  # cum_exportion[reduction_ind, beta_ind] = sum(WH_output$exported[1:diff_t + 1,,])
  disease_burden[reduction_ind, ,] = rowQuantiles(temp_disease_burden, probs = c(0.025, 0.5, 0.975))
  cum_infected[reduction_ind, ,] = rowQuantiles(temp_cum_infected, probs = c(0.025, 0.5, 0.975))
  cum_death[reduction_ind, ,] = rowQuantiles(temp_cum_death, probs = c(0.025, 0.5, 0.975))
  screening_cost[reduction_ind, ,] = rowQuantiles(temp_screening_cost, probs = c(0.025, 0.5, 0.975))
  economy_loss[reduction_ind, ,] = rowQuantiles(temp_economy_loss, probs = c(0.025, 0.5, 0.975))
  ### change meadian to mean
  disease_burden[reduction_ind, ,2] = rowMeans(temp_disease_burden)
  cum_infected[reduction_ind, ,2] = rowMeans(temp_cum_infected)
  cum_death[reduction_ind, ,2] = rowMeans(temp_cum_death)
  screening_cost[reduction_ind, ,2] = rowMeans(temp_screening_cost)
  economy_loss[reduction_ind, ,2] = rowMeans(temp_economy_loss)
  
  # date[reduction_ind, 1] = as.Date(base_line_quarantine_startdate + delta_t)
  
  
}

comp_time = Sys.time() - t_init
comp_time


save(economy_loss, disease_burden, cum_infected, cum_death, param_range, screening_cost, id, file = "vary_contact_reduction.rda")





############### Sensitivity Analysis on different combinations of control policies





  
