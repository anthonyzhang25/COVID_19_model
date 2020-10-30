
cost_treatment = 4124.64
cost_survived = c(131.94, 1317.08, 154.22) # from Appendix, and script calculate_GDP_loss.R
cost_death = c(689680.60, 380892.36, 97043.36) # from Appendix, and script calculate_GDP_loss.R
WTP_ratio = 3 # 3 times China's GDP per capita
cost_survived = cost_survived * WTP_ratio
cost_death = cost_death * WTP_ratio


# total_loss = c(37.48, 1.841, 2.761, 2.978) * 1000
### 091920: updated total economy loss 
total_loss = c(23.46, 8.85, 13.27, 14.32) * 1000 # in million
# disease_burden_mean = c(482.30285, 174.51075,  70.25288,  54.15841)
# disease_burden_mean = c(782.68505, 27.70816, 25.01992,19.08906)
# if (WTP_ratio == 3) disease_burden_mean = c(1920.32486,   50.30923,   35.73182,   25.72194)
# if (WTP_ratio == 1) disease_burden_mean = c(802.17082, 30.40620,  22.07904,  15.90081)

control_policy_cost = total_loss
daily_cost_policy = control_policy_cost/ c(69,38,38,38)
daily_cost_SD = c(daily_cost_policy[2] * 1.2 * 0.3511, daily_cost_policy[2:4])
daily_cost_quarantine = daily_cost_policy - daily_cost_SD
GDP_loss_matrix = rbind(c(daily_cost_policy[1], 0, 0, 0),
                        c(0, daily_cost_SD[2:4]))
# initial population N
# Wuhan: 
N_wh = 11e+06 + 8e+06
# Chongqing: 
N_cq = 30.48e+06   
# Beijing: 
N_bj = 21.54e+06
# Shanghai: 
N_sh = 24.24e+06
N = c(N_wh, N_cq, N_bj, N_sh)
# calculate age distribution for other cities, 
# obtained from https://www.citypopulation.de/php/china-hubei-admin_c.php?adm1id=4201
f_cq = c(3106599 + 4042103, 3715039 + 4308372 + 4927387 + 3722276, 2910105 + 1540808 + 573481)
f_cq = f_cq/sum(f_cq)
f_bj = c(1187149 + 1571366, 5012042 + 3544940 + 3285159 + 2551604, 1267620 + 890379 + 302109)
f_bj = f_bj/sum(f_bj)
f_sh = c(1426078 + 1677953, 5191551 + 4049413 + 3678414 + 3526132, 1803695 + 1078156 + 587804)
f_sh = f_sh/sum(f_sh)

#### from Baidu qianxi:
# exported from Wuhan to other cities from 01-10 to 01-24 https://qianxi.baidu.com/
# Chongqing: 1.27%
# Beijing: 0.86%
# Shanghai: 0.66%

# From December to Jan 10th we use Jan 1st to Jan 9th and average
# Chongqing: 1.1566%   averaged from: 0.85% 1.15% 1.03% 1.04% 1.14% 1.23% 1.24% 1.27% 1.46% 

# Beijing: 1.75%    averaged from: 1.28% 1.90% 1.58% 1.41% 2.21% 2.11% 1.94% 1.83% 1.49%
# Shanghai: 1.2566% averaged from:  1.04%  1.35% 1.19% 1.14% 1.59% 1.33% 1.30% 1.26% 1.11%
exportion_perc_cq = c(0.0115 * 2, 0.0127 * 2, 0.0044) # 1 normal, 2 chunyun 3.after quarantine  #### from calibrated
# exportion_perc_bj = c(0.0175, 0.0086)
exportion_perc_bj = c(0.0086, 0.0086, 0.0024)
# exportion_perc_sh = c(0.0125, 0.0066)
exportion_perc_sh = c(0.0066, 0.0066, 0.0029)

###### input for models: 
time = 130 #
# time = 365 * 2 # simulate two years (to record the date of daily incience go below 1)
# symptomatic_ratio = c(2/60, (301-2 - 200)/(301-2+318-4 -200- 265), (200)/(465))
ratio = c(0.09475153, 0.81885277, 1)
time_to_death = c(20, 20, 11) # from paper
z_0 = 43 * 2 # zoonotic force of infection, 86 cases per day
D_E = 6 # average 6 days of incubation period from Lancet
D_I = 8.4  - D_E # serial interval minus the mean latent period
E_I_beta_ratio = 0
# 
# cost_treatment = 4124.64
# cost_survived = c(131.94, 1317.08, 154.22) # from Appendix, and script calculate_GDP_loss.R
# cost_death = c(689680.60, 380892.36, 97043.36) # from Appendix, and script calculate_GDP_loss.R

#### contact ratio and contact mat estimated from https://royalsocietypublishing.org/doi/suppl/10.1098/rspb.2014.0268
### using script estimate_contact_matrix.R
contact_mat = rbind(c(21.09312, 20.44581, 15.02198),# normal
                    c(21.09312, 20.44581, 15.02198) * 1.0, # cny
                    # c(3, 3, 3))  #  during outbreak
                    c(21.09312, 20.44581, 15.02198) * 0.5)  #  during outbreak
contact_ratio = rbind(c(0.5928983 - (1-0.5928983)/3, 0.8353008 - (1- 0.8353008)/3, 0.29090466 - (1 - 0.29090466)/3), # normal 
                      # c(0, 0, 0)) # during cny
                      c(0.5928983 - (1-0.5928983)/3, 0.8353008 - (1- 0.8353008)/3, 0.29090466 - (1 - 0.29090466)/3)) 
                      # during cny, updated Sep 17th 2020


##### non-pharmaceutical interventions  ####
simulation_start = as.Date("12/01/19", "%m/%d/%y")
simulation_end = simulation_start + time 
chunyun_start = as.Date("01/10/20", "%m/%d/%y")
chunyun_end = chunyun_start + 40 - 1

milestone = as.Date("03/31/20", "%m/%d/%y")
nage = length(f_cq)
SD_end_dates = structure(rep(NA_real_, length(nage)), class="Date")


save.image(file = "model_inputs.rda")