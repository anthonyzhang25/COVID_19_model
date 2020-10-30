wuhan_simulation <- function(init_pop, time, beta, z_0, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu, symptomatic_ratio){
  
simulation_start = time_knot_vec$simulation_start; simulation_end = time_knot_vec$simulation_end
chunyun_start = time_knot_vec$chunyun_start; chunyun_end = time_knot_vec$chunyun_end
quarantine_start = time_knot_vec$quarantine_start; quarantine_end = time_knot_vec$quarantine_end
social_distancing_start = time_knot_vec$social_distancing_start
SD_end_dates = time_knot_vec$SD_end_dates


R0 = matrix(0, time) 
L_wi = matrix(0, time)
L_iw = matrix(0, time)
L_wc = matrix(0, time)
L_cw = matrix(0, time)
z = matrix(0, time)
# contact_rate_mat = matrix(0, time)
# beta_mat = matrix(0, time)
# chunyun_ind = matrix(NA, time)
# quarantine_ind = matrix(NA, time)
# social_distancing_ind = matrix(NA, time)

#### before chunyun (Jan 10)
L_wi_0 = 3633
L_iw_0 = 3546
L_wc_0 = 502013
L_cw_0 = 487310
### during chunyun
L_wc_1 = 717226
L_cw_1 = 810500
L_wi_1 = 3633
L_iw_1 = 3546
### during quarantine
# L_wc_q = 0
# L_cw_q = 0
# L_wi_q = 0
# L_iw_q = 0
### use reduction ratio
r_wc = 1*0.37/8.46 # from baidu qianxi
r_cw = 1*0.43/5.64 # from baidu qianxi
r_wi = 0
r_iw = 0
#### assign values
# calculate age distribution for Wuhan, obtained from https://www.citypopulation.de/php/china-hubei-admin_c.php?adm1id=4201
# 3 age group: Childen & Young adult (0 - 19), Adult (20 -59), the Elderly (60-90)
f = c(655178 + 1193294, 2279487 + 1486647 + 1668609 + 1261163, 718296 + 391084 + 131630)
f = f/sum(f)
# ### debug
# f = c(1/3, 1/3, 1/3)
C11 = matrix(0, time)
C12 = matrix(0, time)
C13 = matrix(0, time)
C21 = matrix(0, time)
C22 = matrix(0, time)
C23 = matrix(0, time)
C31 = matrix(0, time)
C32 = matrix(0, time)
C33 = matrix(0, time)
# beta = matrix(0, time, 3) # Columns by age group




# initialization
t = 1
# R0[t] = R0_0;  
L_wc[t] = L_wc_0; L_cw[t] = L_cw_0; L_wi[t] = L_wi_0; L_iw[t] = L_iw_0;  z[t] = z_0
c_vec = contact_mat[1,]; c_ratio = contact_ratio[1,]
denom = sum(c_vec * (1-c_ratio) * f)
C11[t] = (c_ratio[1] + (1- c_ratio[1])^2 * c_vec[1] * f[1]/denom) * c_vec[1]; C22[t] = (c_ratio[2] + (1- c_ratio[2])^2 * c_vec[2] * f[2]/denom) * c_vec[2]; 
C33[t] = (c_ratio[3] + (1- c_ratio[3])^2 * c_vec[3] * f[3]/denom) * c_vec[3]
C12[t] = (1- c_ratio[1])* (1-c_ratio[2]) * c_vec[2] * f[2]/denom * c_vec[1]; C13[t] = (1- c_ratio[1])* (1-c_ratio[3]) * c_vec[3] * f[3]/denom * c_vec[1];
C21[t] = (1- c_ratio[2])* (1-c_ratio[1]) * c_vec[1] * f[1]/denom * c_vec[2]; C23[t] = (1- c_ratio[2])* (1-c_ratio[3]) * c_vec[3] * f[3]/denom * c_vec[2];
C31[t] = (1- c_ratio[3])* (1-c_ratio[1]) * c_vec[1] * f[1]/denom * c_vec[3]; C32[t] = (1- c_ratio[3])* (1-c_ratio[2]) * c_vec[2] * f[2]/denom * c_vec[3];

# M = matrix(c(C11[t]*f[1]/f[1], C12[t]*f[1]/f[2], C13[t]*f[1]/f[3],
#              C21[t]*f[2]/f[1], C22[t]*f[2]/f[2], C23[t]*f[2]/f[3],
#              C31[t]*f[3]/f[1], C32[t]*f[3]/f[2], C33[t]*f[3]/f[3]),
#            nrow = 3)
# eig = eigen(M)
# # reverse engineer beta from the R0 and gamma 
# # beta[t] = R0[t]*gamma/max(Re(eig$values)) 
# beta = R0[t]/max(Re(eig$values))/D_I
### debug
# cat('beta is ', round(beta*100,2), "\n")



for (t in 2:time){
  today = t + simulation_start - 1
  R0[t] = R0[1];  L_wc[t] = L_wc[1]; L_cw[t] = L_cw[1]; L_wi[t] = L_wi[1]; L_iw[t] = L_iw[1];  z[t] = 0
  c_vec = contact_mat[1,]; c_ratio = contact_ratio[1,]
  
  if (today < as.Date("01/01/20", "%m/%d/%y")){
    z[t] = z_0
  }
  if (today <= chunyun_end && today >= chunyun_start){
    chunyun_ind = 1
    L_wc[t] = L_wc_1; L_cw[t] = L_cw_1; L_wi[t] = L_wi_1; L_iw[t] = L_iw_1
    c_vec = contact_mat[2,]; c_ratio = contact_ratio[2,]
    
  }
  if (today <= quarantine_end && today >= quarantine_start){
    quarantine_ind = 1
    # L_wc[t] = L_wc_q; L_cw[t] = L_cw_q; L_wi[t] = L_wi_q; L_iw[t] = L_iw_q
    L_wc[t] = L_wc[t] * r_wc; L_cw[t] = L_cw[t]*r_cw; L_wi[t] = L_wi[t]*r_wi; L_iw[t] = L_iw[t]*r_iw
  }
  
  ##### SD end differently by age 
  if (today <= SD_end_dates[1] && today >= social_distancing_start){
    #### contact_mat[3, ] is the row for s.d.
    c_vec[1] = contact_mat[3, 1]
  }
  if (today <= SD_end_dates[1] && today >= social_distancing_start){
    #### contact_mat[3, ] is the row for s.d.
    c_vec[2] = contact_mat[3, 2]
  }
  if (today <= SD_end_dates[1] && today >= social_distancing_start){
    #### contact_mat[3, ] is the row for s.d.
    c_vec[3] = contact_mat[3, 3]
  }
  
  
  denom = sum(c_vec * (1-c_ratio) * f)
  C11[t] = (c_ratio[1] + (1- c_ratio[1])^2 * c_vec[1] * f[1]/denom) * c_vec[1]; C22[t] = (c_ratio[2] + (1- c_ratio[2])^2 * c_vec[2] * f[2]/denom) * c_vec[2]; 
  C33[t] = (c_ratio[3] + (1- c_ratio[3])^2 * c_vec[3] * f[3]/denom) * c_vec[3]
  C12[t] = (1- c_ratio[1])* (1-c_ratio[2]) * c_vec[2] * f[2]/denom * c_vec[1]; C13[t] = (1- c_ratio[1])* (1-c_ratio[3]) * c_vec[3] * f[3]/denom * c_vec[1];
  C21[t] = (1- c_ratio[2])* (1-c_ratio[1]) * c_vec[1] * f[1]/denom * c_vec[2]; C23[t] = (1- c_ratio[2])* (1-c_ratio[3]) * c_vec[3] * f[3]/denom * c_vec[2];
  C31[t] = (1- c_ratio[3])* (1-c_ratio[1]) * c_vec[1] * f[1]/denom * c_vec[3]; C32[t] = (1- c_ratio[3])* (1-c_ratio[2]) * c_vec[2] * f[2]/denom * c_vec[3];
  
}


#### 020820:debug for contact mat
# for (t in 1:time){
#   today = t + simulation_start - 1
#   C12[t] = C13[t] = C21[t] = C23[t] = C31[t] = C32[t]  = 0
#   C11[t] = C22[t] = C33[t] = 20
#   if (today <= social_distancing_end && today >= social_distancing_start){
#   C11[t] = C22[t] = C33[t] = 20
#   }
# }
  
  
############ initialize the model, and run
t = 1
nage = length(f)
#### 021220: we added a few data frames
#### to record the rate of leaving in each component
#### 030920: add IA IS 
S = matrix(0, nage, time) 
E = matrix(0, nage, time)
# I = matrix(0, nage, time) 
IA = matrix(0, nage, time) 
IS = matrix(0, nage, time) 
RA = matrix(0, nage, time) 
RS = matrix(0, nage, time) 
D = matrix(0, nage, time) 

N_t = matrix(0, time)
S[, 1] = S[, 1] + init_pop * f
N_t[1] = init_pop
for (t in 2:time){
  today = t + simulation_start - 1
  # npop = S[,t - 1] + E[,t - 1] + IA[,t - 1] + IS[,t - 1] + RA[,t - 1] + RS[, t - 1] + D[,t - 1]
  npop = S[,t - 1] + E[,t - 1] + IA[,t - 1] + IS[,t - 1] + RA[,t - 1] + RS[, t - 1]
  N_t[t] = sum(npop + D[, t - 1])
  current_C = as.matrix(rbind(c(C11[t], C12[t], C13[t]),
                              c(C21[t], C22[t], C23[t]),
                              c(C31[t], C32[t], C33[t])))
  dS = -as.matrix(S[,t - 1])*(beta * current_C %*% as.matrix(IS[,t - 1])/npop + 
        E_I_beta_ratio*beta * current_C %*% as.matrix(E[,t - 1])/npop + 
        1*beta * current_C %*% as.matrix(IA[,t - 1])/npop +   z[t]/sum(npop)) + 
        (L_iw[t] + L_cw[t])/sum(S[,t - 1])*as.matrix(S[,t - 1]) - (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(S[,t - 1])# the derivative of S wrt time
  dE = +as.matrix(S[,t - 1])*(beta * current_C %*% as.matrix(IS[,t - 1])/npop + 
                                E_I_beta_ratio*beta * current_C %*% as.matrix(E[,t - 1])/npop + 
                                1*beta * current_C %*% as.matrix(IA[,t - 1])/npop +   z[t]/sum(npop)) - 
        as.matrix(E[,t - 1])/D_E - (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(E[,t - 1]) 
  # dI = +as.matrix(E[,t - 1])/D_E - as.matrix(I[,t - 1])/D_I -
  #       (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(I[,t - 1]) - mu*as.matrix(I[,t - 1])
  dIA = +as.matrix(E[,t - 1])/D_E *as.matrix(1 - symptomatic_ratio) - as.matrix(IA[,t - 1])/D_I -
        (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(IA[,t - 1]) 
  dIS = +as.matrix(E[,t - 1])/D_E * as.matrix(symptomatic_ratio) - as.matrix(IS[,t - 1])/D_I -
    (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(IS[,t - 1]) - mu*as.matrix(IS[,t - 1])
  dRS = as.matrix(IS[,t - 1])/D_I - (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(RS[,t - 1])
  dRA = as.matrix(IA[,t - 1])/D_I - (L_wi[t] + L_wc[t])/sum(npop) * as.matrix(RA[,t - 1])
  dD = mu * as.matrix(IS[,t - 1])
  
  #### we assume ths_start - ths_window - simulation_start is larger than 1
  
  S[,t] = S[,t - 1] + dS
  E[,t] = E[,t - 1] + dE
  # I[,t] = I[,t - 1] + dI
  IA[,t] = IA[, t - 1] + dIA
  IS[,t] = IS[, t - 1] + dIS
  #R[,t] = R[,t - 1] + dR
  RA[,t] = RA[, t - 1] + dRA
  RS[,t] = RS[, t - 1] + dRS
  D[,t] = D[,t - 1] + dD
  
}

S_mat = t(S); E_mat = t(E); IA_mat = t(IA); IS_mat = t(IS); RA_mat = t(RA); RS_mat = t(RS); D_mat = t(D)
  
###### calculate incidence 

### 030720: by age group
incidence = matrix(0, time)
# debug <- matrix(0, time, nage)
# for (i in 1:(time - 1)){
#   debug[i + 1,] = ((I_mat[i + 1,] + R_mat[i + 1,] + D_mat[i + 1,])
#                       - (I_mat[i,] + R_mat[i,] +D_mat[i,]))
# }

for (i in 1:(time - 1)){
  incidence[i + 1] = (sum(IA_mat[i + 1,] + IS_mat[i + 1,] + RA_mat[i + 1,] + RS_mat[i + 1,] + D_mat[i + 1,])
                                           - sum(IA_mat[i,] + IS_mat[i,] + RA_mat[i,] + RS_mat[i,] + D_mat[i,]))
}

###### calculate exported cases before travel history screening
exported = (array(0, dim = c(time, 3, 6))) # time: time steps 3: by age group  6: S E IA, IS, RA, RS

for (i in 1:time){
  for (age_ind in 1:3){
    exported[i, age_ind, 1] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + IA_mat[i,] +IS_mat[i,] + RA_mat[i,] + RS_mat[i,]))*
      (S_mat[i, age_ind])
    exported[i, age_ind, 2] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + IA_mat[i,] +IS_mat[i,] + RA_mat[i,] + RS_mat[i,]))*
      (E_mat[i, age_ind])
    exported[i, age_ind, 3] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + IA_mat[i,] +IS_mat[i,] + RA_mat[i,] + RS_mat[i,]))*
      (IA_mat[i, age_ind])
    exported[i, age_ind, 4] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + IA_mat[i,] +IS_mat[i,] + RA_mat[i,] + RS_mat[i,]))*
      (IS_mat[i, age_ind])
    exported[i, age_ind, 5] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + IA_mat[i,] +IS_mat[i,] + RA_mat[i,] + RS_mat[i,]))*
      (RA_mat[i, age_ind])
    exported[i, age_ind, 6] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + IA_mat[i,] +IS_mat[i,] + RA_mat[i,] + RS_mat[i,]))*
      (RS_mat[i, age_ind])
  }
}


# exported_E[i, 2] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + I_mat[i,] + R_mat[i,]))*
#   (E_mat[i, 2])
# exported_E[i, 3] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + I_mat[i,] + R_mat[i,]))*
#   (E_mat[i, 3])
# exported_I[i, 1] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + I_mat[i,] + R_mat[i,]))*
#   (I_mat[i, 1])
# exported_I[i, 2] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + I_mat[i,] + R_mat[i,]))*
#   (I_mat[i, 2])
# exported_I[i, 3] = (L_wc[i] /sum(S_mat[i,] + E_mat[i,] + I_mat[i,] + R_mat[i,]))*
#   (I_mat[i, 3])


   return (list("S" = S_mat, "E" = E_mat, "IA" = IA_mat, "IS" = IS_mat, "RA" = RA_mat, "RS" = RS_mat, "D" = D_mat, 'incidence' = incidence, 'exported' = exported, 'N_t' = N_t))
}
