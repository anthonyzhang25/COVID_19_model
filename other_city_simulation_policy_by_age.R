other_city_simulation <- function(init_pop, f, time, beta, D_E, D_I, time_knot_vec, E_I_beta_ratio, contact_mat, contact_ratio, mu, exported, exportion_perc, ths_efficiency, ths_window, symptomatic_ratio){
  ####### for other cities, to simplify our analysis, we do not consider by age group and, in/outbound travels
  ###### no quarantine considered for other cities
  ##### f: the age group distribution
  ###### ths efficiency: travel history screening efficiency
  ##### ths_window: time window for ths. default 14 days
  simulation_start = time_knot_vec$simulation_start; simulation_end = time_knot_vec$simulation_end
  chunyun_start = time_knot_vec$chunyun_start; chunyun_end = time_knot_vec$chunyun_end
  ths_start = time_knot_vec$ths_start; ths_end = time_knot_vec$ths_end
  social_distancing_start = time_knot_vec$social_distancing_start
  # social_distancing_end = time_knot_vec$social_distancing_end
  SD_end_dates = time_knot_vec$SD_end_dates
  
  
  
  R0 = matrix(0, time) 
  z = matrix(0, time)
  # contact_rate_mat = matrix(0, time)
  # beta_mat = matrix(0, time)
  
  #### assume equal age distribution
  # f = c(1/3,1/3,1/3)
  ### 021220: load from the main script
  C11 = matrix(0, time)
  C12 = matrix(0, time)
  C13 = matrix(0, time)
  C21 = matrix(0, time)
  C22 = matrix(0, time)
  C23 = matrix(0, time)
  C31 = matrix(0, time)
  C32 = matrix(0, time)
  C33 = matrix(0, time)
  # travel_perc = matrix(0, time)
  # beta = matrix(0, time, 3) # Columns by age group
  # 
  
  
  
  # initialization
  t = 1
  # R0 = R0_0
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
  
  
  
  
  for (t in 2:time){
    today = t + simulation_start - 1
    # R0[t] = R0[1];
    c_vec = contact_mat[1,]; c_ratio = contact_ratio[1,]
    # travel_perc[t] = exportion_perc[1] # normal
    if (today <= chunyun_end && today >= chunyun_start){
      chunyun_ind = 1
      c_vec = contact_mat[2,]; c_ratio = contact_ratio[2,]
      # travel_perc[t] = exportion_perc[2] # during chunyun/CNY
    }
  
    ##### SD end differently by age 
    if (today <= SD_end_dates[1] && today >= social_distancing_start){
      #### contact_mat[3, ] is the row for s.d.
      c_vec[1] = contact_mat[3, 1]
    }
    if (today <= SD_end_dates[2] && today >= social_distancing_start){
      #### contact_mat[3, ] is the row for s.d.
      c_vec[2] = contact_mat[3, 2]
    }
    if (today <= SD_end_dates[3] && today >= social_distancing_start){
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
  
  
  
  
  
  

  
  
  
  
  ############ initialize the model, and run
  
  
  t = 1
  nage = length(f)
  #### 021220: we added a few data frames
  #### to record the rate of leaving in each component
  S = matrix(0, nage, time); Sw = matrix(0, nage, time); pSw_l = matrix(0, nage, time)
  E = matrix(0, nage, time); Ew = matrix(0, nage, time); pEw_l = matrix(0, nage, time)
  IA = matrix(0, nage, time); IwA = matrix(0, nage, time); pIwA_l = matrix(0, nage, time)
  IS = matrix(0, nage, time); IwS = matrix(0, nage, time); pIwS_l = matrix(0, nage, time)
  RA = matrix(0, nage, time); RwA = matrix(0, nage, time); # pRw_l = matrix(0, nage, time)
  RS = matrix(0, nage, time); RwS = matrix(0, nage, time); # pRw_l = matrix(0, nage, time)
  D = matrix(0, nage, time); Dw = matrix(0, nage, time); # pDw_l = matrix(0, nage, time)
  N_t = matrix(0, time)
  ##### 022820: add the quarantined components progressing
  Sq = Eq = IqA = IqS = RqA = RqS = Dq = array(0, c(nage, time, ths_window))
  # screening population size
  THS_count = matrix(0, nage, 6) # 4: SEIA, IS, RA, RS
  screening_counts = 0
  N_t[1] = init_pop
  S[, 1] = S[, 1] + init_pop * f
  Sw[, 1] = Sw[, 1] + exported[1, ,1] * exportion_perc[1]
  # THS_count[, 1] = 0
  # exported S is non-zero at t = 1 
  # exported format: dim = c(time, 3, 4))) # time: time steps 3: by age group  4: SEIR
  
  for (t in 2:time){
    today = t + simulation_start - 1
    travel_perc = exportion_perc[1]
    ths_ind = 0
    if (today >= chunyun_start && today <= chunyun_end) travel_perc = exportion_perc[2]
    if (today >= quarantine_start && today <= quarantine_end) travel_perc = exportion_perc[3]
    if (today >= ths_start && today <= ths_end) ths_ind = 1
    npop = S[,t - 1] + Sw[,t - 1] + E[,t - 1] + Ew[,t - 1] + IA[,t - 1] + IwA[,t - 1] + 
      IS[,t - 1] + IwS[,t - 1] + RA[,t - 1] + RwA[,t - 1] + RS[,t - 1] + RwS[,t - 1]
    # N_t[t] = sum(npop +  apply(Sq, 1, sum) + apply(Eq, 1, sum) + + apply(IqA, 1, sum) + apply(RqA, 1, sum) +
    #                + apply(IqS, 1, sum) + apply(RqS, 1, sum) + apply(Dq, 1, sum))
    N_t[t] = sum(npop)
    current_C = as.matrix(rbind(c(C11[t], C12[t], C13[t]),
                                c(C21[t], C22[t], C23[t]),
                                c(C31[t], C32[t], C33[t])))
    imported_t = exported[t, ,] *(1 - ths_efficiency * ths_ind) * travel_perc
    quarantined_t = exported[t, ,] *(ths_efficiency * ths_ind) * travel_perc
    screening_counts = screening_counts + sum(quarantined_t)
    Sq[, t - 1, 1] = Sq[, t - 1, 1] + quarantined_t[, 1]; Eq[, t - 1, 1] = Eq[, t - 1, 1] + quarantined_t[, 2] 
    IqA[, t - 1, 1] = IqA[, t - 1, 1] + quarantined_t[, 3]; IqS[, t - 1, 1] = IqS[, t - 1, 1] + quarantined_t[, 4]; 
    RqA[, t - 1, 1] = RqA[, t - 1, 1] + quarantined_t[, 5]; RqS[, t - 1, 1] = RqS[, t - 1, 1] + quarantined_t[, 6]
    ### on the day of ths start
    if (today == ths_start){
      ### cumu remaining at each component/state 
      ### (i.e., stayed at current component after the time window)
      Sw_remaining = as.matrix(Sw[, t - ths_window]) 
      Ew_remaining = as.matrix(Ew[, t - ths_window])
      IwA_remaining = as.matrix(IwA[, t - ths_window])
      IwS_remaining = as.matrix(IwS[, t - ths_window])
      for (delta in ths_window:1){
        # cat(delta, "\n")
        # cat("Sw_remaining = ",Sw_remaining,"\n")
        # cat("Ew_remaining = ",Ew_remaining,"\n")
        # cat("Iw_remaining = ",Iw_remaining,"\n")
        Sw_remaining = Sw_remaining * as.matrix(1 - pSw_l[,t - delta])
        Ew_remaining = Ew_remaining * as.matrix(1 - pEw_l[,t - delta])
        IwA_remaining = IwA_remaining * as.matrix(1 - pIwA_l[,t - delta])
        IwS_remaining = IwS_remaining * as.matrix(1 - pIwS_l[,t - delta])
        
      }
      THS_count[,1] = (Sw[, t - 1] - Sw_remaining) * ths_efficiency; Sw[, t - 1] = Sw[, t - 1] - THS_count[,1]
      THS_count[,2] = (Ew[, t - 1] - Ew_remaining) * ths_efficiency; Ew[, t - 1] = Ew[, t - 1] - THS_count[,2]
      THS_count[,3] = (IwA[, t - 1] - IwA_remaining) * ths_efficiency; IwA[, t - 1] = IwA[, t - 1] - THS_count[,3]
      THS_count[,4] = (IwS[, t - 1] - IwS_remaining) * ths_efficiency; IwS[, t - 1] = IwS[, t - 1] - THS_count[,4]
      THS_count[,5] = (RwA[, t - 1] - RwA[, t - ths_window]) * ths_efficiency; RwA[, t - 1] = RwA[, t- 1] - THS_count[,5]
      THS_count[,6] = (RwS[, t - 1] - RwS[, t - ths_window]) * ths_efficiency; RwS[, t - 1] = RwS[, t- 1] - THS_count[,6]
      # Rw[, t - 1] = Rw[, t - ths_window] + (Rw[, t - 1] - Rw[, t - ths_window]) * (1 - ths_efficiency)
      screening_counts = screening_counts + sum(THS_count)
      Sq[, t - 1, 1] = Sq[, t - 1, 1] + THS_count[, 1]; Eq[, t - 1, 1] = Eq[, t - 1, 1] + THS_count[, 2] 
      IqA[, t - 1, 1] = IqA[, t - 1, 1] + THS_count[, 3]; IqS[, t - 1, 1] = IqS[, t - 1, 1] + THS_count[, 4] 
      RqA[, t - 1, 1] = RqA[, t - 1, 1] + THS_count[, 5]; RqS[, t - 1, 1] = RqS[, t - 1, 1] + THS_count[, 6]
      }
    
    
    # ### on the day of ths quarantine end
    # if (today == ths_start + ths_window){
    #   Sw[, t - 1] = Sw[, t - 1] + Sq[, t - 1]
    #   Ew[, t - 1] = Ew[, t - 1] + Eq[, t - 1]
    #   Iw[, t - 1] = Iw[, t - 1] + Iq[, t - 1]
    #   Rw[, t - 1] = Rw[, t - 1] + Rq[, t - 1]
    #   Dw[, t - 1] = Dw[, t - 1] + Dq[, t - 1]
    # }
    # 
    dS = -as.matrix(S[,t - 1])*(beta * current_C %*% as.matrix(IS[,t - 1] + IwS[,t - 1])/npop + 
          E_I_beta_ratio*beta * current_C %*% as.matrix(E[,t - 1] + Ew[,t - 1] + IA[, t - 1] + IwA[, t - 1])/npop) # the derivative of S wrt time
    pSw_l[,t] =  (beta * current_C %*% as.matrix(IS[,t - 1] + IwS[,t - 1])/npop + 
                E_I_beta_ratio*beta * current_C %*% as.matrix(E[,t - 1] + Ew[,t - 1] + IA[, t - 1] + IwA[, t - 1])/npop)
    dSw = -as.matrix(Sw[,t - 1]) * pSw_l[t] + as.matrix(imported_t[,1]) # the derivative of S wrt time
    
    dE = +as.matrix(S[,t - 1])*(beta * current_C %*% as.matrix(IS[,t - 1] + IwS[,t - 1])/npop + 
          E_I_beta_ratio*beta * current_C %*% as.matrix(E[,t - 1] + Ew[,t - 1] + IA[, t - 1] + IwA[, t - 1])/npop) - as.matrix(E[,t - 1])/D_E
    pEw_l[,t] = 1/D_E
    dEw = +as.matrix(Sw[,t - 1])*(beta * current_C %*% as.matrix(IS[,t - 1] + IwS[,t - 1])/npop + 
          E_I_beta_ratio*beta * current_C %*% as.matrix(E[,t - 1] + Ew[,t - 1] + IA[, t - 1] + IwA[, t - 1])/npop) - as.matrix(Ew[,t - 1])/D_E  +  as.matrix(imported_t[,2])
    
    dIA = as.matrix(E[,t - 1])/D_E * as.matrix(1 - symptomatic_ratio) - as.matrix(IA[,t - 1])/D_I
    dIS = as.matrix(E[,t - 1])/D_E * as.matrix(symptomatic_ratio) - as.matrix(IS[,t - 1])/D_I - mu*as.matrix(IS[,t - 1])
    pIwA_l[,t] = 1/D_I
    pIwS_l[,T] = 1/D_I + mu
    dIwA = as.matrix(Ew[,t - 1])/D_E * as.matrix(1 - symptomatic_ratio) - as.matrix(IwA[,t - 1])/D_I + 
      as.matrix(imported_t[,3])
    dIwS = as.matrix(Ew[,t - 1])/D_E * as.matrix(symptomatic_ratio) - as.matrix(IwS[,t - 1])/D_I - mu*as.matrix(IwS[,t - 1]) + 
      as.matrix(imported_t[,4])
    # dR = I/D_I - (L_wi + L_wc)/npop * R                 # the derivative of R wrt time
    dRA = as.matrix(IA[,t - 1])/D_I
    dRS = as.matrix(IS[,t - 1])/D_I
    dRwA = as.matrix(IwA[,t - 1])/D_I + as.matrix(imported_t[,5])
    dRwS = as.matrix(IwS[,t - 1])/D_I + as.matrix(imported_t[,6])
    dD = mu * as.matrix(IS[,t - 1])
    dDw = mu * as.matrix(IwS[,t - 1])
    
    
    #### we assume ths_start - ths_window - simulation_start is larger than 1
    for (delta_t in 2:ths_window){
      dSq = 0
      dEq = - as.matrix(Eq[, t - 1, delta_t - 1])/D_E 
      dIqA = as.matrix(Eq[, t - 1, delta_t - 1])/D_E * as.matrix(1 - symptomatic_ratio) - as.matrix(IqA[, t - 1, delta_t - 1])/D_I
      dIqS = as.matrix(Eq[, t - 1, delta_t - 1])/D_E  * as.matrix(symptomatic_ratio)- as.matrix(IqS[, t - 1, delta_t - 1])/D_I - mu*as.matrix(IqS[, t - 1, delta_t - 1])
      dRqA = as.matrix(IqA[, t - 1, delta_t - 1])/D_I
      dRqS = as.matrix(IqS[, t - 1, delta_t - 1])/D_I
      dDq = mu * as.matrix(IqS[, t - 1, delta_t - 1])
      Sq[, t, delta_t] = Sq[, t - 1, delta_t - 1] + dSq
      Eq[, t, delta_t] = Eq[, t - 1, delta_t - 1] + dEq
      IqA[, t, delta_t] = IqA[, t - 1, delta_t - 1] + dIqA
      IqS[, t, delta_t] = IqS[, t - 1, delta_t - 1] + dIqS
      RqA[, t, delta_t] = RqA[, t - 1, delta_t - 1] + dRqA
      RqS[, t, delta_t] = RqS[, t - 1, delta_t - 1] + dRqS
      Dq[, t, delta_t] = Dq[, t - 1, delta_t - 1] + dDq
      
    }
    
    
    S[,t] = S[,t - 1] + dS; Sw[,t] = Sw[,t - 1] + dSw + Sq[, t, ths_window]
    E[,t] = E[,t - 1] + dE; Ew[,t] = Ew[,t - 1] + dEw + Eq[, t, ths_window]
    IA[,t] = IA[,t - 1] + dIA; IwA[,t] = IwA[,t - 1] + dIwA + IqA[, t, ths_window];
    IS[,t] = IS[,t - 1] + dIS; IwS[,t] = IwS[,t - 1] + dIwS;
    RA[,t] = RA[,t - 1] + dRA; RwA[,t] = RwA[,t - 1] + dRwA + RqA[, t, ths_window]
    RS[,t] = RS[,t - 1] + dRS; RwS[,t] = RwS[,t - 1] + dRwS + RqS[, t, ths_window]
    D[,t] = D[,t - 1] + dD; Dw[,t] = Dw[,t - 1] + dDw + Dq[, t, ths_window]
    
  }
  
  local_output = as.data.frame(cbind(t(S), t(E), t(IA), t(IS), t(RA), t(RS), t(D)))
  exported_output = as.data.frame(cbind(t(Sw), t(Ew), t(IwA),  t(IwS), t(RwA), t(RwS), t(Dw)))
  # quarantined_output = as.data.frame(cbind(t(Sq), t(Eq), t(Iq), t(Rq), t(Dq)))
  colnames(local_output) = c('S1', 'S2', 'S3', 'E1', 'E2', 'E3', 
                             'IA1', 'IA2', 'IA3', 'IS1', 'IS2', 'IS3', 
                             'RA1', 'RA2', 'RA3', 'RS1', 'RS2', 'RS3',
                             'D1', 'D2', 'D3')
  colnames(exported_output) = c('S1', 'S2', 'S3', 'E1', 'E2', 'E3', 
                                'IA1', 'IA2', 'IA3', 'IS1', 'IS2', 'IS3', 
                                'RA1', 'RA2', 'RA3', 'RS1', 'RS2', 'RS3',
                                'D1', 'D2', 'D3')
  # colnames(quarantined_output) = c('S1', 'S2', 'S3', 'E1', 'E2', 'E3', 
  #                               'I1', 'I2', 'I3', 'R1', 'R2', 'R3',
  #                               'D1', 'D2', 'D3')
  # output = as.data.frame(lsoda(inits, vt, Age_mortality_func, vparameters))
  # S_mat = cbind(local_output$S1, local_output$S2, local_output$S3,
  #               exported_output$S1, exported_output$S2, exported_output$S3)
  # E_mat = cbind(local_output$E1, local_output$E2, local_output$E3,
  #               exported_output$E1, exported_output$E2, exported_output$E3)
  # I_mat = cbind(local_output$I1, local_output$I2, local_output$I3,
  #               exported_output$I1, exported_output$I2, exported_output$I3)
  # R_mat = cbind(local_output$R1, local_output$R2, local_output$R3,
  #               exported_output$R1, exported_output$R2, exported_output$R3)
  # D_mat = cbind(local_output$D1, local_output$D2, local_output$D3,
  #               exported_output$D1, exported_output$D2, exported_output$D3)
  S_mat = cbind(local_output$S1 + exported_output$S1, local_output$S2 + exported_output$S2, local_output$S3 + exported_output$S3)
  E_mat = cbind(local_output$E1 + exported_output$E1, local_output$E2 + exported_output$E2, local_output$E3 + exported_output$E3)
  IA_mat = cbind(local_output$IA1 + exported_output$IA1, local_output$IA2 + exported_output$IA2, local_output$IA3 + exported_output$IA3)
  IS_mat = cbind(local_output$IS1 + exported_output$IS1, local_output$IS2 + exported_output$IS2, local_output$IS3 + exported_output$IS3)
  RA_mat = cbind(local_output$RA1 + exported_output$RA1, local_output$RA2 + exported_output$RA2, local_output$RA3 + exported_output$RA3)
  RS_mat = cbind(local_output$RS1 + exported_output$RS1, local_output$RS2 + exported_output$RS2, local_output$RS3 + exported_output$RS3)
  D_mat = cbind(local_output$D1 + exported_output$D1, local_output$D2 + exported_output$D2, local_output$D3 + exported_output$D3)
  
  
  ###### calculate incidence  
  incidence = matrix(0, dim(local_output[1]))
  # for (i in 1:length(incidence) - 1){
  #   incidence[i + 1] = (sum(I_mat[i + 1,] + R_mat[i + 1,] + D_mat[i + 1,])
  #                     - sum(I_mat[i,] + R_mat[i,] +D_mat[i,]))
  # }
  for (i in 1:(time - 1)){
    incidence[i + 1] = (sum(IA_mat[i + 1,] + IS_mat[i + 1,] + RA_mat[i + 1,] + RS_mat[i + 1,] + D_mat[i + 1,])
                        - sum(IA_mat[i,] + IS_mat[i,] + RA_mat[i,] + RS_mat[i,] + D_mat[i,]))
  }
  
  return (list("S" = S_mat, "E" = E_mat, "IA" = IA_mat, "IS" = IS_mat, "RA" = RA_mat, "RS" = RS_mat, "D" = D_mat, 'THS_count' = THS_count, 'incidence' = incidence, 'screening_counts' = screening_counts, 'N_t' = N_t))
  # return (list('local_output' = local_output, 'exported_output' = exported_output, 'THS_count' = THS_count, 'incidence' = incidence))
}