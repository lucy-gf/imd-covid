
# ** RELOADING DATA **

source("~/Desktop/MSc/Summer Project/R Code/SEIRD model.R")

# ** THE MODEL **

epidemic_prev_deaths <- function(i, u){  # i = IMD decile, u = urban (T)/rural (F)
  
  SEIRD_17 <- function(times, init, parameters) {
    with(as.list(c(init, parameters)), {
      S = matrix(init[1:17], nrow = 17, ncol=1) # Susceptible
      E = matrix(init[18:34], nrow = 17, ncol=1) # Exposed
      Ip = matrix(init[35:51], nrow = 17, ncol=1) # Pre-clinical case
      Ic = matrix(init[52:68], nrow = 17, ncol=1) # Clinical case
      Is = matrix(init[69:85], nrow = 17, ncol=1) # Sub-clinical case
      R = matrix(init[86:102], nrow = 17, ncol=1) # Recovered 
      D = matrix(init[103:119], nrow = 17, ncol=1) # Dead
      New = matrix(init[120:136], nrow = 17, ncol=1) # Cumulative clinical cases
      force_inf <- susc * contact %*% ((Ip + Ic + xi*Is)/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
      dS <- (- force_inf * S)
      dE <-  (force_inf * S) - (infec * E)
      dIp <- (clin_frac * infec * E) - (sympt * Ip)
      dIc <- (sympt * Ip) - (rec_c * Ic)
      dIs <- ((1 - clin_frac) * infec * E) - (rec_s * Is)
      dR <- (rec_c * (1 - mu_c) * Ic) + (rec_s * (1 - mu_s) * Is)
      dD <- (rec_c * mu_c * Ic) + (rec_s * mu_s * Is)
      dNew <- (sympt * Ip) 
      return(list(c(dS, dE, dIp, dIc, dIs, dR, dD, dNew)))
    })
  }
  contact <- final_matrix(i,u) # contact matrix fitted to age structure
  clin_frac <- matrix(final_v2_dataframe$clin_frac[154:170],
                      nrow=17, ncol=1) 
  # age-specific probability of being assigned Ip over Is, dependent on IMD decile
  mu_c <- c(rep(0.0000161/0.29,3), rep(0.0000695/0.21,2), rep(0.000309/0.27,2), 
            rep(0.000844/0.33,2), rep(0.00161/0.40,2), 
            rep(0.00595/0.49,2), rep(0.0193/0.63,2), 0.0428/0.69, mean(c(0.0428, 0.0780))/0.69)
  mu_s = rep(0,17)
  parameters <- c(susc = 0.06, infec = 1/3, sympt = 1/2.1, 
                  rec_c = 1/2.9, rec_s = 1/5, xi = 0.5) 
  # susc the susceptibility to infection 
  # infec rate of becoming infectious from exposed (in either category),
  # sympt rate of moving from preclinical to clinical (developing symptoms), 
  # rec_c the clinical rate of recovery, rec_s the subclinical rate of recovery
  # mu_c, mu_s the probability of death, in clinical and subclinical respectively
  init <- c(S1 = rural_age$Proportion[1 + 17*(i-1) + 170*(1-u)], 
            S2 = rural_age$Proportion[2 + 17*(i-1) + 170*(1-u)], 
            S3 = rural_age$Proportion[3 + 17*(i-1) + 170*(1-u)], 
            S4 = rural_age$Proportion[4 + 17*(i-1) + 170*(1-u)], 
            S5 = rural_age$Proportion[5 + 17*(i-1) + 170*(1-u)], 
            S6 = rural_age$Proportion[6 + 17*(i-1) + 170*(1-u)], 
            S7 = rural_age$Proportion[7 + 17*(i-1) + 170*(1-u)], 
            S8 = rural_age$Proportion[8 + 17*(i-1) + 170*(1-u)] - 1e-3, 
            S9 = rural_age$Proportion[9 + 17*(i-1) + 170*(1-u)], 
            S10 = rural_age$Proportion[10 + 17*(i-1) + 170*(1-u)], 
            S11 = rural_age$Proportion[11 + 17*(i-1) + 170*(1-u)], 
            S12 = rural_age$Proportion[12 + 17*(i-1) + 170*(1-u)],
            S13 = rural_age$Proportion[13 + 17*(i-1) + 170*(1-u)], 
            S14 = rural_age$Proportion[14 + 17*(i-1) + 170*(1-u)], 
            S15 = rural_age$Proportion[15 + 17*(i-1) + 170*(1-u)], 
            S16 = rural_age$Proportion[16 + 17*(i-1) + 170*(1-u)], 
            S17 = rural_age$Proportion[17 + 17*(i-1) + 170*(1-u)], 
            E1 = 0, E2 = 0, E3 = 0, E4 = 0, E5 = 0, E6 = 0, E7 = 0, E8 = 0, E9 = 0, E10 = 0,
            E11 = 0, E12 = 0, E13 = 0, E14 = 0, E15 = 0, E16 = 0, E17 = 0,
            Ip1 = 0, Ip2 = 0, Ip3 = 0, Ip4 = 0, Ip5 = 0, Ip6 = 0, Ip7 = 0, Ip8 = 1e-3, Ip9 = 0, Ip10 = 0,
            Ip11 = 0, Ip12 = 0, Ip13 = 0, Ip14 = 0, Ip15 = 0, Ip16 = 0, Ip17 = 0, # 1e-3 is approximately 1 person in a 1500-person LSOA
            Ic1 = 0, Ic2 = 0, Ic3 = 0, Ic4 = 0, Ic5 = 0, Ic6 = 0, Ic7 = 0, Ic8 = 0, Ic9 = 0, Ic10 = 0,
            Ic11 = 0, Ic12 = 0, Ic13 = 0, Ic14 = 0, Ic15 = 0, Ic16 = 0, Ic17 = 0,
            Is1 = 0, Is2 = 0, Is3 = 0, Is4 = 0, Is5 = 0, Is6 = 0, Is7 = 0, Is8 = 0, Is9 = 0, Is10 = 0,
            Is11 = 0, Is12 = 0, Is13 = 0, Is14 = 0, Is15 = 0, Is16 = 0, Is17 = 0,
            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0, R6 = 0, R7 = 0, R8 = 0, R9 = 0, R10 = 0,
            R11 = 0, R12 = 0, R13 = 0, R14 = 0, R15 = 0, R16 = 0, R17 = 0,
            D1 = 0, D2 = 0, D3 = 0, D4 = 0, D5 = 0, D6 = 0, D7 = 0, D8 = 0, D9 = 0, D10 = 0,
            D11 = 0, D12 = 0, D13 = 0, D14 = 0, D15 = 0, D16 = 0, D17 = 0,
            New1 = 0, New2 = 0, New3 = 0, New4 = 0, New5 = 0, New6 = 0, New7 = 0, New8 = 0, New9 = 0, New10 = 0,
            New11 = 0, New12 = 0, New13 = 0, New14 = 0, New15 = 0, New16 = 0, New17 = 0) 
  times <- 0:365
  
  model <- ode(y=init, times=times, func=SEIRD_17, parms=parameters)
  data <- data.frame(model)
  totals <- data.frame(time = data[,1], 
                       S = rowSums(data[,2:18]),
                       E = rowSums(data[,19:35]),
                       Ip = rowSums(data[,36:52]),
                       Ic = rowSums(data[,53:69]),
                       Is = rowSums(data[,70:86]),
                       R = rowSums(data[,87:103]),
                       D = rowSums(data[,104:120]),
                       New = rowSums(data[,121:137])) #totals of cases etc. across all age groups
  return(totals)
}

epidemic_prev_deaths_age_specific <- function(i, u){  # i = IMD decile, u = urban (T)/rural (F)
  
  SEIRD_17 <- function(times, init, parameters) {
    with(as.list(c(init, parameters)), {
      S = matrix(init[1:17], nrow = 17, ncol=1) # Susceptible
      E = matrix(init[18:34], nrow = 17, ncol=1) # Exposed
      Ip = matrix(init[35:51], nrow = 17, ncol=1) # Pre-clinical case
      Ic = matrix(init[52:68], nrow = 17, ncol=1) # Clinical case
      Is = matrix(init[69:85], nrow = 17, ncol=1) # Sub-clinical case
      R = matrix(init[86:102], nrow = 17, ncol=1) # Recovered 
      D = matrix(init[103:119], nrow = 17, ncol=1) # Dead
      New = matrix(init[120:136], nrow = 17, ncol=1) # Cumulative clinical cases
      force_inf <- susc * contact %*% ((Ip + Ic + xi*Is)/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
      dS <- (- force_inf * S)
      dE <-  (force_inf * S) - (infec * E)
      dIp <- (clin_frac * infec * E) - (sympt * Ip)
      dIc <- (sympt * Ip) - (rec_c * Ic)
      dIs <- ((1 - clin_frac) * infec * E) - (rec_s * Is)
      dR <- (rec_c * (1 - mu_c) * Ic) + (rec_s * (1 - mu_s) * Is)
      dD <- (rec_c * mu_c * Ic) + (rec_s * mu_s * Is)
      dNew <- (sympt * Ip) 
      return(list(c(dS, dE, dIp, dIc, dIs, dR, dD, dNew)))
    })
  }
  contact <- final_matrix(i,u) # contact matrix fitted to age structure
  clin_frac <- matrix(final_v2_dataframe$clin_frac[154:170],
                      nrow=17, ncol=1) 
  # age-specific probability of being assigned Ip over Is, dependent on IMD decile
  mu_c <- c(rep(0.0000161/0.29,3), rep(0.0000695/0.21,2), rep(0.000309/0.27,2), 
            rep(0.000844/0.33,2), rep(0.00161/0.40,2), 
            rep(0.00595/0.49,2), rep(0.0193/0.63,2), 0.0428/0.69, mean(c(0.0428, 0.0780))/0.69)
  mu_s = rep(0,17)
  parameters <- c(susc = 0.06, infec = 1/3, sympt = 1/2.1, 
                  rec_c = 1/2.9, rec_s = 1/5, xi = 0.5) 
  # susc the susceptibility to infection 
  # infec rate of becoming infectious from exposed (in either category),
  # sympt rate of moving from preclinical to clinical (developing symptoms), 
  # rec_c the clinical rate of recovery, rec_s the subclinical rate of recovery
  # mu_c, mu_s the probability of death, in clinical and subclinical respectively
  init <- c(S1 = rural_age$Proportion[1 + 17*(i-1) + 170*(1-u)], 
            S2 = rural_age$Proportion[2 + 17*(i-1) + 170*(1-u)], 
            S3 = rural_age$Proportion[3 + 17*(i-1) + 170*(1-u)], 
            S4 = rural_age$Proportion[4 + 17*(i-1) + 170*(1-u)], 
            S5 = rural_age$Proportion[5 + 17*(i-1) + 170*(1-u)], 
            S6 = rural_age$Proportion[6 + 17*(i-1) + 170*(1-u)], 
            S7 = rural_age$Proportion[7 + 17*(i-1) + 170*(1-u)], 
            S8 = rural_age$Proportion[8 + 17*(i-1) + 170*(1-u)] - 1e-3, 
            S9 = rural_age$Proportion[9 + 17*(i-1) + 170*(1-u)], 
            S10 = rural_age$Proportion[10 + 17*(i-1) + 170*(1-u)], 
            S11 = rural_age$Proportion[11 + 17*(i-1) + 170*(1-u)], 
            S12 = rural_age$Proportion[12 + 17*(i-1) + 170*(1-u)],
            S13 = rural_age$Proportion[13 + 17*(i-1) + 170*(1-u)], 
            S14 = rural_age$Proportion[14 + 17*(i-1) + 170*(1-u)], 
            S15 = rural_age$Proportion[15 + 17*(i-1) + 170*(1-u)], 
            S16 = rural_age$Proportion[16 + 17*(i-1) + 170*(1-u)], 
            S17 = rural_age$Proportion[17 + 17*(i-1) + 170*(1-u)], 
            E1 = 0, E2 = 0, E3 = 0, E4 = 0, E5 = 0, E6 = 0, E7 = 0, E8 = 0, E9 = 0, E10 = 0,
            E11 = 0, E12 = 0, E13 = 0, E14 = 0, E15 = 0, E16 = 0, E17 = 0,
            Ip1 = 0, Ip2 = 0, Ip3 = 0, Ip4 = 0, Ip5 = 0, Ip6 = 0, Ip7 = 0, Ip8 = 1e-3, Ip9 = 0, Ip10 = 0,
            Ip11 = 0, Ip12 = 0, Ip13 = 0, Ip14 = 0, Ip15 = 0, Ip16 = 0, Ip17 = 0, # 1e-3 is approximately 1 person in a 1500-person LSOA
            Ic1 = 0, Ic2 = 0, Ic3 = 0, Ic4 = 0, Ic5 = 0, Ic6 = 0, Ic7 = 0, Ic8 = 0, Ic9 = 0, Ic10 = 0,
            Ic11 = 0, Ic12 = 0, Ic13 = 0, Ic14 = 0, Ic15 = 0, Ic16 = 0, Ic17 = 0,
            Is1 = 0, Is2 = 0, Is3 = 0, Is4 = 0, Is5 = 0, Is6 = 0, Is7 = 0, Is8 = 0, Is9 = 0, Is10 = 0,
            Is11 = 0, Is12 = 0, Is13 = 0, Is14 = 0, Is15 = 0, Is16 = 0, Is17 = 0,
            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0, R6 = 0, R7 = 0, R8 = 0, R9 = 0, R10 = 0,
            R11 = 0, R12 = 0, R13 = 0, R14 = 0, R15 = 0, R16 = 0, R17 = 0,
            D1 = 0, D2 = 0, D3 = 0, D4 = 0, D5 = 0, D6 = 0, D7 = 0, D8 = 0, D9 = 0, D10 = 0,
            D11 = 0, D12 = 0, D13 = 0, D14 = 0, D15 = 0, D16 = 0, D17 = 0,
            New1 = 0, New2 = 0, New3 = 0, New4 = 0, New5 = 0, New6 = 0, New7 = 0, New8 = 0, New9 = 0, New10 = 0,
            New11 = 0, New12 = 0, New13 = 0, New14 = 0, New15 = 0, New16 = 0, New17 = 0) 
  times <- 0:365
  
  model <- ode(y=init, times=times, func=SEIRD_17, parms=parameters)
  data <- data.frame(model)
  return(data)
}

epidemic_age_specific <- function(i, u){  # i = IMD decile, u = urban (T)/rural (F)
  
  SEIRD_17 <- function(times, init, parameters) {
    with(as.list(c(init, parameters)), {
      S = matrix(init[1:17], nrow = 17, ncol=1) # Susceptible
      E = matrix(init[18:34], nrow = 17, ncol=1) # Exposed
      Ip = matrix(init[35:51], nrow = 17, ncol=1) # Pre-clinical case
      Ic = matrix(init[52:68], nrow = 17, ncol=1) # Clinical case
      Is = matrix(init[69:85], nrow = 17, ncol=1) # Sub-clinical case
      R = matrix(init[86:102], nrow = 17, ncol=1) # Recovered 
      D = matrix(init[103:119], nrow = 17, ncol=1) # Dead
      New = matrix(init[120:136], nrow = 17, ncol=1) # Cumulative clinical cases
      force_inf <- susc * contact %*% ((Ip + Ic + xi*Is)/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
      dS <- (- force_inf * S)
      dE <-  (force_inf * S) - (infec * E)
      dIp <- (clin_frac * infec * E) - (sympt * Ip)
      dIc <- (sympt * Ip) - (rec_c * Ic)
      dIs <- ((1 - clin_frac) * infec * E) - (rec_s * Is)
      dR <- (rec_c * (1 - mu_c) * Ic) + (rec_s * (1 - mu_s) * Is)
      dD <- (rec_c * mu_c * Ic) + (rec_s * mu_s * Is)
      dNew <- (sympt * Ip) 
      return(list(c(dS, dE, dIp, dIc, dIs, dR, dD, dNew)))
    })
  }
  contact <- final_matrix(i,u) # contact matrix fitted to age structure
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1) 
  # age-specific probability of being assigned Ip over Is, dependent on IMD decile
  mu_c <- c(rep(0.0000161/0.29,3), rep(0.0000695/0.21,2), rep(0.000309/0.27,2), 
            rep(0.000844/0.33,2), rep(0.00161/0.40,2), 
            rep(0.00595/0.49,2), rep(0.0193/0.63,2), 0.0428/0.69, mean(c(0.0428, 0.0780))/0.69)
  mu_s = rep(0,17)
  parameters <- c(susc = 0.06, infec = 1/3, sympt = 1/2.1, 
                  rec_c = 1/2.9, rec_s = 1/5, xi = 0.5) 
  # susc the susceptibility to infection 
  # infec rate of becoming infectious from exposed (in either category),
  # sympt rate of moving from preclinical to clinical (developing symptoms), 
  # rec_c the clinical rate of recovery, rec_s the subclinical rate of recovery
  # mu_c, mu_s the probability of death, in clinical and subclinical respectively
  init <- c(S1 = rural_age$Proportion[1 + 17*(i-1) + 170*(1-u)], 
            S2 = rural_age$Proportion[2 + 17*(i-1) + 170*(1-u)], 
            S3 = rural_age$Proportion[3 + 17*(i-1) + 170*(1-u)], 
            S4 = rural_age$Proportion[4 + 17*(i-1) + 170*(1-u)], 
            S5 = rural_age$Proportion[5 + 17*(i-1) + 170*(1-u)], 
            S6 = rural_age$Proportion[6 + 17*(i-1) + 170*(1-u)], 
            S7 = rural_age$Proportion[7 + 17*(i-1) + 170*(1-u)], 
            S8 = rural_age$Proportion[8 + 17*(i-1) + 170*(1-u)] - 1e-3, 
            S9 = rural_age$Proportion[9 + 17*(i-1) + 170*(1-u)], 
            S10 = rural_age$Proportion[10 + 17*(i-1) + 170*(1-u)], 
            S11 = rural_age$Proportion[11 + 17*(i-1) + 170*(1-u)], 
            S12 = rural_age$Proportion[12 + 17*(i-1) + 170*(1-u)],
            S13 = rural_age$Proportion[13 + 17*(i-1) + 170*(1-u)], 
            S14 = rural_age$Proportion[14 + 17*(i-1) + 170*(1-u)], 
            S15 = rural_age$Proportion[15 + 17*(i-1) + 170*(1-u)], 
            S16 = rural_age$Proportion[16 + 17*(i-1) + 170*(1-u)], 
            S17 = rural_age$Proportion[17 + 17*(i-1) + 170*(1-u)], 
            E1 = 0, E2 = 0, E3 = 0, E4 = 0, E5 = 0, E6 = 0, E7 = 0, E8 = 0, E9 = 0, E10 = 0,
            E11 = 0, E12 = 0, E13 = 0, E14 = 0, E15 = 0, E16 = 0, E17 = 0,
            Ip1 = 0, Ip2 = 0, Ip3 = 0, Ip4 = 0, Ip5 = 0, Ip6 = 0, Ip7 = 0, Ip8 = 1e-3, Ip9 = 0, Ip10 = 0,
            Ip11 = 0, Ip12 = 0, Ip13 = 0, Ip14 = 0, Ip15 = 0, Ip16 = 0, Ip17 = 0, # 1e-3 is approximately 1 person in a 1500-person LSOA
            Ic1 = 0, Ic2 = 0, Ic3 = 0, Ic4 = 0, Ic5 = 0, Ic6 = 0, Ic7 = 0, Ic8 = 0, Ic9 = 0, Ic10 = 0,
            Ic11 = 0, Ic12 = 0, Ic13 = 0, Ic14 = 0, Ic15 = 0, Ic16 = 0, Ic17 = 0,
            Is1 = 0, Is2 = 0, Is3 = 0, Is4 = 0, Is5 = 0, Is6 = 0, Is7 = 0, Is8 = 0, Is9 = 0, Is10 = 0,
            Is11 = 0, Is12 = 0, Is13 = 0, Is14 = 0, Is15 = 0, Is16 = 0, Is17 = 0,
            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0, R6 = 0, R7 = 0, R8 = 0, R9 = 0, R10 = 0,
            R11 = 0, R12 = 0, R13 = 0, R14 = 0, R15 = 0, R16 = 0, R17 = 0,
            D1 = 0, D2 = 0, D3 = 0, D4 = 0, D5 = 0, D6 = 0, D7 = 0, D8 = 0, D9 = 0, D10 = 0,
            D11 = 0, D12 = 0, D13 = 0, D14 = 0, D15 = 0, D16 = 0, D17 = 0,
            New1 = 0, New2 = 0, New3 = 0, New4 = 0, New5 = 0, New6 = 0, New7 = 0, New8 = 0, New9 = 0, New10 = 0,
            New11 = 0, New12 = 0, New13 = 0, New14 = 0, New15 = 0, New16 = 0, New17 = 0) 
  times <- 0:365
  
  model <- ode(y=init, times=times, func=SEIRD_17, parms=parameters)
  data <- data.frame(model)
  return(data) 
} 
# this is the same as epidemic_age_standard in 'SEIRD model age-standardised.R'

###

pop <- sum(rural_age$tot_pop)/17 # 56550138, England population mid-2020

vec <- data.frame(i=rep(seq(1,10,1),2),
                             u=c(rep(T,10),rep(F,10)),
                             d=rep.int(NA,20)) 

for(k in 1:20){
  vec$d[k] <- epidemic_prev_deaths(vec$i[k],vec$u[k])$D[366]
}

res <- results %>% select(i,u,d) %>% 
  mutate(d_improved = vec$d) %>% 
  mutate(d_prev = d - d_improved) %>% 
  mutate(pop = rural_age$tot_pop[17*(1:20)])

# ggplot(data = res, aes(x = i, y = pop*d_prev, group = u, col = u)) + 
#   geom_line() + theme_minimal()

round(sum(res$pop * res$d_prev)) # 65,163

# deaths if all had health status of IMD10
acc_deaths <- sum(res$pop * res$d) # 405,695
min_deaths <- sum(res$pop * res$d_improved) # 340,532

prev_deaths <- acc_deaths - min_deaths # 65,163
prev_deaths - sum(res$pop * res$d_prev) # equiv

perc_prev <- prev_deaths*100/acc_deaths # 16.06%


## Age-specific results
###### At what ages are the 'prevented' deaths under health equity?

asmrs <- data.frame(i=rep(seq(1,10,1),2),u=c(rep(T,10),rep(F,10)))

asmrs <- cbind(asmrs, rbind(epidemic_prev_deaths_age_specific(1,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(2,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(3,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(4,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(5,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(6,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(7,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(8,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(9,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(10,T)[366,104:120],
                            epidemic_prev_deaths_age_specific(1,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(2,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(3,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(4,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(5,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(6,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(7,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(8,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(9,F)[366,104:120],
                            epidemic_prev_deaths_age_specific(10,F)[366,104:120]))

colnames(asmrs) <- c('IMD','Urban',"Under 1","1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29",
                     "30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                     "70 to 74","75+")

# Pivot wider
asmr <- asmrs %>% pivot_longer(!IMD:Urban) 
asmr$name <- factor(asmr$name, levels=unique(asmr$name))

# ASMR for each IMD and geography with underlying health equity:
ggplot(asmr, aes(x=name, y=value, group=IMD, col=IMD)) +
  geom_line() + theme_minimal() + facet_grid(.~Urban) +
  scale_color_viridis()


# now make another column using epidemic_age_specific

age_spec <- data.frame(i=rep(seq(1,10,1),2),u=c(rep(T,10),rep(F,10)))

age_spec <- cbind(age_spec, rbind(epidemic_age_specific(1,T)[366,104:120],
                                  epidemic_age_specific(2,T)[366,104:120],
                                  epidemic_age_specific(3,T)[366,104:120],
                                  epidemic_age_specific(4,T)[366,104:120],
                                  epidemic_age_specific(5,T)[366,104:120],
                                  epidemic_age_specific(6,T)[366,104:120],
                                  epidemic_age_specific(7,T)[366,104:120],
                                  epidemic_age_specific(8,T)[366,104:120],
                                  epidemic_age_specific(9,T)[366,104:120],
                                  epidemic_age_specific(10,T)[366,104:120],
                                  epidemic_age_specific(1,F)[366,104:120],
                                  epidemic_age_specific(2,F)[366,104:120],
                                  epidemic_age_specific(3,F)[366,104:120],
                                  epidemic_age_specific(4,F)[366,104:120],
                                  epidemic_age_specific(5,F)[366,104:120],
                                  epidemic_age_specific(6,F)[366,104:120],
                                  epidemic_age_specific(7,F)[366,104:120],
                                  epidemic_age_specific(8,F)[366,104:120],
                                  epidemic_age_specific(9,F)[366,104:120],
                                  epidemic_age_specific(10,F)[366,104:120]))

colnames(age_spec) <- c('IMD','Urban',"Under 1","1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29",
                     "30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                     "70 to 74","75+")

age_spe <- age_spec %>% pivot_longer(!IMD:Urban) 

asmr <- asmr %>% mutate(o_deaths = age_spe$value, 
                        prev = o_deaths - value,
                        pop = rural_age$tot_pop,
                        original_deaths = o_deaths*pop,
                        new_deaths = value*pop,
                        actual_prevented = prev*pop) 

# plot of reductions in ASMR 
ggplot(asmr, aes(x=name, y=prev, group=IMD, col=IMD)) +
  geom_line() + theme_minimal() + facet_grid(.~Urban) +
  scale_color_viridis()

# plot of prevented deaths across England
ggplot(asmr, aes(x=name, y=actual_prevented, group=IMD, col=IMD)) +
  geom_line() + theme_minimal() + facet_grid(.~Urban) +
  scale_color_viridis()

# sum over the IMD/geography to get an only age-specific table
vec <- unlist(unname(rural_age %>% group_by(Age) %>% 
  summarise(sum(Population)) %>% select(`sum(Population)`)))

age_prevs <- asmr %>% group_by(name) %>% 
  summarise_each(sum) %>% rename(Age = name) %>% 
  mutate(pop_in_eng = vec)

sum(age_prevs$actual_prevented) # 65,163
sum(age_prevs$actual_prevented[1:14]) # 29,847 aged under 65
sum(age_prevs$actual_prevented[15:17]) # 35,315 aged over 65
sum(age_prevs$actual_prevented[1:14])/sum(age_prevs$actual_prevented) # 45%
sum(age_prevs$original_deaths[1:14])/sum(age_prevs$original_deaths) # 24%
# so u65s make up 45% of the prevented deaths, despite only making 
# up 24% of the original deaths

# ggplot(age_prevs, aes(x=Age, y=actual_prevented)) +
#   geom_bar(stat = "identity") + theme_minimal()  #total prevented deaths in each age group
# ggplot(age_prevs, aes(x=Age, y=actual_prevented/pop_in_eng, group=NA)) +
#   geom_line() + theme_minimal() #proportion of each age group with prevented death
ggplot(age_prevs, aes(x=Age, y=100*actual_prevented/original_deaths, group=NA)) +
  geom_line(lwd=1, col='red') + theme_minimal() +
  ylab("Deaths prevented (%)") + scale_y_continuous(limits=c(0,40)) + 
  theme(strip.text.x = element_text(size = 14, face = "bold"),
                                     axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1,color='black'),
                                     text = element_text(size = 14),
                                     axis.text.y = element_text(color='black'))

#proportion of original deaths prevented in each age group

bar_data <- age_prevs %>% select(Age, new_deaths, actual_prevented) %>% 
  pivot_longer(!Age, names_repair = "unique") %>% rename(type = name)
  
ggplot(bar_data, aes(x=Age, y=value/1000, group=type, fill=type)) +
  geom_bar(stat = "identity", position = 'stack') + theme_minimal() +
  scale_fill_brewer(palette = 'Reds',
                    labels=c('Prevented deaths', 'Deaths under underlying health equity')) + 
  ylab('Deaths (1000s)') +
  scale_y_continuous(limits=c(0,200)) + 
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1,color='black'),
        text = element_text(size = 14),
        axis.text.y = element_text(color='black')) + 
  labs(fill = "")


