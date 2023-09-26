#install.packages("deSolve")
require(deSolve)

# ** RELOADING DATA **

source("~/Desktop/MSc/Summer Project/R Code/loading_data.R")

# ** THE MODEL **

epidemic_length <- function(i, u, a){  # i = IMD decile, u = urban (T)/rural (F)
  
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
  parameters <- c(susc = 0.06, infec = 1/(a*3), sympt = 1/(a*2.1), 
                  rec_c = 1/(a*2.9), rec_s = 1/(a*5), xi = 0.5) 
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
  return(totals) # i have a different code file for a model that outputs age-specific case numbers etc.
}
epidemic_length_as <- function(i, u, a){  # i = IMD decile, u = urban (T)/rural (F)
  
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
  parameters <- c(susc = 0.06, infec = 1/(a*3), sympt = 1/(a*2.1), 
                  rec_c = 1/(a*2.9), rec_s = 1/(a*5), xi = 0.5) 
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

## RESULTS 

epi_data <- data.frame(a = rep(c(0.4,0.5,0.6,0.75,0.8,1,1.5,2,2.5,3),each=10),
                      imd = rep(1:10,20),
                      urban = c(rep(T, 100), rep(F, 100)),
                      cases = rep.int(NA, 200),
                      clin_cases = rep.int(NA, 200),
                      deaths = rep.int(NA, 200),
                      deaths_as = rep.int(NA, 200)) # as = age-standardised
for(k in 1:200){
  vec <- epidemic_length(epi_data$imd[k],epi_data$urban[k],epi_data$a[k])[366,]
  epi_data$cases[k] <- sum(unlist(vec[7:8]))
  epi_data$clin_cases[k] <- unlist(unname(vec[9]))
  epi_data$deaths[k] <- unlist(unname(vec[8]))
  epi_data$deaths_as[k] <- sum(epidemic_length_as(epi_data$imd[k],epi_data$urban[k],epi_data$a[k])[366,104:120]*age_standard$prop[1:17 + 17*(1-epi_data$urban[k])]/rural_age$Proportion[1:17 + 17*(epi_data$imd[k]-1) + 170*(1-epi_data$urban[k])])
}
#write.csv(epi_data, "/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/epi_sensitivity.csv", row.names=F)
epi_data <- read.csv("/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/epi_sensitivity.csv")

ggplot(data=epi_data, aes(x=imd, y=1000*cases, group=a, col=a)) +
  geom_line() + theme_minimal() + scale_color_viridis() + 
  scale_x_continuous(breaks=1:10) +
  facet_grid(.~urban) + ylim(c(0,1000))

ggplot(data=epi_data, aes(x=imd, y=1000*clin_cases, group=a, col=a)) +
  geom_line() + theme_minimal() + scale_color_viridis() + 
  scale_x_continuous(breaks=1:10) + 
  scale_y_continuous(limits=c(0,550),breaks=seq(0,550,by=50)) +
  facet_grid(.~urban)

ggplot(data=epi_data, aes(x=imd, y=1000*deaths, group=a, col=a)) +
  geom_line() + theme_minimal() +scale_color_viridis() + 
  geom_line(aes(x=imd, y=1000*deaths_as, group=a, col=a),lty=2) + 
  scale_x_continuous(breaks=1:10) + scale_y_continuous(limits=c(0,15),breaks=0:15) + 
  facet_grid(.~urban)

epi_data_NA <- epi_data
colnames(epi_data_NA)[7] <- "deaths"
epi_data_2 <- rbind(epi_data_NA[,c(1,2,3,6)],epi_data_NA[,c(1,2,3,7)])
epi_data_2$measure <- c(rep("Crude",200),rep("Age-standardised \nby geography",200))
epi_data_2$measure <- factor(epi_data_2$measure, levels = c("Crude", "Age-standardised \nby geography"))
supp.labs <- c("Urban", "Rural")
names(supp.labs) <- c(T,F)
ggplot(data=epi_data_2[epi_data_2$a %in% c(0.4,0.6,0.75,1,1.5,3),], aes(x=imd, y=1000*deaths, 
                                        group=interaction(a,measure), col=a,
                                        linetype=measure)) +
  geom_line(lwd=0.7) + theme_minimal() + 
  scale_color_viridis(direction=-1, breaks=c(0.4,seq(1,3,0.5)), limits=c(0.4,3)) +
  scale_x_continuous(breaks=1:10) + scale_y_continuous(limits=c(0,15),breaks=0:15) + 
  ylab('Deaths per 1000 population') + xlab('IMD Decile') + 
  labs(color = "Length multiplier",linetype="Measure") +
  scale_linetype_manual(values=c(1,2)) + 
  facet_grid(.~urban,labeller = labeller(urban = supp.labs)) +
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        text=element_text(size=12),
                axis.text.x = element_text(color=1),
                axis.text.y = element_text(color=1),
                panel.grid.minor.x = element_line(color = "white"))














