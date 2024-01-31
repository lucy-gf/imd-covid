#install.packages("deSolve")
require(deSolve)

# ** RELOADING DATA **

source("~/loading_data.R")

# ** THE MODEL **

epidemic_age_constant <- function(i){  # i = IMD decile, u = urban (T)/rural (F)
  
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
      force_inf <- susc * contact %*% ((Ip + Ic + xi*Is)/age_standard_all$prop)
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
  contact <- contact_constant # contact matrix fitted to age structure
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
  init <- c(S1 = age_standard_all$prop[1], 
            S2 = age_standard_all$prop[2], 
            S3 = age_standard_all$prop[3], 
            S4 = age_standard_all$prop[4], 
            S5 = age_standard_all$prop[5], 
            S6 = age_standard_all$prop[6], 
            S7 = age_standard_all$prop[7], 
            S8 = age_standard_all$prop[8] - 1e-3, 
            S9 = age_standard_all$prop[9], 
            S10 = age_standard_all$prop[10], 
            S11 = age_standard_all$prop[11], 
            S12 = age_standard_all$prop[12],
            S13 = age_standard_all$prop[13], 
            S14 = age_standard_all$prop[14], 
            S15 = age_standard_all$prop[15], 
            S16 = age_standard_all$prop[16], 
            S17 = age_standard_all$prop[17], 
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


# results 

r0_constant <- function(i, u){
  
  contact <- contact_constant
  clin_frac <- matrix(final_v2_dataframe[1:17 + 17*(i-1),3],
                      nrow=17, ncol=1) 
  parameters <- c(susc = 0.06, infec = 1/3, sympt = 1/2.1, 
                  rec_c = 1/2.9, rec_s = 1/5, xi = 0.5) 
  
  NGM <- matrix(data=rep.int(NA,17*17), nrow=17, ncol=17) #next-generation matrix
  for(k in 1:17){
    for(j in 1:17){
      NGM[k,j] <- unname(parameters["susc"])*contact[k,j]*
        (clin_frac[j,1]*(1/unname(parameters["sympt"]) + 1/unname(parameters["rec_c"])) + 
           (1-clin_frac[j,1])*(1/unname(parameters["rec_s"]))*unname(parameters["xi"]))
    }
  }
  vals <- eigen(NGM)$values
  real_vals <- Re(vals[Im(vals) == 0])
  R0 <- max(real_vals)
  
  return(R0)
}

results_constant <- data.frame(i=rep(seq(1,10,1),1),
                               u=rep(1,10),
                      d=rep.int(NA,10), # deaths
                      n=rep.int(NA,10), # final size
                      c=rep.int(NA,10), # total clinical cases
                      pcs=rep.int(NA,10), # peak clinical size
                      r0=rep.int(NA,10)) 

for(k in 1:10){
  vec <- unlist(unname(epidemic_age_constant(results_constant$i[k])[366,]))
  results_constant$d[k] <- vec[8]
  results_constant$n[k] <- sum(vec[7:8])
  results_constant$c[k] <- vec[9]
  results_constant$pcs[k] <- max(epidemic_age_constant(results_constant$i[k])$Ic)
  results_constant$r0[k] <- r0_constant(results_constant$i[k])
}
results_constant$d <- unlist(results_constant$d)
results_constant$n <- unlist(results_constant$n)
results_constant$c <- unlist(results_constant$c)
results_constant$pcs <- unlist(results_constant$pcs)
results_constant$prop_dead <- unlist(results_constant$d)/unlist(results_constant$n)

deaths <- ggplot(results_constant, aes(x=i, y=d*1000, group=u, color=u)) + 
  geom_line(lwd=1,col=1) + 
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total deaths per 1000 population") + xlab('IMD decile') + 
  scale_y_continuous(limits=c(0,10), breaks=seq(0,12,by=2))

size <- ggplot(results_constant, aes(x=i, y=n*1000, group=u, color=u)) + 
  geom_line(lwd=1,col=1) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total infections per 1000 population") + xlab('IMD decile') +
  scale_y_continuous(limits=c(700,900))

clin <- ggplot(results_constant, aes(x=i, y=c*1000, group=u, color=u)) + 
  geom_line(lwd=1,col=1) + 
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total clinical cases per 1000 population") + xlab('IMD decile') +
  scale_y_continuous(limits=c(0,500))

ifr <- ggplot(results_constant, aes(x=i, y=prop_dead*100, group=u, color=u)) + 
  geom_line(lwd=1,col=1) + 
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Infection fatality rate (%)") + xlab('IMD decile') +
  scale_y_continuous(limits=c(0,1.2))

pcs <- ggplot(results_constant, aes(x=i, y=pcs*1000, group=u, color=u)) + 
  geom_line(lwd=1,col=1) + 
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Peak clinical size per 1000 population") + xlab('IMD decile') +
  ylim(c(0,45))

r0 <- ggplot(results_constant, aes(x=i, y=r0, group=u, color=u)) + 
  geom_line(lwd=1,col=1) + 
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Basic reproduction number") + xlab('IMD decile') + 
  ylim(c(1,3))

size + clin + pcs + deaths + ifr + r0 + plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')













