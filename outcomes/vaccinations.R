# ** RELOADING DATA **

source("~/SEIRD model.R")
source("~/loading_data.R")

vacc_dec <- 1 - 0.765
cov <- 1

vaccine_SEIRD <- function(i, u, vacc){ 
  # p the proportion of cumulative clinical cases at which schools are closed
  # x the extent to which schools are closed (0 = no closure, 1 = full closure)
  SEIRD_17 <- function(time, init, parameters) {
    with(as.list(c(init, parameters)), {
      S = matrix(init[1:17], nrow = 17, ncol=1)
      E = matrix(init[18:34], nrow = 17, ncol=1)
      Ip = matrix(init[35:51], nrow = 17, ncol=1)
      Ic = matrix(init[52:68], nrow = 17, ncol=1)
      Is = matrix(init[69:85], nrow = 17, ncol=1)
      R = matrix(init[86:102], nrow = 17, ncol=1)
      D = matrix(init[103:119], nrow = 17, ncol=1)
      New = matrix(init[120:136], nrow = 17, ncol=1)
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
  contact <- final_matrix(i,u)
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1)
  if(vacc == T){
    clin_frac[15:17] <- clin_frac[15:17]*(1 - cov) + clin_frac[15:17]*vacc_dec*cov 
  }
  mu_c <- c(rep(0.0000161/0.29,3), rep(0.0000695/0.21,2), rep(0.000309/0.27,2), 
            rep(0.000844/0.33,2), rep(0.00161/0.40,2), 
            rep(0.00595/0.49,2), rep(0.0193/0.63,2), 0.0428/0.69, mean(c(0.0428, 0.0780))/0.69)
  mu_s = rep(0,17)
  parameters <- c(susc = 0.06, infec = 1/3, sympt = 1/2.1, 
                  rec_c = 1/2.9, rec_s = 1/5, xi = 0.5)  
  # susc the susceptibility to infection (leaving as 1 for now, for simplicity)
  # infec rate of becoming infectious (in either category),
  # sympt rate of moving from preclinical to clinical (developing symptoms), 
  # rec_c the clinical rate of recovery, rec_s the subclinical rate of recovery
  # these values are taken from nature paper estimates of disease periods
  # assuming wrongly for now that clinical and subclinical are equally infectious
  # mu_c, mu_s the probability that 'recovery' is death, in clinical and subclinical respectively
  # mu not dependent on IMD or age, as the worsened health is already accounted for in clin_frac
  init <- c(S1 = rural_age$Proportion[1 + 17*(i-1) + 170*(1-u)], S2 = rural_age$Proportion[2 + 17*(i-1) + 170*(1-u)], 
            S3 = rural_age$Proportion[3 + 17*(i-1) + 170*(1-u)], S4 = rural_age$Proportion[4 + 17*(i-1) + 170*(1-u)], 
            S5 = rural_age$Proportion[5 + 17*(i-1) + 170*(1-u)], S6 = rural_age$Proportion[6 + 17*(i-1) + 170*(1-u)], 
            S7 = rural_age$Proportion[7 + 17*(i-1) + 170*(1-u)], S8 = rural_age$Proportion[8 + 17*(i-1) + 170*(1-u)] - 1e-3, 
            S9 = rural_age$Proportion[9 + 17*(i-1) + 170*(1-u)], S10 = rural_age$Proportion[10 + 17*(i-1) + 170*(1-u)], 
            S11 = rural_age$Proportion[11 + 17*(i-1) + 170*(1-u)], S12 = rural_age$Proportion[12 + 17*(i-1) + 170*(1-u)],
            S13 = rural_age$Proportion[13 + 17*(i-1) + 170*(1-u)], S14 = rural_age$Proportion[14 + 17*(i-1) + 170*(1-u)], 
            S15 = rural_age$Proportion[15 + 17*(i-1) + 170*(1-u)], S16 = rural_age$Proportion[16 + 17*(i-1) + 170*(1-u)], 
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
  time <- 0:365
  
  model <- ode(y=init, times=time, func=SEIRD_17, parms=parameters)
  
  data <- data.frame(model)
  totals <- data.frame(time = data[,1], 
                       S = rowSums(data[,2:18]),
                       E = rowSums(data[,19:35]),
                       Ip = rowSums(data[,36:52]),
                       Ic = rowSums(data[,53:69]),
                       Is = rowSums(data[,70:86]),
                       R = rowSums(data[,87:103]),
                       D = rowSums(data[,104:120]),
                       New = rowSums(data[,121:137]))
  return(totals)
}

vaccine_age_standard <- function(i, u, vacc){  # i = IMD decile, u = urban (T)/rural (F)
  
  SEIRD_17 <- function(time, init, parameters) {
    with(as.list(c(init, parameters)), {
      S = matrix(init[1:17], nrow = 17, ncol=1)
      E = matrix(init[18:34], nrow = 17, ncol=1)
      Ip = matrix(init[35:51], nrow = 17, ncol=1)
      Ic = matrix(init[52:68], nrow = 17, ncol=1)
      Is = matrix(init[69:85], nrow = 17, ncol=1)
      R = matrix(init[86:102], nrow = 17, ncol=1)
      D = matrix(init[103:119], nrow = 17, ncol=1)
      New = matrix(init[120:136], nrow = 17, ncol=1)
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
  contact <- final_matrix(i,u)
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1) 
  if(vacc == T){
    clin_frac[15:17] <- clin_frac[15:17]*(1 - cov) + clin_frac[15:17]*vacc_dec*cov 
  }
  mu_c <- c(rep(0.0000161/0.29,3), rep(0.0000695/0.21,2), rep(0.000309/0.27,2), 
            rep(0.000844/0.33,2), rep(0.00161/0.40,2), 
            rep(0.00595/0.49,2), rep(0.0193/0.63,2), 0.0428/0.69, mean(c(0.0428, 0.0780))/0.69)
  mu_s = rep(0,17)
  parameters <- c(susc = 0.06, infec = 1/3, sympt = 1/2.1, 
                  rec_c = 1/2.9, rec_s = 1/5, xi = 0.5)  
  # susc the susceptibility to infection (leaving as 1 for now, for simplicity)
  # infec rate of becoming infectious (in either category),
  # sympt rate of moving from preclinical to clinical (developing symptoms), 
  # rec_c the clinical rate of recovery, rec_s the subclinical rate of recovery
  # these values are taken from nature paper estimates of disease periods
  # assuming wrongly for now that clinical and subclinical are equally infectious
  # mu_c, mu_s the probability that 'recovery' is death, in clinical and subclinical respectively
  # mu not dependent on IMD or age, as the worsened health is already accounted for in clin_frac
  init <- c(S1 = rural_age$Proportion[1 + 17*(i-1) + 170*(1-u)], S2 = rural_age$Proportion[2 + 17*(i-1) + 170*(1-u)], 
            S3 = rural_age$Proportion[3 + 17*(i-1) + 170*(1-u)], S4 = rural_age$Proportion[4 + 17*(i-1) + 170*(1-u)], 
            S5 = rural_age$Proportion[5 + 17*(i-1) + 170*(1-u)], S6 = rural_age$Proportion[6 + 17*(i-1) + 170*(1-u)], 
            S7 = rural_age$Proportion[7 + 17*(i-1) + 170*(1-u)], S8 = rural_age$Proportion[8 + 17*(i-1) + 170*(1-u)] - 1e-3, 
            S9 = rural_age$Proportion[9 + 17*(i-1) + 170*(1-u)], S10 = rural_age$Proportion[10 + 17*(i-1) + 170*(1-u)], 
            S11 = rural_age$Proportion[11 + 17*(i-1) + 170*(1-u)], S12 = rural_age$Proportion[12 + 17*(i-1) + 170*(1-u)],
            S13 = rural_age$Proportion[13 + 17*(i-1) + 170*(1-u)], S14 = rural_age$Proportion[14 + 17*(i-1) + 170*(1-u)], 
            S15 = rural_age$Proportion[15 + 17*(i-1) + 170*(1-u)], S16 = rural_age$Proportion[16 + 17*(i-1) + 170*(1-u)], 
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

# ** RESULTS **

# looking at clinical cases and deaths

vacc_res <- data.frame(i = rep(1:10, 8), u = rep(c(rep(T,20), rep(F,20)), 2),
                       vacc = rep(c(rep(T, 10), rep(F, 10)), 4),
                       as = c(rep(T, 40), rep(F, 40)),
                       cases = NA, deaths = NA)

for(k in 1:nrow(vacc_res)){
  if(vacc_res$as[k] == F){
    totals <- vaccine_SEIRD(vacc_res$i[k], vacc_res$u[k], vacc_res$vacc[k])
    vacc_res$cases[k] <- totals$R[366] + totals$D[366]
    vacc_res$deaths[k] <- totals$D[366]
  }else{
    totals <- vaccine_age_standard(vacc_res$i[k], vacc_res$u[k], vacc_res$vacc[k])
    sum_c <- sum((totals[366,87:103] + totals[366,104:120])*age_standard$prop[1:17 + 17*(1-vacc_res$u[k])]/rural_age$Proportion[1:17 + 17*(vacc_res$i[k]-1) + 170*(1-vacc_res$u[k])])
    sum_d <- sum(totals[366,104:120]*age_standard$prop[1:17 + 17*(1-vacc_res$u[k])]/rural_age$Proportion[1:17 + 17*(vacc_res$i[k]-1) + 170*(1-vacc_res$u[k])])
    vacc_res$cases[k] <- sum_c
    vacc_res$deaths[k] <- sum_d
  }
  if(k %% 5 == 0){
    print(round(k/nrow(vacc_res), 2))
  }
}

supp.labs <- c("Age-standardised", "Crude")
names(supp.labs) <- c(T,F)

cases <- vacc_res %>% ggplot() + 
  geom_line(aes(x=i, y=cases*1000, group=interaction(vacc, u, as), alpha=vacc, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  facet_grid(.~as, labeller = labeller(as = supp.labs)) + 
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total infections per 1000 population") + xlab('IMD decile') +
  scale_alpha_discrete(range = c(0.3, 1), name = "Vaccinations", labels = c("No","Yes")) 

deaths <- vacc_res %>% ggplot() + 
  geom_line(aes(x=i, y=deaths*1000, group=interaction(vacc, u, as), alpha=vacc, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  facet_grid(.~as, labeller = labeller(as = supp.labs)) + 
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total deaths per 1000 population") + xlab('IMD Decile') + 
  ylim(c(0,NA)) +
  scale_alpha_discrete(range = c(0.3, 1), name = "Vaccinations", labels = c("No","Yes")) 

cases + deaths + plot_layout(guides = 'collect', nrow=2)

cases_diff <- vacc_res %>% select(i, u, vacc, as, cases) %>% 
  pivot_wider(names_from = vacc, values_from = cases) %>% 
  mutate(diff = `FALSE` - `TRUE`) %>% ggplot() + 
  geom_line(aes(x=i, y=diff*1000, group=interaction(u, as), col=u, lty=as.factor(as)), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  scale_linetype_manual(values = c(1, 2), name = "Measure", labels = c("Crude","Age-standardised")) + 
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total infections prevented per 1000 population") + xlab('IMD Decile') + 
  ylim(c(0,NA))

deaths_diff <- vacc_res %>% select(i, u, vacc, as, deaths) %>% 
  pivot_wider(names_from = vacc, values_from = deaths) %>% 
  mutate(diff = `FALSE` - `TRUE`) %>% ggplot() + 
  geom_line(aes(x=i, y=diff*1000, group=interaction(u, as), col=u, lty=as.factor(as)), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  scale_linetype_manual(values = c(1, 2), name = "Measure", labels = c("Crude","Age-standardised")) + 
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Deaths prevented per 1000 population") + xlab('IMD decile') + 
  ylim(c(0,NA))

cases_diff + deaths_diff + plot_layout(guides = 'collect')

pop_sizes <- rural_age %>% group_by(IMD, rural) %>% summarise(sum(Population)) %>% 
  rename(pop = `sum(Population)`) %>% mutate(u = grepl('Urban', rural)) %>% 
  select(IMD, u, pop) %>% arrange(u, IMD)

vacc_sizes <- vacc_res %>% mutate(cases_acc = NA, deaths_acc = NA)
for(k in 1:nrow(vacc_sizes)){
  vacc_sizes$cases_acc[k] <- vacc_sizes$cases[k]*(pop_sizes %>% filter(IMD==vacc_sizes$i[k],
                                                                       u==vacc_sizes$u[k]))$pop
  vacc_sizes$deaths_acc[k] <- vacc_sizes$deaths[k]*(pop_sizes %>% filter(IMD==vacc_sizes$i[k],
                                                                       u==vacc_sizes$u[k]))$pop
}

cases_acc <- vacc_sizes %>% ggplot() + 
  geom_line(aes(x=i, y=cases_acc/1000, group=interaction(vacc, u, as), alpha=vacc, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  facet_grid(.~as, labeller = labeller(as = supp.labs)) + 
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total infections, 1000s") + xlab('IMD decile') +
  scale_alpha_discrete(range = c(0.3, 1), name = "Vaccinations", labels = c("No","Yes")) 

deaths_acc <- vacc_sizes %>% ggplot() + 
  geom_line(aes(x=i, y=deaths_acc/1000, group=interaction(vacc, u, as), alpha=vacc, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  facet_grid(.~as, labeller = labeller(as = supp.labs)) + 
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total deaths, 1000s") + xlab('IMD decile') + 
  ylim(c(0,NA)) +
  scale_alpha_discrete(range = c(0.3, 1), name = "Vaccinations", labels = c("No","Yes")) 

cases_acc + deaths_acc + plot_layout(guides = 'collect', nrow=2)

## R0 

fcn_vacc_R0 <- function(i, u, vacc){
  
  contact <- final_matrix(i,u)
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1)
  if(vacc == T){
    clin_frac[15:17] <- clin_frac[15:17]*(1 - cov) + clin_frac[15:17]*vacc_dec*cov 
  }
  parameters <- c(susc = 0.06, infec = 1/(3), sympt = 1/(2.1), 
                  rec_c = 1/(2.9), rec_s = 1/(5), xi = 0.5) 
  
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

vacc_r0 <- data.frame(i = rep(1:10, 4), u = rep(c(T,F), each=20),
                       vacc = rep(c(rep(T, 10), rep(F, 10)), 2),
                       r0 = NA)

for(k in 1:nrow(vacc_r0)){
  vacc_r0$r0[k] <- fcn_vacc_R0(vacc_r0$i[k], vacc_r0$u[k], vacc_r0$vacc[k])
}

vacc_r0 %>% pivot_wider(names_from = vacc, values_from = r0) %>%
  mutate(r0_diff = `FALSE` - `TRUE`) %>% ggplot() + 
  geom_line(aes(x=i, y=r0_diff, group=u, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Change in R0") + xlab('IMD decile') +
  ylim(c(0,NA))

vacc_r0 %>% ggplot() + 
  geom_line(aes(x=i, y=r0, group=interaction(u, vacc), col=u, alpha=vacc), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("R0") + xlab('IMD decile') +
  ylim(c(0,NA))


## VACCINES GIVEN

elderly_pop <- rural_age %>% filter(grepl('65', Age)|grepl('70', Age)|grepl('75', Age)) %>% 
  group_by(IMD, rural) %>% summarise(sum(Population)) %>% 
  rename(pop = `sum(Population)`) %>% mutate(u = grepl('Urban', rural)) %>% 
  select(IMD, u, pop) %>% arrange(u, IMD)

elderly_pop %>% ggplot() +
  geom_line(aes(x=IMD, y=pop, group=u, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Vaccine doses given") + xlab('IMD decile') +
  ylim(c(0,NA))

vacc_doses_elderly <- elderly_pop %>% #pivot_wider(names_from=u, values_from=pop) %>% 
  #mutate(total = `TRUE` + `FALSE`) %>% 
  ggplot() +
  geom_bar(aes(x=IMD, y=pop/1000000, fill=u), stat = "identity", position="stack") + 
  scale_fill_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Vaccine doses given, millions") + xlab('IMD Decile') +
  ylim(c(0,NA))


vacc_sizes %>% filter(as==F) %>% select(i, u, vacc, deaths_acc) %>% 
  pivot_wider(names_from = vacc, values_from = deaths_acc) %>% 
  mutate(diff = `FALSE` - `TRUE`) %>% ggplot() + 
  geom_line(aes(x=i, y=diff, group=u, col=u), lwd=0.8) + 
  scale_color_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total deaths prevented") + xlab('IMD Decile') + 
  ylim(c(0,NA)) 

deaths_prev_bar <- vacc_sizes %>% filter(as==F) %>% select(i, u, vacc, deaths_acc) %>% 
  pivot_wider(names_from = vacc, values_from = deaths_acc) %>% 
  mutate(diff = `FALSE` - `TRUE`) %>% ggplot() + 
  geom_bar(aes(x=i, y=diff, fill=u), stat = "identity", position="stack") + 
  scale_fill_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total deaths prevented") + xlab('IMD Decile') + 
  ylim(c(0,NA)) 

cases_prev_bar <- vacc_sizes %>% filter(as==F) %>% select(i, u, vacc, cases_acc) %>% 
  pivot_wider(names_from = vacc, values_from = cases_acc) %>% 
  mutate(diff = `FALSE` - `TRUE`) %>% ggplot() + 
  geom_bar(aes(x=i, y=diff, fill=u), stat = "identity", position="stack") + 
  scale_fill_brewer(palette='Set1', name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        legend.position = 'none',
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total infections prevented") + xlab('IMD Decile') + 
  ylim(c(0,NA)) 


deaths_diff + deaths_prev_bar + vacc_doses_elderly + plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')


cases_diff + cases_prev_bar + plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')



