# ** RELOADING DATA **

source("~/Desktop/MSc/Summer Project/R Code/SEIRD model.R")

M_school <- socialmixr::contact_matrix(polymod, countries = "United Kingdom", 
                           age.limits = c(0, 1, seq(5,75,5)),
                           return.demography = T, filter = list(cnt_school=1))$matrix
demog2005 <- socialmixr::contact_matrix(polymod, countries = "United Kingdom", 
                            age.limits = c(0, 1, seq(5,75,5)),
                            return.demography = T)$demography #GB population structure 2005
N <- demog2005$population
G_school <- matrix(nrow=17, ncol=17)
for(i in 1:17){
  for(j in 1:17){
    G_school[i,j] = M_school[i,j]*sum(N)/N[j]
  }
}

school_closure_matrix <- function(i, u, x){
  final_matrix <- matrix(nrow=17, ncol=17)
  character <- 'Urban'
  if(u == F){
    character <- 'Rural'
  }
  vec <- rural_age$Population[rural_age$IMD == i & 
                                rural_age$rural == character]
  for(i in 1:17){
    for(j in 1:17){
      final_matrix[i,j] = (G[i,j] - x*G_school[i,j])*vec[j]/sum(vec)
    }
  }
  return(final_matrix)
}


epidemic_intervention <- function(i, u, p, x){ 
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
      contact <- contact0
      if(sum(New) > p){ 
        contact <- contact1
      }
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
  contact0 <- final_matrix(i,u)
  contact1 <- school_closure_matrix(i,u,x)
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1) 
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
# i = IMD
# u = T (urban)/F (rural)
# p = proportion of the population who experience a clinical Covid case to trigger school closure
# x = extent to which schools are closed (x in [0,1])


epidemic_age_standard_school <- function(i, u){  # i = IMD decile, u = urban (T)/rural (F)
  
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
      contact <- contact0
      if(sum(New) > 0.05){ 
        contact <- contact1
      }
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
  contact0 <- final_matrix(i,u)
  contact1 <- school_closure_matrix(i,u,1)
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1) 
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


# epidemic plots:

# matplot(epidemic_intervention(5, T, 0.05, 1)[1:120,c(2,5,6,7,8)], type="l", lty=1, main="SEIR model", xlab="Time",col=2:8)
# legend <- colnames(epidemic_intervention(5, T, 5, 1))[c(2,5,6,7,8)]
# legend("right", legend=legend, col=2:8, lty = 1)
# 
# matplot(epidemic_intervention(8, F, 0.05, 1)[2:121,3], type="l", lty=1, lwd=1, 
#         xlab="Day of outbreak",col=1,ylab='Exposed proportion of population')
# 
# matplot(epidemic_intervention(10, T, 0.03, 1)[2:101,2:8] - epidemic_intervention(10, T, 0.03, 1)[1:100,2:8], type="l", lty=1, lwd=2, main="SEIRD model", 
#         xlab="Time",col=c(2:7,1), 
#         ylab='Proportion')
# legend <- colnames(epidemic_intervention(1, T, 0.03,1))[2:8]
# legend("right", legend=legend, col=c(2:8,1), 
#        lty = 1, lwd=2)


# ** RESULTS **

# Reduction in R0 is in R0s.R

# Reduction in deaths:
deaths <- function(i, u){
  totals <- epidemic(i, u)
  return(totals$D[366])
}
deaths_school <- function(i, u, p, x){
  totals <- epidemic_intervention(i, u, p, x)
  return(totals$D[366])
}

deaths_standardised <- function(i, u){
  totals <- epidemic_age_standard_school(i, u)
  sum <- sum(totals[366,104:120]*age_standard$prop[1:17 + 17*(1-u)]/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
  return(sum)
}

school_deaths <- data.frame(i=rep(1:10,2),
                           u=rep(c(T,F),each=10),
                           school_open=rep.int(NA,20),
                           school_closed=rep.int(NA,20),
                           age_standard_school=rep.int(NA,20))

for(j in 1:20){
  school_deaths$school_open[j] <- deaths(school_deaths$i[j], 
                                             school_deaths$u[j])
  school_deaths$school_closed[j] <- deaths_school(school_deaths$i[j], 
                                                      school_deaths$u[j],
                                                     0.05,1)
  school_deaths$age_standard_school[j] <- deaths_standardised(school_deaths$i[j], 
                                                  school_deaths$u[j])
}
school_deaths$age_standard_school_open <- results_age_standard$d
school_deaths$diff <- school_deaths$school_open - school_deaths$school_closed
school_deaths$standard_diff <- school_deaths$age_standard_school_open - school_deaths$age_standard_school

#write.csv(school_deaths, "/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/school_deaths.csv", row.names=F)

ggplot(data=school_deaths, aes(x=i,y=diff*1000,
                              group=u,col=u)) +
  geom_line(lwd=1) + theme_minimal() + ylim(c(0,0.4)) +
  scale_x_continuous(breaks=1:10) + 
  geom_line(aes(x=i, y=standard_diff*1000, group=u, color=u), lwd=1, lty=2) +
  xlab('IMD Decile') + ylab('Reduction in deaths per 1,000 people') + 
  scale_color_brewer(palette='Set1',name = "Geography", labels = c("Rural","Urban")) +
  theme(text=element_text(size=16),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white"))

ggplot(data=school_deaths, aes(x=i,y=(school_open)*1000,
                                   group=u,col=u)) +
  geom_line(lwd=0.8,lty=2) + 
  geom_line(aes(x=i,y=(school_closed)*1000,group=u,col=u), lwd=1) + 
  theme_minimal() + 
  scale_color_brewer(palette='Set1',name = "Geography", labels = c("Rural","Urban")) +
  ylab('Deaths per 1,000 people') + ylim(c(0,12)) +
  scale_x_continuous(breaks=1:10)





#__________________________________________________________________________________________________
# ** CHANGING THRESHOLD **

n = 10 * 4 * 10 * 2
school_sens <- data.frame(u = rep(c(T,F),each=n/2),
                          imd = rep(1:10,each=n/(2*10)),
                          x = rep(rep(c(0.25,0.5,0.75,1),each=10),n/(4*10)),
                          p = rep(seq(0,0.45,by=0.05),n/10),
                          deaths = rep.int(NA,n))

for(i in 1:800){
  school_sens$deaths[i] <- epidemic_intervention(school_sens$imd[i], school_sens$u[i], school_sens$p[i], school_sens$x[i])$D[366]
}

open_vec <- school_deaths$school_open
school_sens$open <- rep(open_vec, each=40)
school_sens$diff <- school_sens$open - school_sens$deaths

#write.csv(school_sens, "/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/school_sens.csv", row.names=F)
school_sens <- read.csv("/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/school_sens.csv")

# ggplot(school_sens[school_sens$x == 0.5 & school_sens$p == 1,],
#        aes(x = imd, y = diff*1000, group = u, colour = u)) +
#   geom_line(lwd=1) + theme_minimal() + 
#   scale_x_continuous(breaks=1:10) + ylim(c(0,0.5)) + 
#   xlab('IMD Decile') + ylab('Reduction in deaths per 1,000 people') + 
#   scale_color_brewer(palette='Set1',name = "Geography", labels = c("Rural","Urban")) 

ggplot(school_sens[school_sens$p<0.35 & school_sens$x==1,],aes(x = imd, y = diff*1000, group = u, colour = u)) +
  geom_line(lwd=1) + theme_bw() + 
  scale_x_continuous(breaks=c(1:10)) + ylim(c(0,0.5)) + 
  xlab('IMD Decile') + ylab('Reduction in deaths per 1,000 people') + 
  theme(text=element_text(size=12),
        panel.grid.minor.x = element_line(color = "white")) + 
  scale_color_brewer(palette='Set1',name = "Geography", labels = c("Rural","Urban")) + 
  facet_grid(.~p) 



#__________________________________________________________________________________________________
# ** DAYS INSTEAD **

epidemic_intervention_days <- function(i, u, d, x){ 
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
      contact <- contact0
      if(time > d){ 
        contact <- contact1
      }
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
  contact0 <- final_matrix(i,u)
  contact1 <- school_closure_matrix(i,u,x)
  clin_frac <- matrix(final_v2_dataframe$clin_frac[seq(1,17,1) + 17*(i-1)],
                      nrow=17, ncol=1) 
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


matplot(epidemic_intervention_days(8, F, 60, 1)[2:121,3], type="l", lty=1, lwd=1, 
        xlab="Day of outbreak",col=1,ylab='Exposed proportion of population')

n = 5 * 4 * 10 * 2
school_sens_day <- data.frame(u = rep(c(T,F),each=n/2),
                          imd = rep(1:10,each=20),
                          x = rep(rep(c(0.25,0.5,0.75,1),each=5),n/(4*10)),
                          d = rep(seq(10,50,by=10),n/10),
                          deaths = rep.int(NA,n))

for(i in 1:400){
  school_sens_day$deaths[i] <- epidemic_intervention_days(school_sens_day$imd[i], school_sens_day$u[i], school_sens_day$d[i], school_sens_day$x[i])$D[366]
}

open_vec <- school_deaths$school_open
school_sens_day$open <- rep(open_vec, each=20)
school_sens_day$diff <- school_sens_day$open - school_sens_day$deaths
#write.csv(school_sens_day, "/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/school_sens_day.csv", row.names=F)
school_sens_day <- read.csv("/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/school_sens_day.csv")

ggplot(school_sens_day[school_sens_day$x==1,],aes(x = imd, y = diff*1000, group = u, colour = u)) +
  geom_line(lwd=1) + theme_bw() + 
  scale_x_continuous(breaks=c(1:10)) + ylim(c(0,0.5)) + 
  xlab('IMD Decile') + ylab('Reduction in deaths per 1,000 people') + 
  scale_color_brewer(palette='Set1',name = "Geography", labels = c("Rural","Urban")) +
  theme(text=element_text(size=12),panel.grid.minor.x = element_line(color = "white")) + 
  facet_grid(.~d) 








