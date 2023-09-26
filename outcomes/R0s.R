source("~/Desktop/MSc/Summer Project/R Code/loading_data.R")

epidemic_R0 <- function(i, u, a){
  
  contact <- final_matrix(i,u)
  clin_frac <- matrix(final_v2_dataframe[1:17 + 17*(i-1),3],
                      nrow=17, ncol=1) 
  parameters <- c(susc = 0.06, infec = 1/(a*3), sympt = 1/(a*2.1), 
                  rec_c = 1/(a*2.9), rec_s = 1/(a*5), xi = 0.5) 
  
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

R0s <- data.frame(i=rep(seq(1,10,1),2),
                  u=c(rep(T,10),rep(F,10)),
                  r=rep.int(NA,20))
for(j in 1:20){
  R0s$r[j] <- epidemic_R0(R0s$i[j], R0s$u[j], 1)
}
ggplot(R0s, aes(x=i,y=r,group=u,col=u)) + 
  geom_line(lwd=1) + theme_minimal() + scale_x_continuous(breaks=1:10)





# ** R0s AFTER FULL SCHOOL CLOSURES **

epidemic_R0_school_closure <- function(i, u, x){
  
  contact <- school_closure_matrix(i,u,x)
  clin_frac <- matrix(final_v2_dataframe[1:17 + 17*(i-1),3],
                      nrow=17, ncol=1) 
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

r0_data <- data.frame(i=rep(1:10,2),
                      u=rep(c(T,F),each=10),
                      R0=rep.int(NA,20),
                      R0school=rep.int(NA,20))
for(j in 1:20){
  r0_data$R0[j] <- epidemic_R0(r0_data$i[j],r0_data$u[j],1)
  r0_data$R0school[j] <- epidemic_R0_school_closure(r0_data$i[j],r0_data$u[j],1)
}
r0_data$diff <- r0_data$R0 - r0_data$R0school # reduction in R0 by closing schools
r0_data$frac <- r0_data$diff/r0_data$R0 # fractional change from original R0

ggplot(data = r0_data, aes(x=i, y=R0, group=u, col=u)) +
  geom_line(lwd=0.8, lty=3) + scale_color_brewer(palette='Set1',
                                                 name = "Geography", labels = c("Rural","Urban")) +
  geom_line(aes(x=i, y=R0school, group=u, col=u), lwd=1) + 
  theme_minimal() + ylim(c(1,3)) + scale_x_continuous(breaks=1:10) + 
  xlab('IMD Decile') + ylab('R0') +
  theme(text=element_text(size=16),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white"),
        legend.position='none')

ggplot(data = r0_data, aes(x=i, y=diff, group=u, col=u)) +
  geom_line(lwd=1) + scale_color_brewer(palette='Set1',
                                        name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + ylim(c(0,0.6)) + 
  scale_x_continuous(breaks=1:10) + 
  ylab('Reduction in R0') + 
  xlab('IMD Decile') +
  theme(text=element_text(size=16),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white"),
        legend.position='none')


## INVESTIGATING SCHOOL CLOSURES WITH AGE OR CLIN_FRAC CONSTANT

epidemic_R0_school_closure_cfc <- function(i, u, x){
  
  contact <- school_closure_matrix(i,u,x)
  clin_frac <- matrix(c(0.29,0.29,0.29,0.21,0.21,0.27,0.27,0.33,
                        0.33,0.4,0.4,0.49,0.49,0.63,0.63,0.69,0.69),
                      nrow=17, ncol=1) 
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

epidemic_R0_school_closure_ac <- function(i, u, x){
  
  contact <- school_closure_matrix(5,T,x)
  clin_frac <- matrix(final_v2_dataframe[1:17 + 17*(i-1),3],
                      nrow=17, ncol=1) 
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

r0_data <- data.frame(i=rep(1:10,2),
                      u=rep(c(T,F),each=10),
                      neither=rep.int(NA,20),
                      cf=rep.int(NA,20),
                      age=rep.int(NA,20))
for(j in 1:20){
  r0_data$neither[j] <- epidemic_R0_school_closure(r0_data$i[j],r0_data$u[j],0) - 
    epidemic_R0_school_closure(r0_data$i[j],r0_data$u[j],1)
  r0_data$cf[j] <- epidemic_R0_school_closure_cfc(r0_data$i[j],r0_data$u[j],0) -
    epidemic_R0_school_closure_cfc(r0_data$i[j],r0_data$u[j],1)
  r0_data$age[j] <- epidemic_R0_school_closure_ac(r0_data$i[j],r0_data$u[j],0) -
    epidemic_R0_school_closure_ac(r0_data$i[j],r0_data$u[j],1)
}

ggplot(r0_data, aes(x=i, y=neither, group=u, col=u)) + 
  geom_line(lwd=1, col='green') + theme_minimal() + 
  geom_line(aes(x=i, y=cf, group=u, col=u), lwd=1, col='blue') + 
  geom_line(aes(x=i, y=age, group=u, col=u), lwd=1, col='red') +
  scale_y_continuous(lim=c(0,0.65)) + 
  scale_x_continuous(breaks=1:10)













