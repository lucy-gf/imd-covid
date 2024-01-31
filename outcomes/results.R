# ** RESULTS **

deaths <- function(i, u){
  totals <- epidemic(i, u)
  return(totals$D[366])
}
final_size <- function(i, u){
  totals <- epidemic(i, u)
  return(totals$R[366] + totals$D[366])
}
peak <- function(i, u){
  totals <- epidemic(i, u)
  infections <- totals$Ip + totals$Ic + totals$Is
  peak <- max(infections)
  return(match(peak, infections))
}
peak_size <- function(i, u){
  totals <- epidemic(i, u)
  infections <- totals$Ip + totals$Ic + totals$Is
  peak <- max(infections)
  return(peak)
}
peak_clin_size <- function(i, u){
  totals <- epidemic(i, u)
  infections <- totals$Ic
  peak <- max(infections) 
  return(peak)
}
clinical <- function(i, u){
  totals <- epidemic(i, u)
  clinical <- totals$New[366] 
  return(clinical)
}

results <- data.frame(i=rep(seq(1,10,1),2),
                          u=c(rep(T,10),rep(F,10)),
                          d=rep.int(NA,20), # deaths
                          n=rep.int(NA,20), # final size
                          p=rep.int(NA,20), # peak date
                          ps=rep.int(NA,20), # peak size
                          pcs=rep.int(NA,20), # peak clinical size
                          c=rep.int(NA,20)) # total clinical cases

for(k in 1:20){
  results$d[k] <- deaths(results$i[k], results$u[k])
  results$n[k] <- final_size(results$i[k], results$u[k])
  results$p[k] <- peak(results$i[k], results$u[k])
  results$ps[k] <- peak_size(results$i[k], results$u[k])
  results$pcs[k] <- peak_clin_size(results$i[k], results$u[k])
  results$c[k] <- clinical(results$i[k], results$u[k])
}
results$prop_dead <- results$d/results$n



# standardised 
deaths_standardised <- function(i, u){
  totals <- epidemic_age_standard(i, u)
  sum <- sum(totals[366,104:120]*age_standard$prop[1:17 + 17*(1-u)]/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
  return(sum)
}
final_size_standardised <- function(i, u){
  totals <- epidemic_age_standard(i, u)
  sum <- sum(totals[366,87:103]*age_standard$prop[1:17 + 17*(1-u)]/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)]) +
    sum(totals[366,104:120]*age_standard$prop[1:17 + 17*(1-u)]/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
  return(sum)
}
clinical_standardised <- function(i, u){
  totals <- epidemic_age_standard(i, u)
  sum <- sum(totals[366,121:137]*age_standard$prop[1:17 + 17*(1-u)]/rural_age$Proportion[1:17 + 17*(i-1) + 170*(1-u)])
  return(sum)
}

results_age_standard <- data.frame(i=rep(seq(1,10,1),2),
                                   u=c(rep(T,10),rep(F,10)),
                                   d=rep.int(NA,20), # deaths
                                   n=rep.int(NA,20), # final size
                                   c=rep.int(NA,20))

for(k in 1:20){
  results_age_standard$d[k] <- deaths_standardised(results$i[k], results$u[k])
  results_age_standard$n[k] <- final_size_standardised(results$i[k], results$u[k])
  results_age_standard$c[k] <- clinical_standardised(results$i[k], results$u[k])
}
results_age_standard$prop_dead <- results_age_standard$d/results_age_standard$n

death2 <- ggplot(results_age_standard, aes(x=i, y=d*1000, group=u, color=u)) + 
  geom_line(lwd=0.8, lty=2) + scale_color_brewer(palette='Set1',
                                                 name = "Geography", labels = c("Rural","Urban")) +
  geom_line(data=results, aes(x=i, y=d*1000, group=u, color=u), lwd=1, lty=1) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total deaths per 1000 population") + xlab('IMD decile') + 
  scale_y_continuous(limits=c(0,12), breaks=seq(0,12,by=2))

cases2 <- ggplot(results_age_standard, aes(x=i, y=n*1000, group=u, color=u)) + 
  geom_line(lwd=0.8, lty=2) + scale_color_brewer(palette='Set1',
                                                 name = "Geography", labels = c("Rural","Urban")) +
  geom_line(data=results, aes(x=i, y=n*1000, group=u, color=u), lwd=1, lty=1) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total infections per 1000 population") + xlab('IMD decile') + 
  scale_y_continuous(limits=c(700,900))

clincases2 <- ggplot(results_age_standard, aes(x=i, y=c*1000, group=u, color=u)) + 
  geom_line(lwd=0.8, lty=2) + scale_color_brewer(palette='Set1',
                                                 name = "Geography", labels = c("Rural","Urban")) +
  geom_line(data=results, aes(x=i, y=c*1000, group=u, color=u), lwd=1, lty=1) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Total clinical cases per 1000 population") + xlab('IMD decile') +
  scale_y_continuous(limits=c(0,500))

ifr2 <- ggplot(results, aes(x=i, y=prop_dead*100, group=u, color=u)) + 
  geom_line(lwd=1,lty=1) + scale_color_brewer(palette='Set1',
                                                name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Infection fatality ratio (%)") + xlab('IMD decile') +
  scale_y_continuous(limits=c(0,1.5))

peak_clin_size_plot <- ggplot(results, aes(x=i, y=pcs*1000, group=u, color=u)) + 
  geom_line(lwd=1) + scale_color_brewer(palette='Set1',
                                        name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Peak clinical size per 1000 population") + xlab('IMD decile') +
  scale_y_continuous(limits=c(0,50))

R0_plot <- ggplot(R0s, aes(x=i,y=r,group=u,col=u)) + 
  geom_line(lwd=1) + scale_color_brewer(palette='Set1',name = "Geography", labels = c("Rural","Urban")) +
  theme_minimal() + scale_x_continuous(breaks = 1:10) +
  theme(axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1),
        panel.grid.minor.x = element_line(color = "white")) +
  ylab("Basic reproduction number") + xlab('IMD decile') +
  scale_y_continuous(limits=c(1.5,3))

cases2 + clincases2 + peak_clin_size_plot +
  death2 + ifr2 + R0_plot + plot_layout(nrow=2,guides='collect') +
  plot_annotation(tag_levels = 'a', tag_prefix = '(',
                  tag_suffix = ')  ')





#__________________________________________________________________________________________

## EPIDEMIC TRAJECTORIES ##


death_df <- data.frame(imd=rep(c(rep(1:10,each=200)),2),
                       u=rep(c(T,F),each=2000),
                       day=rep(1:200,20))

death_vec <- c(epidemic(1, T)[1:200,8],epidemic(2, T)[1:200,8],
               epidemic(3, T)[1:200,8],epidemic(4, T)[1:200,8],
               epidemic(5, T)[1:200,8],epidemic(6, T)[1:200,8],
               epidemic(7, T)[1:200,8],epidemic(8, T)[1:200,8],
               epidemic(9, T)[1:200,8],epidemic(10, T)[1:200,8],
               epidemic(1, F)[1:200,8],epidemic(2, F)[1:200,8],
               epidemic(3, F)[1:200,8],epidemic(4, F)[1:200,8],
               epidemic(5, F)[1:200,8],epidemic(6, F)[1:200,8],
               epidemic(7, F)[1:200,8],epidemic(8, F)[1:200,8],
               epidemic(9, F)[1:200,8],epidemic(10, F)[1:200,8])

death_df$D <- death_vec

death_urban <- ggplot(data = death_df[death_df$u==T,], aes(x=day, y=D*1000, group=interaction(imd,u), col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_color_viridis(option='A', name='IMD Decile') + 
  ylab('Deaths per 1,000 population') + xlab('Day') + scale_y_continuous(limits=c(0,11),breaks=1:11) + 
  ggtitle('Urban death rate')

death_rural <- ggplot(data = death_df[death_df$u==F,], aes(x=day, y=D*1000, group=interaction(imd,u), col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_color_viridis(option='A', name='IMD Decile') + 
  ylab('Deaths per 1,000 population') + xlab('Day') + scale_y_continuous(limits=c(0,11),breaks=1:11) + 
  ggtitle('Rural death rate')

death_urban + death_rural + plot_layout(guides='collect')

clin_inf_vec <- c(epidemic(1, T)[1:200,5],epidemic(2, T)[1:200,5],
               epidemic(3, T)[1:200,5],epidemic(4, T)[1:200,5],
               epidemic(5, T)[1:200,5],epidemic(6, T)[1:200,5],
               epidemic(7, T)[1:200,5],epidemic(8, T)[1:200,5],
               epidemic(9, T)[1:200,5],epidemic(10, T)[1:200,5],
               epidemic(1, F)[1:200,5],epidemic(2, F)[1:200,5],
               epidemic(3, F)[1:200,5],epidemic(4, F)[1:200,5],
               epidemic(5, F)[1:200,5],epidemic(6, F)[1:200,5],
               epidemic(7, F)[1:200,5],epidemic(8, F)[1:200,5],
               epidemic(9, F)[1:200,5],epidemic(10, F)[1:200,5])

death_df$C <- clin_inf_vec

ggplot(data = death_df[death_df$day %in% c(1:140),], aes(x=day, y=C, group=imd, col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_color_viridis() + facet_grid(.~u)

cum_clin_inf_vec <- c(epidemic(1, T)[1:200,9],epidemic(2, T)[1:200,9],
                  epidemic(3, T)[1:200,9],epidemic(4, T)[1:200,9],
                  epidemic(5, T)[1:200,9],epidemic(6, T)[1:200,9],
                  epidemic(7, T)[1:200,9],epidemic(8, T)[1:200,9],
                  epidemic(9, T)[1:200,9],epidemic(10, T)[1:200,9],
                  epidemic(1, F)[1:200,9],epidemic(2, F)[1:200,9],
                  epidemic(3, F)[1:200,9],epidemic(4, F)[1:200,9],
                  epidemic(5, F)[1:200,9],epidemic(6, F)[1:200,9],
                  epidemic(7, F)[1:200,9],epidemic(8, F)[1:200,9],
                  epidemic(9, F)[1:200,9],epidemic(10, F)[1:200,9])

death_df$cum_C <- cum_clin_inf_vec
death_df$imd <- as.numeric(death_df$imd)
supp.labs <- c("Urban", "Rural")
names(supp.labs) <- c(T,F)
ggplot(data = death_df[death_df$day %in% c(1:150),], aes(x=day, y=cum_C*1000, group=imd, col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_y_continuous(limits=c(0,450),breaks=seq(0,450,by=50)) + 
  scale_color_distiller(palette = "RdPu",breaks=seq(1,10,1),
                        lab=c('Most deprived','2','3','4','5','6','7','8','9','Least deprived')) +
  facet_grid(.~u,labeller = labeller(u = supp.labs)) +
  labs(color = "IMD Decile") + 
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        text = element_text(size = 14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1)) + 
  ylab('Cumulative clinical cases per 1000 population') + xlab('Day of epidemic')

subclin_inf_vec <- c(epidemic(1, T)[1:200,6],epidemic(2, T)[1:200,6],
                  epidemic(3, T)[1:200,6],epidemic(4, T)[1:200,6],
                  epidemic(5, T)[1:200,6],epidemic(6, T)[1:200,6],
                  epidemic(7, T)[1:200,6],epidemic(8, T)[1:200,6],
                  epidemic(9, T)[1:200,6],epidemic(10, T)[1:200,6],
                  epidemic(1, F)[1:200,6],epidemic(2, F)[1:200,6],
                  epidemic(3, F)[1:200,6],epidemic(4, F)[1:200,6],
                  epidemic(5, F)[1:200,6],epidemic(6, F)[1:200,6],
                  epidemic(7, F)[1:200,6],epidemic(8, F)[1:200,6],
                  epidemic(9, F)[1:200,6],epidemic(10, F)[1:200,6])

death_df$sub <- subclin_inf_vec

ggplot(data = death_df[death_df$u==T & death_df$day %in% c(1:140),], aes(x=day, y=sub, group=imd, col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_color_viridis()

ggplot(data = death_df[death_df$day %in% c(1:140),], aes(x=day, y=sub, group=interaction(imd, u), col=u)) + 
  geom_line(lwd=1) + theme_minimal()

recov_vec <- c(epidemic(1, T)[1:200,7],epidemic(2, T)[1:200,7],
                     epidemic(3, T)[1:200,7],epidemic(4, T)[1:200,7],
                     epidemic(5, T)[1:200,7],epidemic(6, T)[1:200,7],
                     epidemic(7, T)[1:200,7],epidemic(8, T)[1:200,7],
                     epidemic(9, T)[1:200,7],epidemic(10, T)[1:200,7],
                     epidemic(1, F)[1:200,7],epidemic(2, F)[1:200,7],
                     epidemic(3, F)[1:200,7],epidemic(4, F)[1:200,7],
                     epidemic(5, F)[1:200,7],epidemic(6, F)[1:200,7],
                     epidemic(7, F)[1:200,7],epidemic(8, F)[1:200,7],
                     epidemic(9, F)[1:200,7],epidemic(10, F)[1:200,7])

death_df$R <- recov_vec

ggplot(data = death_df[death_df$u==F & death_df$day %in% c(1:200),], aes(x=day, y=R, group=interaction(imd,u), col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_color_viridis()


death_df$inf <- rep.int(NA,nrow(death_df))

inf_vec <- c(rowSums(epidemic(1, T)[1:200,4:6]),rowSums(epidemic(2, T)[1:200,4:6]),
             rowSums(epidemic(3, T)[1:200,4:6]),rowSums(epidemic(4, T)[1:200,4:6]),
             rowSums(epidemic(5, T)[1:200,4:6]),rowSums(epidemic(6, T)[1:200,4:6]),
             rowSums(epidemic(7, T)[1:200,4:6]),rowSums(epidemic(8, T)[1:200,4:6]),
             rowSums(epidemic(9, T)[1:200,4:6]),rowSums(epidemic(10, T)[1:200,4:6]),
             rowSums(epidemic(1, F)[1:200,4:6]),rowSums(epidemic(2, F)[1:200,4:6]),
             rowSums(epidemic(3, F)[1:200,4:6]),rowSums(epidemic(4, F)[1:200,4:6]),
             rowSums(epidemic(5, F)[1:200,4:6]),rowSums(epidemic(6, F)[1:200,4:6]),
             rowSums(epidemic(7, F)[1:200,4:6]),rowSums(epidemic(8, F)[1:200,4:6]),
             rowSums(epidemic(9, F)[1:200,4:6]),rowSums(epidemic(10, F)[1:200,4:6]))

death_df$inf <- inf_vec

death_df$imd <- as.numeric(death_df$imd)

ggplot(data = death_df[death_df$day %in% c(1:150),], aes(x=day, y=inf*1000, group=imd, col=imd)) + 
  geom_line(lwd=1) + theme_minimal() + scale_y_continuous(limits=c(0,160),breaks=seq(0,160,by=20)) + 
  facet_grid(.~u,labeller = labeller(u = supp.labs)) +
  scale_color_distiller(palette = "RdPu",breaks=seq(1,10,1),
                        lab=c('Most deprived','2','3','4','5','6','7','8','9','Least deprived')) +
  labs(color = "IMD Decile") + 
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        text = element_text(size = 14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1)) + 
  ylab('Infections per 1000 population') + xlab('Day of epidemic') 




