
deciles <- read.csv("/datasets/healthdatatables20211.csv",skip=4)
# found at:
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/healthandwellbeing/datasets/generalhealthbyagesexanddeprivationenglandandwales
deciles <- deciles[,1:11] #removing empty columns

#cleaning the data
deciles$Age <- as.character(deciles$Age)
deciles$Age <- factor(deciles$Age, levels=unique(deciles$Age))
deciles$Health.Status <- as.character(deciles$Health.Status)
deciles$Health.Status <- factor(deciles$Health.Status, levels=unique(deciles$Health.Status))
deciles$Age.specific.Percentage <- as.numeric(deciles$Age.specific.Percentage)
deciles$Age.specific.Percentage[is.na(deciles$Age.specific.Percentage)] <- 0 #setting values that were originally [low] to 0
deciles$IMD.Decile <- as.numeric(deciles$IMD.Decile)
deciles$Count <- as.numeric(gsub(",", "", deciles$Count))
deciles$Population <- as.numeric(gsub(",", "", deciles$Population))


# ** PLOTS ** 
supp.labs2 <- c("Most deprived decile", "Least deprived decile")
names(supp.labs2) <- c(1,10)
ggplot(deciles[deciles$IMD.Decile %in% c(1,10) & deciles$Sex =='Persons',] , aes(x = Age, y = Age.specific.Percentage, group = Health.Status, fill=Health.Status)) +
  geom_area(position='stack') + scale_fill_brewer(palette = "RdPu",name = "Health status") + theme_minimal() + 
  ylab(label='Population in health status (%)') + 
  facet_grid(.~IMD.Decile,labeller = labeller(IMD.Decile = supp.labs2)) + 
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1,color='black'),
        text = element_text(size = 14),
        axis.text.y = element_text(color='black'))


# ** CATEGORISING HEALTH INTO BINARY **

# making [75,79), [80,84), [85,89), [90,infty) into [75,infty)

binary_health <- data.frame('Age'=rep(deciles$Age[1:17],10),
                            'IMD.Decile'=as.numeric(rep(seq(1,10,1),each=17)))

# renaming the open bracket age group
binary_health$Age <- as.character(binary_health$Age)
binary_health[17*1:10, 1] <- "75+" 
deciles$Age <- as.character(deciles$Age)

binary_health$health_prev <- rep.int(NA, 170)

for(i in 1:170){
  if(binary_health$Age[i] != '75+'){
   binary_health$health_prev[i] <- as.numeric(sum(deciles$Count[deciles$Age == binary_health$Age[i] &
                                                                  deciles$IMD.Decile == binary_health$IMD.Decile[i] &
                                                                  deciles$Sex == 'Persons' &
                                                                    deciles$Health.Status %in% c('Very good','Good')]))/
     (as.numeric(deciles$Population[deciles$Age == binary_health$Age[i] &
                                                   deciles$IMD.Decile == binary_health$IMD.Decile[i] &
                                                   deciles$Sex == 'Persons' &
                                                   deciles$Health.Status %in% c('Very good')]))
  }
  if(binary_health$Age[i] == '75+'){
    binary_health$health_prev[i] <- as.numeric(sum(deciles$Count[deciles$Age %in% c('75 to 79','80 to 84','85 to 89','90+') &
                                                                                  deciles$IMD.Decile == binary_health$IMD.Decile[i] &
                                                                                  deciles$Sex == 'Persons' &
                                                                                  deciles$Health.Status %in% c('Very good','Good')]))/
                                           as.numeric(sum(deciles$Population[deciles$Age %in% c('75 to 79', '80 to 84', '85 to 89', '90+') &
                                 deciles$IMD.Decile == binary_health$IMD.Decile[i] &
                                 deciles$Sex == 'Persons' &
                                 deciles$Health.Status == 'Very good']))
  }
}

binary_health$Age <- as.character(binary_health$Age)
binary_health$Age <- factor(binary_health$Age, levels=unique(binary_health$Age))

#library(ggtext)

by_decile_plot <- ggplot(binary_health, aes(x=Age, y=health_prev, group=IMD.Decile, color=IMD.Decile)) +
  geom_line() + 
  scale_color_distiller(palette = "RdPu",breaks=seq(1,10,1),
                        lab=c('Most deprived','2','3','4','5','6','7','8','9','Least deprived')) + ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(color = "IMD decile") + 
  theme(plot.title = element_markdown(linewidth = 16), axis.title.y = element_markdown()) +
  ylab("Health prevalence") + xlab("Age group") + 
  theme(panel.background = element_rect(fill = "white",color=NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.07),
        panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.1),
        panel.grid.minor.y = element_line(color = "#A8BAC4", size = 0.1),
        axis.ticks.length.y = unit(0, "mm"), 
        axis.ticks.length.x = unit(0, "mm"),
        legend.key.size = unit(1, 'cm'),
        text = element_text(size = 14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1)); by_decile_plot

# ** OVERALL HEALTH PREVALENCE **

overall <- data.frame('decile'=1:10,
                    'prev'=rep.int(NA,10))
for(i in 1:10){
  overall$prev[i] <- sum(as.numeric(deciles$Count[deciles$IMD.Decile == i &
                                                  deciles$Sex == 'Persons' &
                                                  deciles$Health.Status %in% c('Very good','Good')]))/
    sum(as.numeric(deciles$Population[deciles$IMD.Decile == i &
                                        deciles$Sex == 'Persons' &
                                        deciles$Health.Status %in% c('Very good')]))
  
}

# ggplot(overall, aes(x=decile, y=prev)) + 
#   geom_line(lwd=1) + ylim(c(0,1)) + scale_colour_brewer(palette='RdPu') + 
#   labs(title = "Overall health prevalence",
#        subtitle = '(1 most deprived, 10 least deprived)') + 
#   theme(plot.title = element_markdown(size = 18), axis.title.y = element_markdown()) +
#   ylab("Health prevalence") + xlab("IMD Decile") + scale_x_continuous(breaks=seq(1,10,1)) +
#   theme(panel.background = element_rect(fill = "white",color=NA),
#         panel.grid = element_blank(),
#         panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.07),
#         panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.1),
#         axis.ticks.length.y = unit(0, "mm"), 
#         axis.ticks.length.x = unit(0, "mm"),
#         legend.key = element_rect(fill = NA))


#________________________________________________________________________________________________

# ** INCLUDING AVERAGE UK HEALTH **

all_UK <- data.frame('Age' = binary_health$Age[1:17],
                            'IMD.Decile' = rep("All UK",17))

all_UK$health_prev <- rep.int(NA, 17)

for(i in 1:17){
  if(binary_health$Age[i] != '75+'){
  all_UK$health_prev[i] <- sum(deciles$Count[deciles$Age == all_UK$Age[i] & deciles$Sex == 'Persons' & deciles$Health.Status %in% c('Very good','Good')])/
    sum(deciles$Population[deciles$Age == all_UK$Age[i] & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(binary_health$Age[i] == '75+'){
    all_UK$health_prev[i] <- sum(deciles$Count[deciles$Age %in% c('75 to 79', '80 to 84', '85 to 89', '90+') & 
                                                 deciles$Sex == 'Persons' & 
                                                 deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('75 to 79', '80 to 84', '85 to 89', '90+') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
}

binary_health_with_av <- rbind(binary_health,all_UK)

binary_health_with_av$IMD.Decile <- as.character(binary_health_with_av$IMD.Decile)
binary_health_with_av$IMD.Decile <- factor(binary_health_with_av$IMD.Decile, levels=unique(binary_health_with_av$IMD.Decile))

ggplot(binary_health_with_av, aes(x=Age, y=health_prev, group=IMD.Decile, color=IMD.Decile)) +
  geom_line() + scale_color_manual(values=c(rep('black',10),'red')) + theme_minimal() + ylim(c(0,1))
# Average UK health corresponds to health of 5th IMD decile

# __________________________________________________________________________________________


# ** CLIN_FRAC MAP **

v2 <- data.frame('age' = c('0-9','10-19','20-29','30-39',
                           '40-49','50-59','60-69','70+'),
                 'clin_frac' = c(0.29,0.21,0.27,0.33,0.40,0.49,0.63,0.69),
                 'health' = rep.int(NA,8))
ggplot(v2, aes(x=age, y=clin_frac, group=NA)) + 
  geom_line() + theme_minimal() + ylim(c(0,1))
bump <- 0.29 - 0.21 # increase in clin_frac in children under 10

for(i in 1:8){
  if(v2$age[i] == '0-9'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('Under 1', '1 to 4', '5 to 9') & 
                                                 deciles$Sex == 'Persons' & 
                                                 deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('Under 1', '1 to 4', '5 to 9') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '10-19'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('10 to 14', '15 to 19') & 
                                        deciles$Sex == 'Persons' & 
                                        deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('10 to 14', '15 to 19') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '20-29'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('20 to 24', '25 to 29') & 
                                        deciles$Sex == 'Persons' & 
                                        deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('20 to 24', '25 to 29') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '30-39'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('30 to 34', '35 to 39') & 
                                        deciles$Sex == 'Persons' & 
                                        deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('30 to 34', '35 to 39') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '40-49'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('40 to 44', '45 to 49') & 
                                        deciles$Sex == 'Persons' &
                                        deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('40 to 44', '45 to 49') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '50-59'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('50 to 54', '55 to 59') & 
                                        deciles$Sex == 'Persons' &
                                        deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('50 to 54', '55 to 59') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '60-69'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('60 to 64', '65 to 69') & 
                                        deciles$Sex == 'Persons' & 
                                        deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('60 to 64', '65 to 69') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
  if(v2$age[i] == '70+'){
    v2$health[i] <- sum(deciles$Count[deciles$Age %in% c('70 to 74', '75 to 79', '80 to 84', '85 to 89', '90+') & 
                                                 deciles$Sex == 'Persons' & 
                                                 deciles$Health.Status %in% c('Very good','Good')])/
      sum(deciles$Population[deciles$Age %in% c('70 to 74', '75 to 79', '80 to 84', '85 to 89', '90+') & deciles$Sex == 'Persons' & deciles$Health.Status == 'Very good'])
  }
} # adding health prevalences for 10-year age brackets

ggplot(v2[2:8,], aes(x=health, y=clin_frac, group=NA)) + 
  geom_point() + theme_minimal() + ylim(c(0,1)) + 
  xlim(c(0,1)) # plotting except without 0-10

# setting max and min domain values of the loess smoother
upper <- min(v2$health[v2$clin_frac==0.21],na.rm=T)
lower <- max(v2$health[v2$clin_frac==0.69],na.rm=T) 

loess <- loess(clin_frac ~ health, data=v2[2:8,], span=1)
func_loess <- function(i){
  return(predict(loess, i))
}

v2_loess <- function(x){
  ans <- min(func_loess(x),0.69)
  if(x <= lower){
    ans <- 0.69
  }
  if(x >= upper){
    ans <- 0.21
  }
  return(ans)
}

x <- seq(0,1,0.001)
v2_dataframe <- data.frame('health' = x, 'cf' = rep.int(NA,1001))
for(i in 1:1001){
  v2_dataframe$cf[i] <- v2_loess(v2_dataframe$health[i])
}
ggplot(v2_dataframe, aes(x=health, y=cf)) +
  geom_line(col='red',lty=2,lwd=1) + 
  geom_point(data=v2[2:8,], aes(x=health, y=clin_frac),col=1,pch=16,size=2.5) + 
  ylim(c(0,1)) +
  theme(plot.title = element_markdown(size = 18), 
        axis.title.x = element_markdown(size=13),
        axis.title.y = element_markdown(size=13)) +
  ylab("Clinical fraction") + xlab("Health prevalence") + 
  theme(panel.background = element_rect(fill = "white",color=NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.1),
        panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.1),
        panel.grid.minor.x = element_line(color = "#A8BAC4", size = 0.07),
        panel.grid.minor.y = element_line(color = "#A8BAC4", size = 0.07),
        axis.ticks.length.y = unit(0, "mm"), 
        axis.ticks.length.x = unit(0, "mm"),
        text = element_text(size = 14),
        axis.text.x = element_text(color = 1),
        axis.text.y = element_text(color = 1)) +
  geom_label(data = v2[2:4,], aes(x=health, y=clin_frac, label=age), hjust=1.15, vjust=0.9,
             label.size=0, size=4.5) + 
  geom_label(data = v2[5:8,], aes(x=health, y=clin_frac, label=age), hjust=0, vjust=-0.2,
             label.size=0, size=4.5)

final_v2 <- function(i,a){
  CH <- binary_health$health_prev[binary_health$Age==a & binary_health$IMD.Decile==i]
  CF <- v2_loess(CH)
  return(CF)
}

final_v2_dataframe <- data.frame('Ageband'=rep(c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29",
                                                              "30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                                                              "70 to 74","75+"),10),
                                              'IMD'=rep(seq(1,10,1),each=17),'clin_frac' = rep.int(NA,170))

for(i in 1:170){
  final_v2_dataframe$clin_frac[i] <- final_v2(final_v2_dataframe$IMD[i],
                                              final_v2_dataframe$Ageband[i])
  if(final_v2_dataframe$Ageband[i] %in% c('Under 1', '1 to 4', '5 to 9')){
    final_v2_dataframe$clin_frac[i] <- final_v2_dataframe$clin_frac[i] + bump
  }
}

final_v2_dataframe$Ageband <- factor(final_v2_dataframe$Ageband, levels=unique(final_v2_dataframe$Ageband))

ggplot(final_v2_dataframe, aes(x=Ageband, y = clin_frac, group = IMD, color=IMD)) + 
  geom_line() + 
  scale_color_distiller(palette = "RdPu",breaks=seq(1,10,1),
                        lab=c('Most deprived','2','3','4','5','6','7','8','9','Least deprived')) + 
  geom_line(data=final_v2_dataframe[final_v2_dataframe$IMD==2,],aes(x=Ageband, y = clin_frac, group = IMD, color=IMD)) +
  geom_line(data=final_v2_dataframe[final_v2_dataframe$IMD==1,],aes(x=Ageband, y = clin_frac, group = IMD, color=IMD)) +
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_continuous(breaks=seq(0,1,0.1), limits = c(0,0.75)) +
  labs(color = "IMD Decile") + 
  theme(plot.title = element_markdown(size = 18), 
        axis.title.x = element_markdown(size = 14),
        axis.title.y = element_markdown(size = 14)) +
  ylab("Clinical fraction") + xlab("Age group") + 
  theme(panel.background = element_rect(fill = "white",color=NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.1),
        panel.grid.major.y = element_line(color = "#A8BAC4", size = 0.1),
        panel.grid.minor.y = element_line(color = "#A8BAC4", size = 0.1),
        axis.ticks.length.y = unit(0, "mm"), 
        axis.ticks.length.x = unit(0, "mm"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=14),
        axis.text.x = element_text(color=1),
        axis.text.y = element_text(color=1))

    


