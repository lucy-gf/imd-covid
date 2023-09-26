library(dplyr)
#install.packages("readxl")
library(readxl)

# general health in each LSOA in England downloaded from 
# https://www.ons.gov.uk/datasets/TS037/editions/2021/versions/1/filter-outputs/a7d10c27-59cc-42e1-9f5e-21cd3d60f49e#get-data

lsoa_health <- read.csv("/Users/lucy/Desktop/MSc/Summer Project/datasets/TS037-2021-1-filtered-2023-06-08T11_17_48Z.csv")

# sub-domains of deprivation scores for each LSOA:
# File 2 at https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019

lsoa_doms <- read_xlsx("/Users/lucy/Desktop/MSc/Summer Project/datasets/File_2_-_IoD2019_Domains_of_Deprivation.xlsx",sheet=2)

lsoa_data <- data.frame(code = lsoa_doms$`LSOA code (2011)`,
                        name = lsoa_doms$`LSOA name (2011)`,
                        imd = lsoa_doms$`Index of Multiple Deprivation (IMD) Rank (where 1 is most deprived)`,
                        income = lsoa_doms$`Income Rank (where 1 is most deprived)`,
                        employment = lsoa_doms$`Employment Rank (where 1 is most deprived)`,
                        education = lsoa_doms$`Education, Skills and Training Rank (where 1 is most deprived)`,
                        health = lsoa_doms$`Health Deprivation and Disability Rank (where 1 is most deprived)`,
                        crime = lsoa_doms$`Crime Rank (where 1 is most deprived)`,
                        housing = lsoa_doms$`Barriers to Housing and Services Rank (where 1 is most deprived)`,
                        living = lsoa_doms$`Living Environment Rank (where 1 is most deprived)`,
                        health_prev = rep.int(NA,nrow(lsoa_doms)))

for(i in 1:nrow(lsoa_data)){
  lsoa_data$health_prev[i] = sum(lsoa_health$Observation[lsoa_health$Lower.layer.Super.Output.Areas.Code == lsoa_data$code[i] &
                                                    lsoa_health$General.health..6.categories..Code %in% c(1,2)])/
    sum(lsoa_health$Observation[lsoa_health$Lower.layer.Super.Output.Areas.Code == lsoa_data$code[i] &
                                lsoa_health$General.health..6.categories..Code %in% c(-8,1,2,3,4,5)])
    
}

#write.csv(lsoa_data, "/Users/lucy/Desktop/MSc/Summer Project/datasets/lsoa_data.csv", row.names=F)
lsoa_data <- read.csv("/Users/lucy/Desktop/MSc/Summer Project/MY DATASETS/lsoa_data.csv")

imd <- ggplot(lsoa_data, aes(x=imd, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("IMD rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
income <- ggplot(lsoa_data, aes(x=income, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Income rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
employment <- ggplot(lsoa_data, aes(x=employment, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Employment rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
education <- ggplot(lsoa_data, aes(x=education, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Education, skills and training rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
health <- ggplot(lsoa_data, aes(x=health, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Health and disability rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
crime <- ggplot(lsoa_data, aes(x=crime, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Crime rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
housing <- ggplot(lsoa_data, aes(x=housing, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Barriers to housing and services rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))
living <- ggplot(lsoa_data, aes(x=living, y=health_prev*100,color=imd)) +
  geom_point(size=0.75) + theme_minimal() + ylim(c(50,100)) +
  scale_color_viridis(limits=c(0,33719)) +
  ylab("Health prevalence (%)") + xlab("Living environment rank") + 
  labs(color = "IMD rank") + xlim(c(0,34753))

# library(patchwork)
png(filename="/Users/lucy/Desktop/LSHTM/Summer Project/Rough First Draft LaTeX/imd-components.png",
    width = 200, height = 225, units='mm', res = 300)
imd + income + employment + education + health + crime + 
  housing + living + plot_layout(nrow=4,ncol=2,guides='collect')
dev.off()







# ** CHECKING URBAN-RURAL STRATIFICATION **

lsoa_data$rural <- rep.int(NA,nrow(lsoa_data))
for(i in 1:nrow(lsoa_data)){
  lsoa_data$rural[i] <- rural$`Rural Urban Classification 2011 (2 fold)`[rural$`Lower Super Output Area 2011 Code` == lsoa_data$code[i]]
}

ggplot(lsoa_data[lsoa_data$rural=='Rural',], aes(x=imd, y=health_prev, color=rural)) +
  geom_point(size=0.75,shape=20) + theme_minimal() + ylim(c(0.5,1))

ggplot(lsoa_data[lsoa_data$rural=='Urban',], aes(x=crime, y=imd, color=rural)) +
  geom_point(size=0.75,shape=20) + theme_minimal()



