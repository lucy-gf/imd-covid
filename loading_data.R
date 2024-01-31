#install.packages("deSolve")
require(deSolve)

# loading necessary data

final_v2_dataframe <- read.csv("/clin_frac.csv")
final_v2_dataframe$Ageband <- factor(final_v2_dataframe$Ageband, levels=unique(final_v2_dataframe$Ageband))

rural_age <- read.csv("/rural_age.csv")
rural_age$Age <- factor(rural_age$Age, levels=unique(rural_age$Age))

G <- read.csv("/G.csv")
final_matrix <- function(i, u){
  final_matrix <- matrix(nrow=17, ncol=17)
  character <- 'Urban'
  if(u == F){
    character <- 'Rural'
  }
  vec <- rural_age$Population[rural_age$IMD == i & 
                                rural_age$rural == character]
  for(i in 1:17){
    for(j in 1:17){
      final_matrix[i,j] = G[i,j]*vec[j]/sum(vec)
    }
  }
  return(final_matrix)
}

urban_IMD_distr <- read.csv("/urban_IMD_distr.csv")

results <- read_csv("/results.csv")
results_age_standard <- read_csv("/results_age_standard.csv")

age_standard <- data.frame(age=rep(rural_age$Age[1:17],2),
                           urban=c(rep(T,17),rep(F,17)),
                           prop=rep.int(NA,34))
for(i in 1:17){
  age_standard$prop[i] <- sum(rural_age$Population[(1:10)*17 + (i-17)])/
    sum(rural_age$tot_pop[(1:10)*17 + (i-17)])
}
for(i in 18:34){
  age_standard$prop[i] <- sum(rural_age$Population[(11:20)*17 + (i-34)])/
    sum(rural_age$tot_pop[(11:20)*17 + (i-34)])
}
age_standard$age <- factor(age_standard$age, levels=unique(age_standard$age))

age_standard_all <- data.frame(age=rural_age$Age[1:17],
                               prop=rep.int(NA,17))
for(i in 1:17){
  age_standard_all$prop[i] <- sum(rural_age$Population[(1:20)*17 + (i-17)])/
    sum(rural_age$tot_pop[(1:20)*17 + (i-17)])
}

age_standard_all$age <- factor(age_standard_all$age, levels=unique(age_standard_all$age))
age_standard_all$prop <- as.numeric(age_standard_all$prop)

contact_constant <- matrix(nrow=17, ncol=17)
vec <- age_standard_all$prop
for(i in 1:17){
  for(j in 1:17){
    contact_constant[i,j] = G[i,j]*vec[j]/sum(vec)
  }
}


