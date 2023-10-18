#install.packages('socialmixr')
#library(socialmixr)

data(polymod)

M <- socialmixr::contact_matrix(polymod, countries = "United Kingdom", 
               age.limits = c(0, 1, seq(5,75,5)),
               return.demography = T)$matrix

# heatmap(M, Colv = NA, Rowv = NA, scale='none',
#         xlab="Contact age", ylab="Participant age", main="Contact matrix")

## GGPLOT FORMATTING

ages <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29",
          "30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
          "70 to 74","75+")

data <- data.frame(Contact = rep(ages, each=17),
                   Participant = rep(ages, 17),
                   Var1 = expand.grid(M)) 

data$Contact <- factor(data$Contact, levels=unique(data$Contact))
data$Participant <- factor(data$Participant, levels=unique(data$Participant))

ggplot(data, aes(Contact, Participant, fill= Var1)) + 
  geom_tile() + scale_fill_gradient(low="grey99", high="red") + labs(fill = "Contacts") +
  theme(axis.text.x = element_text(angle = 90,color='black'),
        axis.text.y = element_text(color='black'),
        text = element_text(size = 14)) 

# ** REPARAMETRISING TO NEW AGE STRUCTURE **

# ** INTRINSIC CONNECTIVITY MATRIX ** 

demog2005 <- socialmixr::contact_matrix(polymod, countries = "United Kingdom", 
                    age.limits = c(0, 1, seq(5,75,5)),
                    return.demography = T)$demography #GB population structure 2005
N <- demog2005$population
G <- matrix(nrow=17, ncol=17)

for(i in 1:17){
  for(j in 1:17){
    G[i,j] = M[i,j]*sum(N)/N[j]
  }
}

#write.csv(G, "/G.csv", row.names=F)

rural_age <- read.csv("/rural_age.csv")
rural_age$Age <- as.character(rural_age$Age)
rural_age$Age <- factor(rural_age$Age, levels=unique(rural_age$Age))

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

# heatmap(final_matrix(1,T), Colv = NA, Rowv = NA,
#         xlab="Contact age", ylab="Participant age", main="Contact matrix", scale='none')


## GGPLOT FORMAT 

final_matrix_plot <- function(i, u){
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
  plot <- data.frame(Contact = rep(ages, each=17),
                         Participant = rep(ages, 17),
                         Var1 = expand.grid(final_matrix)) 
  plot$Contact <- factor(plot$Contact, levels=unique(plot$Contact))
  plot$Participant <- factor(plot$Participant, levels=unique(plot$Participant))
  return(plot)
}

# ggplot(final_matrix_plot(1,T), aes(Contact, Participant, fill= Var1)) + 
#   geom_tile() + scale_fill_gradient(low="grey99", high="red",limits=c(0,9)) + labs(fill = "Contacts") +
#   theme(axis.text.x = element_text(angle = 90)) 


# ** MATRICES THAT EXCLUDE SCHOOL CONTACTS **

M_school <- contact_matrix(polymod, countries = "United Kingdom", 
                    age.limits = c(0, 1, seq(5,75,5)),
                    return.demography = T, filter = list(cnt_school=1))$matrix

# heatmap(M_school, Colv = NA, Rowv = NA, scale='none',
#         xlab="Contact age", ylab="Participant age", main="Contact matrix")

G_school <- matrix(nrow=17, ncol=17)

for(i in 1:17){
  for(j in 1:17){
    G_school[i,j] = M_school[i,j]*sum(N)/N[j]
  }
}

# heatmap(G_school, Colv = NA, Rowv = NA,scale='none',
#         xlab="Contact age", ylab="Participant age", main="Contact matrix")

school_closure_matrix <- function(i, u, x=1){
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
} # x the extent to which schools are closed

# heatmap(school_closure_matrix(10,F), Colv = NA, Rowv = NA,
#         xlab="Contact age", ylab="Participant age", main="Contact matrix", scale='none')
