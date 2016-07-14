## R-Script - modeling species richness
## author: Javier Lopatin
## mail:javierlopatin@gmail.com
## Manuscript: Comparing Generalized Linear Models and random forest to model vascular plant species richness using LiDAR data in a natural forest in central Chile
## last changes: 12/11/2015

library(lme4)
library(hier.part)
library(splines)
library(MASS)
library(randomForest)

##### set working directory
setwd("direction/to/your/folder")

#### Load data
dat <- read.table("Richness_model.csv", header=T, sep=",", dec=".") 
attach(dat)
summary(dat)


###############################
## Prepare Bootstrap samples
###############################

set.seed(550)

# create empty lists in which subsets can be stored
train <- list()
validation <- list()

# set the bootstrap parameters
N = length(dat[,1]) # N° of observations
B = 500             # N° of bootstrap iterations

# start loop
for(i in 1:B){
  
  # create random numbers with replacement to select samples from each group
  idx = sample(1:N, N, replace=TRUE)
  
  # select subsets of the five groups based on the random numbers
  train[[i]] <- dat[idx,]
  validation[[i]] <- dat[-idx,]
  }

#################################################
## start regression modelling with Random Forest
#################################################
# create empty lists in which the models accuracies can be stored
# Obs = observed variable
# Pred = predictors
# r2 = squared Pearson's correlation
# Nrmse = normalized root mean squared error
# imp = variable importance
# bias = bias of the model

# Total
Obs.rf<-list()
Pred.rf<-list()
r2.rf<-list()
rmse.rf<-list()
Nrmse.rf<-list()
imp.rf<-list()
bias.rf<-list()

# Tree
Obs.rf.A<-list()
Pred.rf.A<-list()
r2.rf.A<-list()
rmse.rf.A<-list()
Nrmse.rf.A<-list()
imp.rf.A<-list()
bias.rf.A <- list()

# Shrub
Obs.rf.AR<-list()
Pred.rf.AR<-list()
r2.rf.AR<-list()
rmse.rf.AR<-list()
Nrmse.rf.AR<-list()
imp.rf.AR<-list()
bias.rf.AR <- list()

# Herbs
Obs.rf.H<-list()
Pred.rf.H<-list()
r2.rf.H<-list()
rmse.rf.H<-list()
Nrmse.rf.H<-list()
imp.rf.H<-list()
bias.rf.H<-list()

# Run RF
for(i in 1:B){
  
  TRAIN <- train[[i]] 
  VALIDATION <- validation[[i]]
  len<-length( VALIDATION[,1])
  
  # store and select the observations
  obs <- VALIDATION[,2]
  Obs.rf[[i]]<-obs
  
  #select the predictors for the train and validation
  variables<-na.roughfix(TRAIN[,6:17])
  variables.t<-na.roughfix(VALIDATION[,6:17])

  # run the RF model using ntrees=500 (selected by bibliography) and mtry=7 (selected after initial tuning procidure)
  RF_total<-randomForest(TRAIN[,2] ~ ., data= variables, ntrees=500, na.action= na.roughfix,importance=F, do.trace=100, mtry=7)
  
  # predict richness values using the 
  Pred<-predict(RF_total,variables.t) 
  
  # store the model accuracies
  Pred<-Pred[1:len]
  Pred.rf[[i]]<-Pred
  r2.rf[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.rf[[i]]<-s1
  Nrmse.rf[[i]]<-(s1/(max(obs)-min(obs)))*100
  imp.rf[[i]]<-importance(RF_total, type=1)
  lm = lm(Pred ~ obs-1)
  bias.rf[[i]] <-1-coef(lm)
  
  # starting here this process is repeated for tree, shrub and herb richness
  # tree
  obs <- VALIDATION[,3]
  Obs.rf.A[[i]]<-obs
  RF_tree<-randomForest(TRAIN$Tree_richness ~ ., data= variables, ntrees=500, na.action= na.roughfix,importance=TRUE, do.trace=100, mtry=12)
  Pred<-predict(RF_tree,variables.t) 
  Pred<-Pred[1:len]
  Pred.rf.A[[i]]<-Pred
  r2.rf.A[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.rf.A[[i]]<-s1
  Nrmse.rf.A[[i]]<-(s1/(max(obs)-min(obs)))*100
  imp.rf.A[[i]]<-importance(RF_tree, type=1)
  lm = lm(Pred ~ obs-1)
  bias.rf.A[[i]] <-1-coef(lm)
  
  # shrub
  obs <- VALIDATION[,4]
  Obs.rf.AR[[i]]<-obs
  RF_shrub<-randomForest(TRAIN$Shrub_richness ~ ., data= variables, ntrees=500, na.action= na.roughfix,importance=TRUE, do.trace=100, mtry=7)
  Pred<-predict(RF_shrub,variables.t) 
  Pred<-Pred[1:len]
  Pred.rf.AR[[i]]<-Pred
  r2.rf.AR[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.rf.AR[[i]]<-s1
  Nrmse.rf.AR[[i]]<-(s1/(max(obs)-min(obs)))*100
  imp.rf.AR[[i]]<-importance(RF_shrub, type=1)
  lm = lm(Pred ~ obs-1)
  bias.rf.AR[[i]] <-1-coef(lm)
  
  # herb
  obs <- VALIDATION[,5]
  Obs.rf.H[[i]]<-obs
  RF_herb<-randomForest(TRAIN$Herb_richness ~ ., data= variables, ntrees=500, na.action= na.roughfix,importance=TRUE, do.trace=100, mtry=12)
  Pred<-predict(RF_herb,variables.t) 
  Pred<-Pred[1:len]
  Pred.rf.H[[i]]<-Pred
  r2.rf.H[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.rf.H[[i]]<-s1
  Nrmse.rf.H[[i]]<-(s1/(max(obs)-min(obs)))*100
  imp.rf.H[[i]]<-importance(RF_herb, type=1)
  lm = lm(Pred ~ obs-1)
  bias.rf.H[[i]] <-1-coef(lm)
}


######################################
## start regression modelling with GLM
######################################
# create empty lists in which the models accuracies can be stored
# Total
ID.nb<-list()
Obs.nb<-list()
Pred.nb<-list()
r2.nb<-list()
rmse.nb<-list()
Nrmse.nb<-list()
imp.nb<-list()
bias.nb <- list()

# Tree
ID.nb.A<-list()
Obs.nb.A<-list()
Pred.nb.A<-list()
r2.nb.A<-list()
rmse.nb.A<-list()
Nrmse.nb.A<-list()
imp.nb.A<-list()
bias.nb.A <- list()

# Shrub
ID.nb.AR<-list()
Obs.nb.AR<-list()
Pred.nb.AR<-list()
r2.nb.AR<-list()
rmse.nb.AR<-list()
Nrmse.nb.AR<-list()
imp.nb.AR<-list()
bias.nb.AR <- list()

# Herb
ID.nb.H<-list()
Obs.nb.H<-list()
Pred.nb.H<-list()
r2.nb.H<-list()
rmse.nb.H<-list()
Nrmse.nb.H<-list()
imp.nb.H<-list()
bias.nb.H <- list()

# Run GLM
for(i in 1:B){
  TRAIN <- train[[i]] 
  VALIDATION <- validation[[i]]
  len<-length( VALIDATION[,1])
  
  # store and select the observations
  ID<-VALIDATION$ID
  ID.nb [[i]]<-ID
  obs <- VALIDATION[,2]
  Obs.nb[[i]]<-obs
  
  # run the GLM using Negative Binomial family. Three variables were selected using previous tuning procidere 
  GLM_total <- glm(Total_richness ~  one_mean + DTM_1_mean + slope_1m_std, data=TRAIN,  family=negative.binomial(theta=1 , link="log"))
  
  # predict richness values 
  Pred<-stats:::predict(GLM_total, newdata=VALIDATION, type="response")
  
  # store the model accuracies
  Pred.nb[[i]]<-Pred
  r2.nb[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.nb[[i]]<-s1
  Nrmse.nb[[i]]<-(s1/(max(obs)-min(obs)))*100
  hp <- hier.part(VALIDATION$Total_richness, VALIDATION[,c("DTM_1_mean", "slope_1m_std", "norm_H_1_mean", "Asp_1m", "TWI_1m", "one_mean", "one_std", "homogeneity_1", "contrast_1", "dissimilarity_1", "entropy_1", "second_moment_1" )],  family=negative.binomial(theta=1 , link="log"))
  imp.nb[[i]]<-hp$I.perc
  lm = lm(Pred ~ obs-1)
  bias.nb[[i]] <-1-coef(lm)
  
  # starting here this process is repeated for tree, shrub and herb richness
  # tree
  ID<-VALIDATION$ID
  ID.nb.A[[i]]<-ID
  obs <- VALIDATION[,3]
  Obs.nb.A[[i]]<-obs
  GLM_tree <- glm(Tree_richness ~  one_mean + DTM_1_mean + slope_1m_std, data=TRAIN,  family=negative.binomial(theta=1 , link="log"))
  Pred<-stats:::predict(GLM_tree, newdata=VALIDATION, type="response")
  Pred.nb.A[[i]]<-Pred
  r2.nb.A[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.nb.A[[i]]<-s1
  Nrmse.nb.A[[i]]<-(s1/(max(obs)-min(obs)))*100
  hp <- hier.part(VALIDATION$Total_richness, VALIDATION[,c("DTM_1_mean", "slope_1m_std", "norm_H_1_mean", "Asp_1m", "TWI_1m", "one_mean", "one_std", "homogeneity_1", "contrast_1", "dissimilarity_1", "entropy_1", "second_moment_1" )],  family=negative.binomial(theta=1 , link="log"))
  imp.nb.A[[i]]<-hp$I.perc
  lm = lm(Pred ~ obs-1)
  bias.nb.A[[i]] <-1-coef(lm)
  
  # shrub
  ID<-VALIDATION$ID
  ID.nb.AR[[i]]<-ID
  obs <- VALIDATION[,4]
  Obs.nb.AR[[i]]<-obs
  GLM_shrub <- glm(Shrub_richness ~  one_mean + DTM_1_mean + TWI_1m, data=TRAIN,  family=negative.binomial(theta=1 , link="log"))
  Pred<-stats:::predict(GLM_shrub, newdata=VALIDATION, type="response")
  Pred.nb.AR[[i]]<-Pred
  r2.nb.AR[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.nb.AR[[i]]<-s1
  Nrmse.nb.AR[[i]]<-(s1/(max(obs)-min(obs)))*100
  hp <- hier.part(VALIDATION$Total_richness, VALIDATION[,c("DTM_1_mean", "slope_1m_std", "norm_H_1_mean", "Asp_1m", "TWI_1m", "one_mean", "one_std", "homogeneity_1", "contrast_1", "dissimilarity_1", "entropy_1", "second_moment_1" )],  family=negative.binomial(theta=1 , link="log"))
  imp.nb.AR[[i]]<-hp$I.perc
  lm = lm(Pred ~ obs-1)
  bias.nb.AR[[i]] <-1-coef(lm)
  
  # herb
  ID<-VALIDATION$ID
  ID.nb.H[[i]]<-ID
  obs <- VALIDATION[,5]
  Obs.nb.H[[i]]<-obs
  GLM_herb <- glm(Herb_richness ~  one_mean + DTM_1_mean + slope_1m_std, data=TRAIN,  family=negative.binomial(theta=1 , link="log"))
  Pred<-stats:::predict(GLM_herb, newdata=VALIDATION, type="response")
  Pred.nb.H[[i]]<-Pred
  r2.nb.H[[i]]<-(cor(Pred, obs, method="pearson"))^2
  s1<-sqrt(mean((obs-Pred)^2))
  rmse.nb.H[[i]]<-s1
  Nrmse.nb.H[[i]]<-(s1/(max(obs)-min(obs)))*100
  hp <- hier.part(VALIDATION$Total_richness, VALIDATION[,c("DTM_1_mean", "slope_1m_std", "norm_H_1_mean", "Asp_1m", "TWI_1m", "one_mean", "one_std", "homogeneity_1", "contrast_1", "dissimilarity_1", "entropy_1", "second_moment_1" )],  family=negative.binomial(theta=1 , link="log"))
  imp.nb.H[[i]]<-hp$I.perc
  lm = lm(Pred ~ obs-1)
  bias.nb.H[[i]] <- 1-coef(lm)
}

# export the variable importances
# RF
write.table(imp.rf, file="importancia.rf.csv")
write.table(imp.rf.A, file="importancia.rf.A.csv")
write.table(imp.rf.AR, file="importancia.rf.AR.csv")
write.table(imp.rf.H, file="importancia.rf.H.csv")
# GLM
write.table(imp.nb, file="importancia.nb.csv")
write.table(imp.nb.A, file="importancia.nb.A.csv")
write.table(imp.nb.AR, file="importancia.nb.AR.csv")
write.table(imp.nb.H, file="importancia.nb.H.csv")

# merge all accuracies together
BOOTS_ACC <- data.frame(unlist(r2.rf), unlist(r2.rf.A), unlist(r2.rf.AR), unlist(r2.rf.H), 
                        unlist(r2.nb), unlist(r2.nb.A), unlist(r2.nb.AR), unlist(r2.nb.H),
                        unlist(Nrmse.rf), unlist(Nrmse.rf.A), unlist(Nrmse.rf.AR), unlist(Nrmse.rf.H),
                        unlist(Nrmse.nb), unlist(Nrmse.nb.A), unlist(Nrmse.nb.AR), unlist(Nrmse.nb.H),
                        unlist(bias.rf), unlist(bias.rf.A), unlist(bias.rf.AR), unlist(bias.rf.H), 
                        unlist(bias.nb), unlist(bias.nb.A), unlist(bias.nb.AR), unlist(bias.nb.H))
colnames(BOOTS_ACC) <- c("r2.rf.Total","r2.rf.Tree", "r2.rf.Shrub", "r2.rf.Herb",
                         "r2.nb.Total","r2.nb.Tree", "r2.nb.Shrub", "r2.nb.Herb",
                         "Nrmse.rf.Total","Nrmse.rf.Tree", "Nrmse.rf.Shrub", "Nrmse.rf.Herb",
                         "Nrmse.nb.Total","Nrmse.nb.Tree", "Nrmse.nb.Shrub", "Nrmse.nb.Herb",
                         "bias.rf.Total","bias.rf.Tree", "bias.rf.Shrub", "bias.rf.Herb",
                         "bias.nb.Total","bias.nb.Tree", "bias.nb.Shrub", "bias.nb.Herb")
# export the results
write.table(BOOTS_ACC, file="BOOTS_ACC.csv")

## Residuals
# RF
res.rf<- unlist(Pred.rf)-unlist(Obs.rf)
res.rf.A<- unlist(Pred.rf.A)-unlist(Obs.rf.A)
res.rf.AR<- unlist(Pred.rf.AR)-unlist(Obs.rf.AR)
res.rf.H<- unlist(Pred.rf.H)-unlist(Obs.rf.H)
#GLM
res.nb<- unlist(Pred.nb)-unlist(Obs.nb)
res.nb.A<- unlist(Pred.nb.A)-unlist(Obs.nb.A)
res.nb.AR<- unlist(Pred.nb.AR)-unlist(Obs.nb.AR)
res.nb.H<- unlist(Pred.nb.H)-unlist(Obs.nb.H)

# median accuracies of the models
MED <- list()
medians <- for (i in 1:length(BOOTS_ACC[1,])){
  a <- median(BOOTS_ACC[,i])
  a <- round(a, 2)
  MED[[i]] <- a
}
MED <- as.data.frame(unlist(MED), names(BOOTS_ACC))
MED


save.image("Richness.RData")

