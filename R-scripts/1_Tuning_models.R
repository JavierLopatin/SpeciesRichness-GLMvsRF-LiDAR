## R-Script - Model tuning
## author: Javier Lopatin & Klara Dolos
## mail: javierlopatin@gmail.com; klara.dolos@kit.de
## Manuscript: Comparing Generalized Linear Models and random forest to model vascular plant species richness using LiDAR data in a natural forest in central Chile
## last changes: 12/11/2015

library(lme4)
library(hier.part)
library(splines)
library(MASS)
library(caret)

##### set working directory
setwd("direction/to/your/folder")

### Load data ####
dat <- read.table("Richness_model.csv", header=T, sep=",", dec=".")   
head(dat)
summary(dat)

###############################################################
### Tuning the Random Forest models using the caret package ###
###############################################################

# select the predictors
predictors <- dat[,6:17]

# define parameter tuning method
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, returnResamp ="all")

# tuning for total richness
tuning_total <- train(predictors, dat$Total_richness, method = "rf", trControl = fitControl)
tuning_total

# tuning for tree richness
tuning_tree <- train(predictors, dat$Tree_richness, method = "rf", trControl = fitControl)
tuning_tree

# tuning for shrub richness
tuning_shrub <- train(predictors, dat$Shrub_richness, method = "rf", trControl = fitControl)
tuning_shrub

# tuning for herb richness
tuning_herb <- train(predictors, dat$Herb_richness, method = "rf", trControl = fitControl)
tuning_herb

#######################
### Tuning the GLMs ###
#######################

### Explore data using the total richness dataset ###
x11(width=10, height=10)
par(mfrow=c(3,3))
plot(dat$Total_richness ~ dat$one_mean)
plot(dat$Total_richness ~ dat$DTM_1_mean)
plot(dat$Total_richness ~ dat$slope_1m_std)
plot(dat$Total_richness ~ dat$Asp_1m)
plot(dat$Total_richness ~ dat$DTM)
plot(dat$Total_richness ~ dat$homogeneity_1)

### pairs(dat)
### Do not use variables correlated with r > 0.7 at the same time. This is only a rule of thumb. 
round(cor(dat),2)   

### Find model family for Species Number ###
hist(dat$Total_richness)

### Simple (non-adequate) linear regression.
m1 <- glm(Total_richness ~ one_mean, data=dat, family=gaussian(link="identity"))
par(mfrow=c(2,2))
plot(m1) 
summary(m1)
par(mfrow=c(1,2))
plot(resid(m1))
hist(resid(m1))

### At least log-link transformation, which sometimes helps to model count data.
m2 <- glm(Total_richness ~ one_mean, data=dat, family=gaussian(link="log"))
par(mfrow=c(2,2))
plot(m2)  
summary(m2)
par(mfrow=c(1,2))
plot(resid(m2))
hist(resid(m2))

### Test Poisson family with log-link
m3 <- glm(Total_richness ~ one_mean, data=dat, family=poisson(link="log"))
par(mfrow=c(2,2))
plot(m3)   
summary(m3)
par(mfrow=c(1,2))
plot(resid(m3))
hist(resid(m3))

### Test Quasi-Poisson family with log-link
m4 <- glm(Total_richness ~ one_mean, data=dat, family=quasipoisson(link="log"))
par(mfrow=c(2,2))
plot(m4)   
summary(m4)
par(mfrow=c(1,2))
plot(resid(m4))
hist(resid(m4))

### Test Negative binomial family with theta=1. QQplot: Appears to be the best one.
m5 <- glm(Total_richness ~ one_mean, data=dat, family=negative.binomial(theta = 1, link="log"))
par(mfrow=c(2,2))
plot(m5)   
summary(m5)
par(mfrow=c(1,2))
plot(resid(m5))
hist(resid(m5))

# see the assumed error distribution using different theta
m5.nb <- glm.nb(Total_richness ~ one_mean, data=dat)
m5.nb$theta
theta.est <- round(m5.nb$theta,0)
x <- 1:40
par(mfrow=c(1,1))
plot(x, dnbinom(x, size=1, mu=mean(dat$Total_richness), log = FALSE), ylim=c(0,0.1), type="l", main="Assumed error distribution", ylab="Density", xlab="Species number")
lines(x, dnorm(x, mean=mean(dat$Total_richness), sd=sd(dat$Total_richness)), lty=2)
lines(x,dnbinom(x, size=theta.est, mu=mean(dat$Total_richness), log = FALSE), lty=3)
legend("topright", legend=c("Normal distribution", "Negative binomial (theta=1)", paste("Negative binomial (theta=", theta.est, ")", sep="")), lty=c(2,1,3))


### Calculate bias of prediction
lm(predict(m1, type="link") ~ dat$Total_richness -1)$coef
lm(predict(m5, type="link") ~ log(dat$Total_richness) -1)$coef  ### m5 does not have a bias at all at the link scale, which is relevant for the model algorithm. 

lm(predict(m5, type="response") ~ (dat$Total_richness) -1)$coef ### Back transformed: This bias is relavant for predictions.
lm(predict(m1, type="response") ~ (dat$Total_richness) -1)$coef 

### Variable selection using hirarquical partitioning ###
names(dat)
hp <- hier.part(dat$Total_richness, dat[,c("DTM_1_mean", "slope_1m_std", "norm_H_1_mean", "Asp_1m", "TWI_1m", "one_mean", "one_std", "homogeneity_1", "contrast_1", "dissimilarity_1", "entropy_1", "second_moment_1" )],  family=negative.binomial(theta=1 , link="log"))
hp$I.perc
o <-  order(hp$I.perc$I, decreasing=T)
rownames(hp$I.perc)[o]

COR<-data.frame(dat$one_mean, dat$DTM_1_mean, dat$slope_1m_std, dat$TWI_1m, dat$entropy_1, dat$second_moment_1)
COR<-cor(COR)
COR

### Variable selection as usual
formulas <- expand.grid(DTM_1_mean =c(0,1), one_mean=c(0,1), slope_1m_std=c(0,1), norm_H_1_mean=c(0,1), Asp_1m=c(0,1) , TWI_1m=c(0,1) , TWI_1m=c(0,1), homogeneity_1=c(0,1) , contrast_1=c(0,1) , dissimilarity_1=c(0,1), entropy_1=c(0,1), second_moment_1=c(0,1))
formulas$form <- apply(formulas, 1, function (x) paste("Total_richness ~ ",paste(names(x)[x==1], collapse=" + ") , sep=""))

models<- list()
models[[1]] <- glm(Total_richness ~ 1, data=dat, family=negative.binomial(theta = 1, link="log"))
for(i in 2:nrow(formulas)){
  models[[i]] <- glm(as.formula(formulas$form[i]), data=dat, family=negative.binomial(theta = 1, link="log"))
}
save(models, file="models.RData")

formulas$AIC <- unlist(lapply(models, function(x) AIC(x)))

formulas$expl.dev <- unlist(lapply(models, function(x) (summary(x)$null.deviance - summary(x)$deviance) /summary(x)$null.deviance))

formulas$o.dev <- order(formulas$expl.dev, decreasing=T)
formulas$o.aic <- order(formulas$AIC, decreasing=T)
formulas[o.dev,c("form", "AIC", "expl.dev")][1:10,]

plot(formulas$AIC[formulas$o.aic])
plot(formulas$AIC[formulas$o.aic][1:100])

plot(formulas$expl.dev[formulas$o.dev][1:100])

max(formulas$expl.dev)
max(formulas$AIC)

best100 <- formulas[formulas$AIC < (max(formulas$AIC)+10) & formulas$o.dev < 100,]
nrow(best100)
plot(best100$expl.dev, best100$AIC)

sort(rowSums(best100[,1:12]), decreasing=T)
sort(colSums(best100[,1:12]), decreasing=T)   

best10 <- formulas[formulas$AIC < (max(formulas$AIC)+10) & formulas$o.dev < 10,]
nrow(best10)
plot(best10$expl.dev, best10$AIC)


best10 <- best10[order(best10$AIC, decreasing=T),]
best10

sort(rowSums(best10[,1:12]), decreasing=T)
sort(colSums(best10[,1:12]), decreasing=T)   

### Best model is the one with all variables. However, the singly parameters are not significant. Better way of model selection is thus to drop non-significant patameters.
summary(models[[as.numeric(rownames(best10))[1]]])

### Seems that the function "step" does a good job in this case. 
m <- step(models[[as.numeric(rownames(best10))[1]]])
summary(m)

par(mfrow=c(2,2))
plot(m)  

### Include quadratic terms for each variable and all of them
m3 <- update(m, . ~ . + I(one_mean^2) + I(DTM_1_mean^2) + I(slope_1m_std^2)+I(Asp_1m^2)+I(homogeneity_1^2))
summary(m3)   ### Quadratic terms do not seem to be required.
AIC(m)
AIC(m3)
m3 <- step(m3)
summary(m3)
anova(m3, m, test="LRT")  ### Do not differ significantly. Stay with the simpler model.

par(mfrow=c(2,2))
plot(predict(m, type="response") ~ dat$one_mean)
plot(predict(m, type="response") ~ dat$DTM_1_mean)
plot(predict(m, type="response") ~ dat$slope_1m_std)
plot(predict(m, type="response") ~ dat$Total_richness)
abline(0,1)

### Negative binomial model, natural splines ####
ms <- glm(Total_richness ~ ns(one_mean,2) + ns(DTM_1_mean,2) + ns(slope_1m_std,2)+ ns(Asp_1m,2)+ ns(homogeneity_1,2), data=dat, family=negative.binomial(theta=1 , link="log"))
summary(ms)
expl.dev <- round((summary(m)$null.deviance - summary(m)$deviance) /summary(m)$null.deviance, 2); expl.dev

par(mfrow=c(2,2))
plot(ms)

AIC(m)
AIC(ms)  ### AIC with natural splines is not better. Stay with simpler models.
anova(m, ms, test="LRT")  ### No significant difference among models. -> Stay witht he simpler one.

### Also response curves do not show the need for a more complex model.

par(mfrow=c(2,2))
plot(predict(ms, type="response") ~ dat$one_mean)
plot(predict(ms, type="response") ~ dat$DTM_1_mean)
plot(predict(ms, type="response") ~ dat$slope_1m_std)
plot(predict(ms, type="response") ~ dat$Total_richness)
abline(0,1)


### Conclusion: GLM (neg.binom with theta=1, log-link) and just the linear terms for the selected variables is best model for total richness. 
