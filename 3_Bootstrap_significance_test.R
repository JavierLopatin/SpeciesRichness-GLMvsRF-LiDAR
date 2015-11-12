## R-Script - Bootstrap test for sigfificative differences
## author: Klara Dolos
## mail: klara.dolos@kit.de
## Manuscript: Comparing Generalized Linear Models and random forest to model vascular plant species richness using LiDAR data in a natural forest in central Chile
## last changes: 12/11/2015

library(randomForest)
library(boot)
library(MASS)

##### set working directory
setwd("direction/to/your/folder")


#### Set the bootstrap test for sigfificative differences
to_boot_all <- function(data, i){
  reponse_names <- c("Tolat_richness", "Shrub_richness", "Tree_richness", "Herb_richness")
  diff_r2 <- numeric()
  diff_bias <- numeric()
  diff_rmse <- numeric()
  for ( i in 1:length(reponse_names)){
    data$response <- data[,reponse_names[i]]
    # best GLM model
    m1 <- glm(response ~ DTM_1_mean + slope_1m_std + norm_H_1_mean, data=data, family=negative.binomial(theta=1 , link="log"))    
    # best RF model
    m2 <- randomForest(response ~ DTM_1_mean + slope_1m_std + norm_H_1_mean + norm_H_1_mean + Asp_1m + TWI_1m + TWI_1m + one_mean + one_std + homogeneity_1 + contrast_1 + dissimilarity_1 + entropy_1 + second_moment_1, data=data,  ntrees=500, na.action= na.roughfix,importance=F, do.trace=100, mtry=7)
    
    ### compute the differences between GLM and RF
    ### r2 of GLM should be large. So, if r2(m1) - r2(m2) is positive, if GLM is better.
    diff_r2[i] <- cor(data$response, predict(m1))^2 - cor(data$response, predict(m2))^2
    
    ### RMSE of GLM should be small. So, if rmse(m2) - rmse(m1) is positive, if GLM is better. Caution, I change the order!
    diff_rmse[i] <- sqrt(mean((data$response-predict(m2, type="response"))^2)) - sqrt(mean((data$response-predict(m1, type="response"))^2))
    
    ### Bias of GLM should be small. So, if bias(m2) - bias(m1) is positive, if GLM is better. Caution, I change the order!
    diff_bias[i] <- (1 - lm(data$response ~ predict(m2, type="response") -1)$coef)- (1 - lm(data$response ~ predict(m1, type="response") -1)$coef)  
    }
    # store the differences 
   c(diff_r2, diff_rmse, diff_bias)
}


### Execute ####

# load data
dat <- read.table("Richness_model.csv", header=T, sep=",", dec=".")   ### Dataset: Javier
head(dat)

# apply boot test and prepare the results
boot.result <- boot(data=dat, statistic=to_boot_all, R=500,  stype = "i")  ### Takes about 2-4 min
boot.result
boot.values <- boot.result$t
colnames(boot.values) <-  paste(rep(c("Tolat_richness", "Shrub_richness", "Tree_richness", "Herb_richness"), 3), rep(c("r2", "RMSE", "bias"), each=4), sep="_")
par(mfrow=c(3,4), mar=c(2,3,3,1))
for(i in 1:ncol(boot.values)){
  hist(boot.values[,i], main=colnames(boot.values)[i], col="grey", border="white", xlab="", ylab="")
  abline(v=quantile(boot.values[,i], probs=c(0.05, 0.95)), col=c("blue", "green"))
  abline(v=0, col=c("black"))
  box()
}
### The black "zero-line" needs to be left of the blue "alpha-line". The green line is just the upper quantile.

### Range of the values
apply(boot.values, 2, function(x) range(x))

### Lower quantile needs to be larger than zero.
test.result <- (t(apply(boot.values, 2, function(x) quantile(x, probs = c(0.01, 0.05)) > 0))); test.result
# save results
write.table(test.result, "bootstrap_test_result.txt")
