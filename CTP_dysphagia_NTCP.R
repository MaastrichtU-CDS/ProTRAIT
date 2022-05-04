#######################################################
#####install all the packages and libraries I need#####
#######################################################
library(dplyr)
library("readxl")
library(openxlsx)
library("data.table")
library(tidyverse)
library(rms)
library(Hmisc)
library(knitr)
library(broom)
library(pander)
library(ggbeeswarm)
library(gridExtra)
library(grid)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(viridis)
library(mlbench)
library(pROC)
library(ModelGood)
library(rms)
library(dplyr)


#import the dataset wwith the DVH data and toxicities
dvhdataunifiednew<- read_excel("path\\naldataset.xlsx")



dvhdataunifiednew$dysphagia_baseline<-ifelse(dvhdataunifiednew$dysphagia_baseline <=1, 0,
                                         ifelse(dvhdataunifiednew$dysphagia_baseline == 2 , 0.9382,
                                                ifelse(dvhdataunifiednew$dysphagia_baseline >= 3 , 1.2900 , dvhdataunifiednew$dysphagia_baseline)))



dvhdataunifiednew$Tumorlokatie <-ifelse(dvhdataunifiednew$Tumorlokatie  == "mondholte", 0,
                                   ifelse(dvhdataunifiednew$Tumorlokatie == "farynxtumor" , -0.6281,
                                          ifelse(dvhdataunifiednew$Tumorlokatie == "larynxtumor" , -0.7711 ,dvhdataunifiednew$Tumorlokatie)))




dvhdataunifiednew$dysphagia_baseline<-as.double(dvhdataunifiednew$dysphagia_baseline)
dvhdataunifiednew$Tumorlokatie<-as.double(dvhdataunifiednew$Tumorlokatie)
dvhdataunifiednew$`Dmean_PCM_inferior_(Gy)`<-as.double(dvhdataunifiednew$`Dmean_PCM_inferior_(Gy)`)
dvhdataunifiednew$`Dmean_PCM_superior_(Gy)`<-as.double(dvhdataunifiednew$`Dmean_PCM_superior_(Gy)`)
dvhdataunifiednew$`Dmean_PCM_medius_(Gy)`<-as.double(dvhdataunifiednew$`Dmean_PCM_medius_(Gy)`)
dvhdataunifiednew$`Dmean_mondholte_(Gy)`<-as.double(dvhdataunifiednew$`Dmean_mondholte_(Gy)`)




#lets calculate the linear predictor for the NTCP model LIPP2 for dypshagia more than 2nd grade 6 months

#lets calculate the linear predictor for the NTCP model LIPP2 for dypshagia more than 2nd grade 6 months
dvhdataunifiednew$linear_factor <- (-4.0536
                                       +0.030*(dvhdataunifiednew$`Dmean_mondholte_(Gy)`)
                                       +0.0236*(dvhdataunifiednew$`Dmean_PCM_superior_(Gy)`)
                                       +0.0095*(dvhdataunifiednew$`Dmean_PCM_medius_(Gy)`)
                                       +0.0133*(dvhdataunifiednew$`Dmean_PCM_inferior_(Gy)`)
                                       +offset(dvhdataunifiednew$dysphagia_baseline)
                                       +offset(dvhdataunifiednew$Tumorlokatie))


#convert the numeric values to character that will be transformed to a factor
table(dvhdataunifiednew$dysphagia_sixmonths)
dvhdataunifiednew$dysphagia_sixmonths<- ifelse(dvhdataunifiednew$dysphagia_sixmonths >=2,1, 0)


#####################################################################################
######################create the four different models of the CTP####################
#####################################################################################

# original model (intercept = 0 and beta = 1)
m0 <- glm(dysphagia_sixmonths~ -1 + offset(linear_factors), data = dvhdataunified, family = binomial)

# Recalibration in the Large (estimate intercept, beta = 1)
m1 <- glm(dysphagia_sixmonths ~ 1 +offset(linear_factors), data = dvhdataunified, family = binomial)

# Logisitic Calibration (estimate both)
m2 <- glm(dysphagia_sixmonths ~1 +  linear_factors,  data = dvhdataunified, family = binomial)

# model revision###
formula_4 <- as.formula(
  "dysphagia_sixmonths ~  1+ `Dmean_mondholte_(Gy)` + `Dmean_PCM_superior_(Gy)` +  dysphagia_baseline "
)

m3 <- glm(formula_4, data = dvhdataunified, family = binomial)


#####################################################################################
######################produce the AUC curves for the m1 and m4#######################
#####################################################################################

roc(dysphagia_sixmonths~m0$fitted.values, data = dvhdataunified, plot = TRUE,lwd = 2, legacy.axes = FALSE, print.auc=TRUE, ci=TRUE, main = "ROC curves of the original and revised NTCP model", col= "blue")
roc(dysphagia_sixmonths~m3$fitted.values, data = dvhdataunified, plot = TRUE,lwd = 2, add= TRUE,legacy.axes = FALSE,ci=TRUE, print.auc=TRUE, col = "red",print.auc.y=0.45)




#####################################################################################
######################calculate the sensitivity and specificity of the models########
#####################################################################################



threshold=0.5
predicted_values<-ifelse(predict(m1,type="response")>threshold,1,0)
actual_values<-m0$y
conf_matrix<-table(actual_values, predicted_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)


threshold=0.5
predicted_values<-ifelse(predict(m3,type="response")>threshold,1,0)
actual_values<-m4$y
conf_matrix<-table(actual_values, predicted_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)

threshold=0.5
predicted_values<-ifelse(predict(m2,type="response")>threshold,1,0)
actual_values<-m2$y
conf_matrix<-table(actual_values, predicted_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)


threshold=0.5
predicted_values<-ifelse(predict(m1,type="response")>threshold,1,0)
actual_values<-m3$y
conf_matrix<-table(actual_values, predicted_values)
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)




#####################################################################################
######################produce the  calibration curves of the different models########
#####################################################################################

calPlot2(list(m0,m1,m2,m3),data=dvhdataunified,legend=FALSE,xlab = "Predicted probability" ,ylab = "Observed probability",cex.axis=0.2)
# plotl()
legend(x = 0.8, y = 0.4, # Coordinates          # Position
       legend = c("original model", "Re-calibration in the large", "Logistic Recalibration","Model revision"),  
       lty = c(1),
       bty = "n",
       inset = c(-0.45, 0),# Line types
       cex = 0.8, # Change legend size
       
       col = c("black","red","green", "lightblue"),         
       lwd = 2,xpd = TRUE,horiz = FALSE)




############################################################################################
######################produce the  calibration curves for the supplementary material########
############################################################################################

val.prob(fitted(m0),dvhdataunified$dysphagia_sixmonths)
val.prob(fitted(m1),dvhdataunified$dysphagia_sixmonths)
val.prob(fitted(m2),dvhdataunified$dysphagia_sixmonths)
val.prob(fitted(m3),dvhdataunified$dysphagia_sixmonths)


#####################################################################
########################Closed testing procedure#####################
#####################################################################


dvhdataunifiednew<- select(dvhdataunified,
                         'Dmean_mondholte_(Gy)',
                         'Dmean_PCM_superior_(Gy)',
                         dysphagia_sixmonths,
                         dysphagia_baseline)


coefs<-as.vector(c(-3.303, 0.024, 0.024, 0.967))
X<-as.matrix(dvhdataunified$Dmean_mondholte,dvhdataunified$Dmean_PCM_superior,dvhdataunified$dysphagia_baseline)
y<-as.vector(dvhdataunified$dysphagia_sixmonths)


ClosedTest <- function(coefs, X, y){
  require(rms)
  if(class(X)=="data.frame"){
    X <- data.matrix(X)
  }
  if(ncol(X)!=(length(coefs)-1)){
    stop("Number of predictors not equal to the number of coefficients")
  }
  n_coefs <- length(coefs)
  lp_old <- X %*% as.matrix(coefs[2:n_coefs])
  # Calculate updated model intercept
  intercept <- lrm.fit(y = y, offset = lp_old)$coefficients
  coefs_int <- c(intercept , coefs[2:n_coefs])
  # Calculate coefficients after recalibration
  recal <- lrm.fit(x = lp_old, y = y)$coefficients
  coefs_recal <- c(recal[1], recal[2] * coefs[2:n_coefs])
  # Calculate coefficients after model revision
  coefs_refit <- lrm.fit(x = X, y = y)$coefficients
  # Calculate the log-likelihood of the different models
  lp <- cbind(1, X) %*% coefs
  ll_original <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_int
  ll_intercept <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_recal
  ll_recalibration <- sum(y * lp - log(1 + exp(lp)))
  lp <- cbind(1, X) %*% coefs_refit
  ll_revision <- sum(y * lp - log(1 + exp(lp)))
  # Calculate difference in log-likelihood for testing of the models
  dev_original <- -2 * ll_original + 2 * ll_revision
  dev_intercept <- -2 * ll_intercept + 2 * ll_revision
  dev_recalibration <- -2 * ll_recalibration + 2 * ll_revision
  # See if difference in model fit was significant
  test1 <- (1 - pchisq(dev_original, ncol(X) + 1)) < 0.05
  test2 <- (1 - pchisq(dev_intercept, ncol(X))) < 0.05
  test3 <- (1 - pchisq(dev_recalibration, ncol(X) - 1)) < 0.05
  # See which model is chosen, index_test indicates the chosen model
  # 1. Original model
  # 2. Model with updated intercept
  # 3. Recalibrated model
  # 4. Revised model
  test_original <- 1 * (!test1)
  test_intercept <- 2 * ((test1)&(!test2))
  test_recalibration <- 3 * ((test1)&(test2)&(!test3))
  test_revision <- 4 * ((test1)&(test2)&(test3))
  index_test <- (test_original + test_intercept + test_recalibration +
                   test_revision)
  coefs_result <- rbind(coefs, coefs_int, coefs_recal, coefs_refit)
  
  # Output of the function
  new_coefs <- coefs_result[index_test, ]
  model <- c("Original Model", "Model with updated intercept",
             "Recalibrated model", "Model Revision")[index_test]
  cat("Method chosen by closed test procedure:\n", model, "\n",
      "Resulting coefficients:\n", new_coefs, "\n")
  res <- list(model = model, coefs = coefs_result)
  return(res)
} 