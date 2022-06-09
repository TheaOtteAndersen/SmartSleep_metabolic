
###### Analysis of associations between metabolic biomarkers and smartphone activity ######


#This script contains analyses of the multiply imputed data sets from the SmartSleep project.
#The multiple imputation results are combined using Rubin's rules.

library(dplyr)
library(matrixStats)
library(MASS)
library(mitml)
library(dataReporter)
library(stringr)
library(numbers)
library(gamlss)
library(mice)
library(miceadds)
library(Publish)
library(ggplot2)
library(ggtern)

expit = function(x) exp(x)/(1+exp(x))

estimate.pooler <- function(coef,sd){
  n_row <- nrow(coef)
  n_col <- ncol(coef)
  
  coefs <- rowMeans(coef)
  sds <- sqrt(rowMeans(sd^2)+(1+1/n_col)*rowSums((coef-coefs)^2)/(n_col-1))
  df <- data.frame("estimate"=coefs,"sd"=sds,"lower.CI"=coefs-1.96*sds,"upper.CI"=coefs+1.96*sds)
  return(df)
}

#Reading in the data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")

## load tracking data 
subject_tracking_six_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_clusters.csv")
subject_tracking_four_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_four_clusters.csv")

## Collecting the two clusterings in one file

subject_tracking_clusters <- left_join(subject_tracking_six_clusters,subject_tracking_four_clusters[,c("userid","cluster","cluster1prob","cluster2prob","cluster3prob","cluster4prob","description","state0prob","state1prob","state2prob","state3prob")],by="userid")
subject_tracking_clusters <- rename(subject_tracking_clusters,cluster1prob=cluster1prob.x,cluster2prob=cluster2prob.x,cluster3prob=cluster3prob.x,cluster4prob=cluster4prob.x,
                                    state0prob=state0prob.x,state1prob=state1prob.x,state2prob=state2prob.x,state3prob=state3prob.x,cluster=cluster.x,description=description.x)

## load baseline data
base_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv")
base_data$mobileUseNight <- factor(base_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
base_data$mobileUseBeforeSleep <- factor(base_data$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))

## load followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")
CSS$mobileUseNight <- factor(CSS$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
CSS$mobileUseBeforeSleep <- factor(CSS$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))

## load population sample
pop_data <-read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv")
pop_data$mobileUseNight <- factor(pop_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
pop_data$mobileUseBeforeSleep <- factor(pop_data$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))

## load clinical data (survey and clinical data)
clin_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample/imp_clinical.csv")
clin_clinical <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/SmartSleep Clinical/Data/Rådata/SmartSleepClinicalData.csv")

clin_data$mobileUseNight <- factor(clin_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
clin_data$mobileUseBeforeSleep <- factor(clin_data$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))


# --------------------------------------------------------------------------- ##

#Baseline data with self-reports

## if no mobile phone = NA
base_data$pmpuScale[base_data$mobilephone=="No mobile phone"] <- NA

## tjek risk profiles
publish(univariateTable( ~ mobileUseBeforeSleep,data=base_data, column.percent=TRUE))
base_data$mobileUseBeforeSleep <- factor(base_data$mobileUseBeforeSleep, levels = c("Never", "Every month or less", "Once a week", "2-4 times per week", "5-7 times per week"))
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))

## bmi
base_data$bmi30 <- (base_data$bmi>=30)
base_data$bmi25 <- (base_data$bmi>=25)

#save(base_data,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/base_data.RData")

base_data_mids <- as.mids(base_data,.imp="imputation")


# --------------------------------------------------------------------------- ##
#Followup sample - using quartile levels from baseline sample

table(CSS$mobileUseBeforeSleep, useNA="always")
table(CSS$mobileUseNight, useNA="always")

## merge survey and tracking data 
CSS_track <- inner_join(CSS,subject_tracking_clusters,by="userid")
CSS_track$sample_weights <- as.numeric(CSS_track$sample_weights)

#save(CSS_track,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/CSS_track.RData")

CSS_track_mids<-as.mids(CSS_track,.imp="imputation",.id="userid")

# --------------------------------------------------------------------------- ##
#Merging base and followup

#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

bmi_followup <- rename(inner_join(CSS,base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

## difference mellem follow-up og baseline
bmi_followup$difference <- bmi_followup$bmi.fu-bmi_followup$bmi.base
mean(bmi_followup$difference[!is.na(bmi_followup$difference)])

bmi_followup$basebmi25=(bmi_followup$bmi.base>=25)
bmi_followup$basebmi30=(bmi_followup$bmi.base>=30)

bmi_followup$bmi25change <- as.numeric((bmi_followup$bmi.fu>=25)!=(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeUp <- as.numeric((bmi_followup$bmi.fu>=25)>(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeDown <- as.numeric((bmi_followup$bmi.fu>=25)<(bmi_followup$bmi.base>=25))

bmi_followup$bmi30change <- as.numeric((bmi_followup$bmi.fu>=30)!=(bmi_followup$bmi.base>=30))
bmi_followup$bmi30changeUp <- as.numeric((bmi_followup$bmi.fu>=30)>(bmi_followup$bmi.base>=30))
bmi_followup$bmi30changeDown <- as.numeric((bmi_followup$bmi.fu>=30)<(bmi_followup$bmi.base>=30))

bmi_followup$sample_weights <- bmi_followup$sample_weights.y
bmi_followup$followup_time <- (as.Date(str_c(substr(bmi_followup$EndDate,7,10),"-",substr(bmi_followup$EndDate,4,5),"-",substr(bmi_followup$EndDate,1,2)))-as.Date(str_c(substr(bmi_followup$responseDate.y,7,10),"-",substr(bmi_followup$responseDate.y,4,5),"-",substr(bmi_followup$responseDate.y,1,2))))/365.25

bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

# --------------------------------------------------------------------------- ##
#Population sample

##
table(pop_data$mobileUseNight)
table(pop_data$mobileUseBeforeSleep)

## merge tracking and survey data for population sample
pop_track <- inner_join(pop_data,subject_tracking_clusters,by="userid")
pop_track$sample_weights<-as.numeric(pop_track$sample_weights)

## Tracking clusters som én numerisk variabel
pop_track$track_severity <- (pop_track$cluster %in% c("Cluster 1"))*1+(pop_track$cluster %in% c("Cluster 2","Cluster 3"))*2+(pop_track$cluster %in% c("Cluster 5","Cluster 6"))*3+(pop_track$cluster %in% c("Cluster 4"))*4

#save(pop_track,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/pop_track.RData")

pop_track_mids<-as.mids(pop_track,.imp="imputation",.id="userid")

## Clinical Data

## merge survey and clinical data
clinical_sample <- rename(inner_join(clin_data,rename(clin_clinical,PNR=cpr),by="PNR"),bmi.self=bmi.x , bmi.clinical=bmi.y)
## merge with tracking data
clinical_sample <- inner_join(clinical_sample,subject_tracking_clusters,by="userid")



#### -------------------------------- ####

#### Write out data files for gamlss bootstrap

#### -------------------------------- ####


# --------------------------------------------------------------------------- #
boot_path <- "S:/SUND-IFSV-SmartSleep/Christoffer/gamlssBootstrap/"
N_imp <- 25
M <- 10

# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##
#Looking at general patterns

#Calculated bmi from imputed height and weight vs directly imputed bmi - which one to use? At the moment we use the directly imputed bmi.

hist(CSS$bmi,breaks=40)
hist(log(CSS$bmi),breaks=40)
hist(log(log(CSS$bmi)),breaks=40)
hist(sqrt(CSS$bmi),breaks=40)
hist(CSS$age,breaks=40)
hist(CSS$height,breaks=40)
hist(CSS$weight,breaks=40)


bmi_ids <- CSS$CS_ID[which(is.na(CSS$bmi))]
cor(CSS$bmi[CSS$CS_ID%in% bmi_ids][115:length(CSS$bmi[CSS$CS_ID%in% bmi_ids])],(CSS$weight/((CSS$height/100)^2))[CSS$CS_ID%in% bmi_ids][115:length(CSS$bmi[CSS$CS_ID%in% bmi_ids])])
plot(CSS$bmi[CSS$CS_ID%in% bmi_ids][115:length(CSS$bmi[CSS$CS_ID%in% bmi_ids])],(CSS$weight/((CSS$height/100)^2))[CSS$CS_ID%in% bmi_ids][115:length(CSS$bmi[CSS$CS_ID%in% bmi_ids])],ylim=c(0,50))


hist(base_data$bmi,breaks=40,xlim=c(0,50))
hist(log(base_data$bmi),breaks=40)
hist(log(log(base_data$bmi)),breaks=40)
hist(sqrt(base_data$bmi),breaks=40)
hist(base_data$age,breaks=40)
hist(base_data$height,breaks=40)
hist(base_data$weight,breaks=40)


# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#Regression analyses


# --------------------------------------------------------------------------- ##
## Cross-sectional analysis of night-time smarpthone use and BMI continous
# --------------------------------------------------------------------------- ##

#### BASE POPULATION

## test kontinuert bmi
#plot(fitted(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))),residuals(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))))
#hist(residuals(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))),xlim=c(-20,20),breaks=200)
#hist(simulate(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1)))$sim_1,breaks=40) #The bell-shape is not that well suited

#Conclusion: The lm is not by itself appropriate for describing the distribution of BMI.
#Alternative: Pretty good fit. A general family of models.


#Confidence intervals and estimates for continuous BMI:

## MobileUseNight and bmi in baseline sample:

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ mobileUseNight+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseNight","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG) #May use BCS instead of BCCG which corresponds to using a t distribution instead of normal. This can fit heavier tails, though in this case a very large df is fitted, meaning that there is not much difference.
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_baseNight <- miceadds::pool_mi(qhat = coefs, u = vcovs)
#pool_inf_baseNight$qbar
#pool_inf_baseNight$ubar
#pool_inf_baseNight$ba
#pool_inf_baseNight$pval


#Seems that we can get stable contrats of the mean (taking in varying medians), in spite of skewness.

#One slightly hacky way to achieve this may be to take one of the fitted models created by fit() and replace the stored coefficients with the final pooled estimates. I haven't done detailed testing but it seems to be working on this simple example:
m$mu.coefficients <- pool_inf_baseNight$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]


#Confidence intervals:
lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[2,5],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[2,1],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[2,6],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[3,5],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[3,1],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[3,6],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[4,5],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[4,1],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[4,6],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

## estimates og 95%CI for mobileUseNight og continous BMI
confints_base_Night <- cbind(c(lowerCat2,lowerCat3,lowerCat4),c(estCat2,estCat3,estCat4),c(upperCat2,upperCat3,upperCat4))

# --------------------------------------------------------------------------- ##
# MobileUseBeforeSleep and BMI continous in the baseline Citizen Science Sample

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ mobileUseBeforeSleep+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseBeforeSleep","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG) #May use BCS instead of BCCG which corresponds to using a t distribution instead of normal. This can fit heavier tails, though in this case a very large df is fitted, meaning that there is not much difference.
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_baseBefore <- miceadds::pool_mi(qhat = coefs, u = vcovs)
#pool_inf_baseBefore$qbar
#pool_inf_baseBefore$ubar
#pool_inf_baseBefore$ba
#pool_inf_baseBefore$pval


#Seems that we can get stable contrats of the mean (taking in varying medians), in spite of skewness.

#One slightly hacky way to achieve this may be to take one of the fitted models created by fit() and replace the stored coefficients with the final pooled estimates. I haven't done detailed testing but it seems to be working on this simple example:
m$mu.coefficients <- pool_inf_baseBefore$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]

#Confidence intervals:
lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[2,5],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[2,1],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[2,6],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[3,5],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value - integrate(function(y) y*dBCCG(x=y,mu=1,sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[3,1],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[3,6],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[4,5],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[4,1],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[4,6],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[5,5]+1,sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[5,1],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseBefore)[5,6],sigma=exp(pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

## 95%CI for mobileUseBefore sleep og continuos BMI
confints_base_Before <- cbind(c(lowerCat2,lowerCat3,lowerCat4,lowerCat5),c(estCat2,estCat3,estCat4,estCat5),c(upperCat2,upperCat3,upperCat4,upperCat5))

# --------------------------------------------------------------------------- ##
## test for trend (night-time smartphone use and BMI continuous) in baseline Citizen Science sample (table 2 in paper)
# --------------------------------------------------------------------------- ##
#Because the skewness is constant in covariates, the shift between median and mean is constant in covariates.
#We can thus simply make a test for trend on the median to get a p value, as median equality <=> mean equality.

# Smartphone use during sleep period and BMI continuous:

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ as.numeric(mobileUseNight)+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseNight","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_baseTrendNight <- miceadds::pool_mi(qhat = coefs, u = vcovs)

#Going to means:

m$mu.coefficients <- pool_inf_baseTrendNight$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]

#interval
lowerCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendNight)[2,5],sigma=exp(pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendNight)[2,1],sigma=exp(pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendNight)[2,6],sigma=exp(pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

summary(pool_inf_baseTrendNight)$p[2]

# --------------------------------------------------------------------------- ##
## test for trend: Smartphone use Before sleep and BMI continous (table 2 in paper)

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseBeforeSleep","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_baseTrendBefore <- miceadds::pool_mi(qhat = coefs, u = vcovs)

#Going to means:

m$mu.coefficients <- pool_inf_baseTrendBefore$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]

#interval
lowerCatTrendBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendBefore)[2,5],sigma=exp(pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCatTrendBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendBefore)[2,1],sigma=exp(pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCatTrendBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendBefore)[2,6],sigma=exp(pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendBefore$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

summary(pool_inf_baseTrendBefore)$p[2]

# Collected trend intervals
confints_baseTrend <-rbind(c(lowerCatTrendNight,estCatTrendNight,upperCatTrendNight,summary(pool_inf_baseTrendNight)$p[2]),
                           c(lowerCatTrendBS,estCatTrendBS,upperCatTrendBS,summary(pool_inf_baseTrendBefore)$p[2]))


# --------------------------------------------------------------------------- ##
## cross-sectional associations between night-time smartphone use and bmi (25, 30) in baseline Citizen Science sample (table 2 in paper)
# --------------------------------------------------------------------------- ##

## Using the mice package with mids objects

# smartphone use during the sleep period:
mod25 <- with(base_data_mids,glm(bmi25~(mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod30 <- with(base_data_mids,glm(bmi30~(mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))

## test for trend
## BMI > 25
TEST2 <- with(base_data_mids,glm((bmi>=25)~((as.numeric(mobileUseNight))+age+gender+education+occupation), weights=sample_weights,family=binomial))
test2Night <- summary(pool(TEST2), conf.int = T)

## BMI > 30
TEST <- with(base_data_mids,glm((bmi>=30)~((as.numeric(mobileUseNight))+age+gender+education+occupation), weights=sample_weights,family=binomial))
testTNight <- summary(pool(TEST), conf.int = T)

## OR for BMI>30
model30Night <- summary(pool(mod30),conf.int = T)
cbind(exp(model30Night$estimate),
exp(model30Night$`2.5 %`),
exp(model30Night$`97.5 %`))

## OR for BMI >25
model25Night <- summary(pool(mod25), conf.int=T)
cbind(exp(model25Night$estimate),
exp(model25Night$`2.5 %`),
exp(model25Night$`97.5 %`))

# --------------------------------------------------------------------------- ##
# Smartphone use Before sleep and BMI (25 or 30)
mod25 <- with(base_data_mids,glm(bmi25~(mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod30 <- with(base_data_mids,glm(bmi30~(mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))

## test for trend
## BMI > 25
TEST2 <- with(base_data_mids,glm((bmi>=25)~(as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))
test2Before <- summary(pool(TEST2), conf.int = T)

## BMI >30
TEST <- with(base_data_mids,glm((bmi>=30)~(as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))
testTBefore <- summary(pool(TEST), conf.int = T)

## OR for BMI>30
model30Before <- summary(pool(mod30),conf.int = T)
cbind(exp(model30Before$estimate),
      exp(model30Before$`2.5 %`),
      exp(model30Before$`97.5 %`))

## OR for BMI >25
model25Before <- summary(pool(mod25), conf.int=T)
cbind(exp(model25Before$estimate),
      exp(model25Before$`2.5 %`),
      exp(model25Before$`97.5 %`))


# --------------------------------------------------------------------------- ##
## longitudinal analysis of risk scores of smartphone behavior and changes in BMI
# --------------------------------------------------------------------------- ##
#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

# night-time smartphone use and changes in BMI (25 or 30):

# --------------------------------------------------------------------------- ##
## from below 25 to above 25

## smartphone use before sleep onset:
model25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseBeforeSleep.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
model_summary25 <- summary(pool(with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseBeforeSleep.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))), conf.int = T)

exp(cbind(model_summary25$estimate[2:5],model_summary25$`2.5 %`[2:5],model_summary25$`97.5 %`[2:5]))

## smartphone use during the seep period
model25Night <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
model_summary25Night <- summary(pool(with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))), conf.int = T)

exp(cbind(model_summary25Night$estimate[2:4],model_summary25Night$`2.5 %`[2:4],model_summary25Night$`97.5 %`[2:4]))

## test for trend (smartphone use before sleep)
test25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (as.numeric(mobileUseBeforeSleep.y):followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
testT25 <- summary(pool(test25), conf.int = T)
#anova(test25,model25)

## smartphone use during the sleep period
#model25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseNight.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
#model_summary25<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseNight.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))), conf.int = T)
#exp(cbind(model_summary25$estimate[1:4],model_summary25$`2.5 %`[1:4],model_summary25$`97.5 %`[1:4]))

## test for trend (smartphone use during the sleep period)

#test25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (as.numeric(mobileUseNight.y))+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial)
#testT25 <- summary(pool(test25), conf.int = T)

test25Night <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (as.numeric(mobileUseNight.y):followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
testT25Night <- summary(pool(test25Night), conf.int = T)

#anova(test25,model25)
# --------------------------------------------------------------------------- ##

## from below 30 to above 30

## smartphone use before sleep
model30 <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseBeforeSleep.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
model_summary30 <- summary(pool(with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseBeforeSleep.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))), conf.int = T)
exp(cbind(model_summary30$estimate[2:5],model_summary30$`2.5 %`[2:5],model_summary30$`97.5 %`[2:5]))

## test for trend
test30 <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (as.numeric(mobileUseBeforeSleep.y):followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
testT30 <- summary(pool(test30), conf.int=T)
#anova(test30,model30)

## smartphone use during the sleep period
model30Night <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
model_summary30Night <- summary(pool(with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))), conf.int = T)
exp(cbind(model_summary30Night$estimate[2:4],model_summary30Night$`2.5 %`[2:4],model_summary30Night$`97.5 %`[2:4]))

## test for trend
test30Night <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (as.numeric(mobileUseNight.y):followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
testT30Night <- summary(pool(test30Night), conf.int=T)

# --------------------------------------------------------------------------- ##
## BMI continous

## Modelling numeric difference in bmi between baseline and followup
hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))

plot(fitted(lm(difference~(mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),
     residuals(lm(difference~(mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))))
hist(residuals(lm(difference~(mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),breaks=50)

## smartphone use during sleep period and changes in BMI:
#m <- lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*followup_time-selfScoreCat.y-age.y-gender.y-education.y-occupation.y,weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","selfScoreCat.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights")]))
m <- lm(bmi.fu~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","mobileUseNight.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_diff_Night <- summary(pool(with(bmi_followup_mids,lm(bmi.fu~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))), conf.int = T)

cbind(model_summary_diff_Night$estimate[18:20],model_summary_diff_Night$`2.5 %`[18:20],model_summary_diff_Night$`97.5 %`[18:20])

## test for trend (mobileUseNight)
test_numNight <- with(bmi_followup_mids,lm(bmi.fu~((as.numeric(mobileUseNight.y)):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))
test_TnumNight <- summary(pool(test_numNight), conf.int=T)


## smartphone use Before Sleep and changes in BMI (OBS: gentagelse fra ovenstående??!?)

## from below 25 to above 25
#model25Before <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
#model_summary25Before<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))), conf.int = T)
#exp(cbind(model_summary25$estimate[1:4],model_summary25$`2.5 %`[1:4],model_summary25$`97.5 %`[1:4]))


## test for trend
#test25Before <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (as.numeric(mobileUseBeforeSleep.y)+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
#testT25Before <- summary(pool(test25Before), conf.int = T)
#anova(test25,model25)

## from below 30 to above 30
#model30Before <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
#model_summary30Before<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))), conf.int = T)
#exp(cbind(model_summary30$estimate[1:4],model_summary30$`2.5 %`[1:4],model_summary30$`97.5 %`[1:4]))
# test for trend
#test30Before <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (as.numeric(mobileUseBeforeSleep.y)+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
#testT30Before <- summary(pool(test30Before), conf.int = T)
#anova(test30,model30)

## smartphone use before sleep onset and bmi continous
## Modelling numeric difference in bmi between baseline and followup
hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))

plot(fitted(lm(difference~(mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),
     residuals(lm(difference~(mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))))
hist(residuals(lm(difference~(mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),breaks=50)


#m <- lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*followup_time-selfScoreCat.y-age.y-gender.y-education.y-occupation.y,weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","selfScoreCat.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights")]))

m <- lm(bmi.fu~((mobileUseBeforeSleep.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","mobileUseBeforeSleep.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_diff_Before <- summary(pool(with(bmi_followup_mids,lm(bmi.fu~((mobileUseBeforeSleep.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))), conf.int = T)

cbind(model_summary_diff_Before$estimate[18:21],model_summary_diff_Before$`2.5 %`[18:21],model_summary_diff_Before$`97.5 %`[18:21])

## test for trend (mobileUseBefore)
test_numBefore <- with(bmi_followup_mids,lm(bmi.fu~((as.numeric(mobileUseBeforeSleep.y)):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))
test_TnumBefore <- summary(pool(test_numBefore), conf.int=T)


#Generally:
#Wald intervals with Robust=T are better for misspecified models.
#Profile likelihood intervals are better for models that are close to correct (mostly so for smaller sample sizes).


# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#Tracking data: Population sample (random sample) - analyses of numerical bmi and of indicators

hist(pop_track$bmi,breaks=50,xlim=c(0,50))

ggplot(pop_track, aes(x = factor(cluster))) +
  geom_bar()

#analyses

## regression analysis of clusters of night-time smartphone use and overweight/obesity

## Continuous Outcome

## Smartphone use during the sleep period and BMI continous in population sample:

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){ #Slow
  m <- gamlss(bmi~(mobileUseNight+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("bmi","mobileUseNight","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTNight <- miceadds::pool_mi(qhat = coefs, u = vcovs)


lowerNight2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerNight3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[3,5]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[3,1]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[3,6]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerNight4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[4,5]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[4,1]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNight)[4,6]+10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoTNight <- cbind(c(lowerNight2,lowerNight3,lowerNight4),
                                   c(estNight2,estNight3,estNight4),
                                   c(upperNight2,upperNight3,upperNight4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNight$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


#Trend:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(as.numeric(mobileUseNight)+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("mobileUseNight","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTNightTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNightTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatNight <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNightTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNightTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackNoTTrendNight <- rbind(c(lowerCatNight,estCatNight,upperCatNight) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value)

## BeforeSleep

#Without adjustment for tracking:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){ #Slow
  m <- gamlss(bmi~(mobileUseBeforeSleep+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("bmi","mobileUseBeforeSleep","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTBefore <- miceadds::pool_mi(qhat = coefs, u = vcovs)


lowerBS2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[3,5]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[3,1]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[3,6]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[4,5]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[4,1]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[4,6]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[5,5]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[5,1]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBefore)[5,6]+10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackNoTBefore <- cbind(c(lowerBS2,lowerBS3,lowerBS4,lowerBS5),
                                    c(estBS2,estBS3,estBS4,estBS5),
                                    c(upperBS2,upperBS3,upperBS4,upperBS5))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBefore$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

#Trend:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(as.numeric(mobileUseBeforeSleep)+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("mobileUseBeforeSleep","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTBeforeTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerCatBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBeforeTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatBS <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBeforeTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTBeforeTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoTTrendBefore <- rbind(c(lowerCatBS,estCatBS,upperCatBS) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTBeforeTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value)


## Maximal posterior probability assignment 6 clusters

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}
pool_inf_PopTrackNoS_mp <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[5,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[5,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[5,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[6,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[6,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[6,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoS_mpSix <- cbind(c(lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6),
                              c(estClust2,estClust3,estClust4,estClust5,estClust6),
                              c(upperClust2,upperClust3,upperClust4,upperClust5,upperClust6))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
## prediction

m <- gamlss(bmi~(cluster+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmipop_sixmax <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("cluster","bmi","age","sex","education","occupation","sample_weights","imputation")])

MSEbmipopsixmax <- mean((predbmipop_sixmax - pop_track$bmi[pop_track$imputation!=0])^2)


## Maximal posterior probability assignment 4 clusters

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster.y+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster.y","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS_mp <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoS_mpFour <- cbind(c(lowerClust2,lowerClust3,lowerClust4),
                              c(estClust2,estClust3,estClust4),
                              c(upperClust2,upperClust3,upperClust4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
## prediction

m <- gamlss(bmi~(cluster.y+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster.y","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmipop_fourmax <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("cluster.y","bmi","age","sex","education","occupation","sample_weights","imputation")])

MSEbmipopfourmax <- mean((predbmipop_fourmax - pop_track$bmi[pop_track$imputation!=0])^2)


# --------------------------------------------------------------------------- ##
### Binary Outcomes for population sample
# --------------------------------------------------------------------------- ##

## BMI > 25

## Maximal posterior probability assignment

# Six clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No_mpSix <- summary(pool(Random25No), conf.int = T)
m <- glm((bmi>=25) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random25No)$pooled$estimate
predpopbin25_maxsix <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin25_predmaxsix <- mean((expit(predpopbin25_maxsix)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)

# Four clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No_mpFour <- summary(pool(Random25No), conf.int = T)
m <- glm((bmi>=25) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random25No)$pooled$estimate
predpopbin25_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin25_predmaxfour <- mean((expit(predpopbin25_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)



## SelfScoreCat and BMI>25 (no adjustment for tracking)
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25NoTNight <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25NoTNight <- summary(pool(Random25NoTNight), conf.int=T)
cbind(exp(modelRandom25NoTNight$estimate),
exp(modelRandom25NoTNight$`2.5 %`),
exp(modelRandom25NoTNight$`97.5 %`))

summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseBeforeSleep+age+sex+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25NoTBefore <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseBeforeSleep+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25NoTBefore <- summary(pool(Random25NoTBefore), conf.int=T)
cbind(exp(modelRandom25NoTBefore$estimate),
      exp(modelRandom25NoTBEfore$`2.5 %`),
      exp(modelRandom25NoTBefore$`97.5 %`))

#test for trend (selfscoreCat uden cluster)
Random25NoTestNight <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseNight)+age+sex+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25NoTestNight), conf.int=T)

Random25NoTestBefore <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseBeforeSleep)+age+sex+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25NoTestBefore), conf.int=T)

# --------------------------------------------------------------------------- ##
## BMI > 30

## BMI >30 

## Maximal posterior probability assignment

# Six clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No_mpSix <- summary(pool(Random30No), conf.int = T)
m <- glm((bmi>=30) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random30No)$pooled$estimate
predpopbin30_maxsix <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin30_predmaxsix <- mean((expit(predpopbin30_maxsix)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)

# Four clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No_mpFour <- summary(pool(Random30No), conf.int = T)
m <- glm((bmi>=30) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random30No)$pooled$estimate
predpopbin30_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin30_predmaxfour <- mean((expit(predpopbin30_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)


## no adjustment for tracking
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30NoTNight <- with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30NoTNight <- summary(pool(Random30NoTNight), conf.int=T)
cbind(exp(modelRandom30NoTNight$estimate),
      exp(modelRandom30NoTNight$`2.5 %`),
      exp(modelRandom30NoTNight$`97.5 %`))

summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (mobileUseBeforeSleep+age+sex+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30NoTBefore <- with(pop_track_mids,glm((bmi>=30) ~ (mobileUseBeforeSleep+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30NoTBefore <- summary(pool(Random30NoTBefore), conf.int=T)
cbind(exp(modelRandom30NoTBefore$estimate),
      exp(modelRandom30NoTBEfore$`2.5 %`),
      exp(modelRandom30NoTBefore$`97.5 %`))

#test for trend (selfscoreCat uden cluster)
Random30NoTestNight <- with(pop_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseNight)+age+sex+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random30NoTestNight), conf.int=T)

Random30NoTestBefore <- with(pop_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseBeforeSleep)+age+sex+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random30NoTestBefore), conf.int=T)


###############################################################################
###############################################################################
###############################################################################

#Analysis of the clinical sample data - interest in biomarkers

# --------------------------------------------------------------------------- ##
## BMI
table(clinical_sample$bmi.clinical)
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
summary(clinical_sample$bmi.clinical, useNA="always")

## BMI clinical
table(clinical_sample$bmi.clinical)
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable( ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))



# --------------------------------------------------------------------------- ##
## descriptive of clinical sample
## age
publish(univariateTable(mobileUseNight ~ age.x,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ age.x,data=clinical_sample, column.percent=TRUE))
clinical_sample$age.x <- as.numeric(clinical_sample$age.x)

## BMI
publish(univariateTable(mobileUseNight ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

#The subjects are scoring in the high end. Is this an issue or a characteristic of the data?

#Introducing interesting derived variables

clinical_sample$bmi <- as.numeric(clinical_sample$bmi.clinical)
clinical_sample$bmi25 <- as.numeric(clinical_sample$bmi.clinical>=25)
clinical_sample$bmi30 <- as.numeric(clinical_sample$bmi.clinical>=30)

## age at clinical examination (OBS. NOGET ER GALT I DENNE VARIABLE)
table(clinical_sample$age.y, useNA="always")
clinical_sample$age<- as.numeric(str_c(substr(clinical_sample$age.y,1,1),substr(clinical_sample$age.y,2+(mod(nchar(clinical_sample$age.y),4)==1),2+(mod(nchar(clinical_sample$age.y),4)==1)),".",
                 substr(clinical_sample$age.y,3+(mod(nchar(clinical_sample$age.y),4)!=3),3+(mod(nchar(clinical_sample$age.y),4)!=3))))


# --------------------------------------------------------------------------- ##

## systolic blood pressure
clinical_sample$sbp<-rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ sbp,data=clinical_sample, column.percent=TRUE))

# diastolic blood pressure
clinical_sample$dbp<-rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ dbp,data=clinical_sample, column.percent=TRUE))

## hip waist ratio

clinical_sample$ratiowaisthip <- as.numeric(clinical_sample$ratiowaisthip)
publish(univariateTable(selfScoreCat ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))


## bmi clinical
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable(selfScoreCat ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))


#hdl, ldl, vldl, t_cholesterol, triglycerid, hba1c, (glucose), waist, hip, ratio waist hip, systolic bp og distolic bp 1-3: Ift. selvrapporteringer og tracking clusters


#Looks into data

hist(as.numeric(clinical_sample$hdl),breaks=20)
hist(as.numeric(clinical_sample$ldl),breaks=20)
hist(as.numeric(clinical_sample$t_cholesterol),breaks=20)
hist(as.numeric(clinical_sample$triglycerids),breaks=40)
hist(as.numeric(clinical_sample$hba1c),breaks=40)
hist(as.numeric(clinical_sample$glucose),breaks=40)
hist(as.numeric(clinical_sample$waist),breaks=20)
hist(as.numeric(clinical_sample$hip),breaks=20)
hist(as.numeric(clinical_sample$ratiowaisthip),breaks=20)
hist(as.numeric(clinical_sample$bmi.clinical),breaks=20)
hist(rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T),breaks=20)
hist(rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T),breaks=20)

## HDL
clinical_sample$hdl <- as.numeric(clinical_sample$hdl)
publish(univariateTable(selfScoreCat ~ hdl,data=clinical_sample, column.percent=TRUE))

## LDL
clinical_sample$ldl <- as.numeric(clinical_sample$ldl)
publish(univariateTable( selfScoreCat~ ldl,data=clinical_sample, column.percent=TRUE))

## VLDL
clinical_sample$vldl <- as.numeric(clinical_sample$vldl)
publish(univariateTable( selfScoreCat~ vldl,data=clinical_sample, column.percent=TRUE))

## total cholesterol
clinical_sample$t_cholesterol <- as.numeric(clinical_sample$t_cholesterol)
publish(univariateTable(selfScoreCat ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))

## triglycerides
clinical_sample$triglycerids <- as.numeric(clinical_sample$triglycerids)
publish(univariateTable(selfScoreCat ~ triglycerids,data=clinical_sample, column.percent=TRUE))

## hba1c
clinical_sample$hba1c <- as.numeric(clinical_sample$hba1c)
publish(univariateTable(selfScoreCat ~ hba1c,data=clinical_sample, column.percent=TRUE))


# --------------------------------------------------------------------------- ##

#Models - multiple testing issue if we are going to 'pick and choose' which responses we would like to look at.
table(clinical_sample$age.x)

#Transforming to mids for modelling and inference

clinical_sample$cluster.y <- factor(clinical_sample$cluster.y,levels = c("Cluster 3", "Cluster 2", "Cluster 4", "Cluster 1"))

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")


## Maximal posterior probability assignment: Six clusters

dbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(dbp) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

vldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

df_ints_mpSix <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_mpSix) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")


## Maximal posterior probability assignment: Four clusters

dbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(dbp) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

vldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

## estimater for four clusters i clinical sample
df_ints_mpFour <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_mpFour) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")



## night-time smartphone use 

predict.person <- data.frame("age"=18,education="",occupation="Employed")

## mobileUseBeforeSleep and biomarkers

dbp_intS <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intS <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

## Self assesment: Night

predict.person <- data.frame("age"=18,education="",occupation="Employed")

dbp_intsNight <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseNight+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intsNight <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseNight+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])


df_intsNight <- data.frame(rbind(hdl_intsNight[2,],ldl_intsNight[2,],vldl_intsNight[2,],t_chol_intsNight[2,],tri_intsNight[2,],hba1c_intsNight[2,],dbp_intsNight[2,],sbp_intsNight[2,],wh_intsNight[2,],glu_intsNight[2,],bmi_intsNight[2,]),
                      rbind(hdl_intsNight[3,],ldl_intsNight[3,],vldl_intsNight[3,],t_chol_intsNight[3,],tri_intsNight[3,],hba1c_intsNight[3,],dbp_intsNight[3,],sbp_intsNight[3,],wh_intsNight[3,],glu_intsNight[3,],bmi_intsNight[3,]),
                      rbind(hdl_intsNight[4,],ldl_intsNight[4,],vldl_intsNight[4,],t_chol_intsNight[4,],tri_intsNight[4,],hba1c_intsNight[4,],dbp_intsNight[4,],sbp_intsNight[4,],wh_intsNight[4,],glu_intsNight[4,],bmi_intsNight[4,]))

colnames(df_intsNight) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsNight) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsNight

## mobileUseNight and biomarkers
dbp_intS <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseNight+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intS <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseNight+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_intsS <- data.frame(rbind(hdl_intS[2,],ldl_intS[2,],vldl_intS[2,],t_chol_intS[2,],tri_intS[2,],hba1c_intS[2,],dbp_intS[2,],sbp_intS[2,],wh_intS[2,],glu_intS[2,],bmi_intS[2,]),
                       rbind(hdl_intS[3,],ldl_intS[3,],vldl_intS[3,],t_chol_intS[3,],tri_intS[3,],hba1c_intS[3,],dbp_intS[3,],sbp_intS[3,],wh_intS[3,],glu_intS[3,],bmi_intS[3,]),
                       rbind(hdl_intS[4,],ldl_intS[4,],vldl_intS[4,],t_chol_intS[4,],tri_intS[4,],hba1c_intS[4,],dbp_intS[4,],sbp_intS[4,],wh_intS[4,],glu_intS[4,],bmi_intS[4,]))

colnames(df_intsS) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsS) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsS


## Smartphone use Before sleep and biomarkers

dbp_intsBefore <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intsBefore <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_intsBefore <- data.frame(rbind(hdl_intsBefore[2,],ldl_intsBefore[2,],vldl_intsBefore[2,],t_chol_intsBefore[2,],tri_intsBefore[2,],hba1c_intsBefore[2,],dbp_intsBefore[2,],sbp_intsBefore[2,],wh_intsBefore[2,],glu_intsBefore[2,],bmi_intsBefore[2,]),
                                rbind(hdl_intsBefore[3,],ldl_intsBefore[3,],vldl_intsBefore[3,],t_chol_intsBefore[3,],tri_intsBefore[3,],hba1c_intsBefore[3,],dbp_intsBefore[3,],sbp_intsBefore[3,],wh_intsBefore[3,],glu_intsBefore[3,],bmi_intsBefore[3,]),
                                rbind(hdl_intsBefore[4,],ldl_intsBefore[4,],vldl_intsBefore[4,],t_chol_intsBefore[4,],tri_intsBefore[4,],hba1c_intsBefore[4,],dbp_intsBefore[4,],sbp_intsBefore[4,],wh_intsBefore[4,],glu_intsBefore[4,],bmi_intsBefore[4,]))

colnames(df_intsBefore) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsBefore) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsBefore

