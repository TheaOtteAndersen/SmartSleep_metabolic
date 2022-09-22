
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

N_imp = 25

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
#subject_tracking_six_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_clusters.csv") ## forkert navn p? csv-fil?
## Er det den korrekt fil??
subject_tracking_six_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_clusters_four_and_six.csv")

subject_tracking_four_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_four_clusters.csv")


## Collecting the two clusterings in one file

subject_tracking_clusters <- left_join(subject_tracking_six_clusters,subject_tracking_four_clusters[,c("userid","cluster","cluster1prob","cluster2prob","cluster3prob","cluster4prob","description","state0prob","state1prob","state2prob","state3prob")],by="userid")
subject_tracking_clusters <- rename(subject_tracking_clusters,cluster1prob=cluster1prob.x,cluster2prob=cluster2prob.x,cluster3prob=cluster3prob.x,cluster4prob=cluster4prob.x,
                                    state0prob=state0prob.x,state1prob=state1prob.x,state2prob=state2prob.x,state3prob=state3prob.x,cluster=cluster.x,description=description.x)

## load baseline data
base_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv")
base_data$mobileUseNight <- factor(base_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))


## load followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")
CSS$mobileUseNight <- factor(CSS$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))


## load population sample
pop_data <-read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv")
pop_data$mobileUseNight <- factor(pop_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))


## load clinical data (survey and clinical data)
#load("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample/full_imp_clinical.RData")
#clin_data <- full_imp_clinical
clin_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample/imp_clinical.csv")
clin_clinical <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/SmartSleep Clinical/Data/R?data/SmartSleepClinicalData.csv")
## Er dette det rigtige clinical data? Det er ikke alle fra clin_data der er i clin_clinical og omvendt?
unique(clin_data$PNR[!clin_data$PNR %in% clin_clinical$cpr])
unique(clin_clinical$cpr[!clin_clinical$cpr %in% clin_data$PNR])

clin_data$mobileUseNight <- factor(clin_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))

# --------------------------------------------------------------------------- ##

#Baseline data with self-reports

## if no mobile phone = NA
base_data$pmpuScale[base_data$mobilephone=="No mobile phone"] <- NA

## night-time smartphone use
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))

## weight
publish(univariateTable( ~ weight,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ bmi,data=base_data, column.percent=TRUE))
table(base_data$height[base_data$imputation==0])

## bmi
base_data$bmi30 <- (base_data$bmi>=30)
base_data$bmi25 <- (base_data$bmi>=25)

base_data_mids <- as.mids(base_data,.imp="imputation")


# --------------------------------------------------------------------------- ##
#Followup sample - using quartile levels from baseline sample

table(CSS$mobileUseNight, useNA="always")

## merge survey and tracking data 
CSS_track <- left_join(CSS,subject_tracking_clusters,by="userid")
CSS_track$sample_weights <- as.numeric(CSS_track$sample_weights)

CSS_track_mids<-as.mids(CSS_track,.imp="imputation",.id="userid")

# --------------------------------------------------------------------------- ##
#Merging base and followup

#BMI followup difference - match with emailAddress or CS_ID
table(CSS$weight)
table(base_data$weight)

#y: base, x: followup

bmi_followup <- rename(inner_join(CSS,base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

#weight_followup <- rename(inner_join(CSS,base_data,by=c("CS_ID","imputation")),weight.base=weight.y,weight.fu=weight.x)

## weight difference between follow-up and baselnie
bmi_followup$differenceWeight <- bmi_followup$weight.y-bmi_followup$weight.x
mean(bmi_followup$differenceWeight[!is.na(bmi_followup$differenceWeight)])


## store forskelle i weight from baseline to follow-up
weight <- subset(bmi_followup, select=c(weight.x, weight.y, differenceWeight))
weight <- weight[bmi_followup$imputation==0,]

## difference in BMI mellem follow-up og baseline
bmi_followup$difference <- bmi_followup$bmi.fu-bmi_followup$bmi.base
mean(bmi_followup$difference[!is.na(bmi_followup$difference)])
table(bmi_followup$difference)


BMI <- subset(bmi_followup, select=c(bmi.fu, bmi.base, difference))

## naming sample_weights and constructing followup_time
bmi_followup$sample_weights <- bmi_followup$sample_weights.y
bmi_followup$followup_time <- (as.Date(str_c(substr(bmi_followup$EndDate,7,10),"-",substr(bmi_followup$EndDate,4,5),"-",substr(bmi_followup$EndDate,1,2)))-as.Date(str_c(substr(bmi_followup$responseDate.y,7,10),"-",substr(bmi_followup$responseDate.y,4,5),"-",substr(bmi_followup$responseDate.y,1,2))))/365.25

## subset according to non-missing weight in the unimputed data
bmi_followup_complete <- subset(bmi_followup, imputation == 0 & !(is.na(height.x) | is.na(height.y)))
bmi_followup_complete <- subset(bmi_followup_complete, imputation == 0 & !(is.na(weight.x) | is.na(weight.y)))

bmi_followup_complete <- subset(bmi_followup, userid %in% bmi_followup_complete$userid) # 1768 are left from the total of 1885.


## difference in BMI and weight
table(bmi_followup_work$differenceWeight)
table(bmi_followup_work$difference)

## difference in bmi and weight according to sex
bmi_followupMen <- subset(bmi_followup_work, sex.x=="Man")
bmi_followupWomen <- subset(bmi_followup_work, sex.x=="Woman")

## mean difference in men and women (BMI and weight)
## BMI
mean(bmi_followupMen$difference[!is.na(bmi_followupMen$difference)])
mean(bmi_followupWomen$difference[!is.na(bmi_followupWomen$difference)])

## Weight
mean(bmi_followupMen$differenceWeight[!is.na(bmi_followupMen$differenceWeight)])
mean(bmi_followupWomen$differenceWeight[!is.na(bmi_followupWomen$differenceWeight)])

BMI <- subset(bmi_followup_work, select=c(bmi.fu, bmi.base, difference))
weight <- subset(bmi_followup_work, select=c(weight.x, weight.y, differenceWeight))


## Assigning mids objects
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")
bmi_followup_complete_mids <- as.mids(bmi_followup_complete,.imp="imputation")

# --------------------------------------------------------------------------- ##
#Population sample

##
table(pop_data$mobileUseNight, useNA="always")

## merge tracking and survey data for population sample
pop_track <- left_join(pop_data,subject_tracking_clusters,by="userid")
pop_track$sample_weights<-as.numeric(pop_track$sample_weights)

## omkategoriser 4 clusters
table(pop_track$description.y, pop_track$cluster.y)

prop.table(table(pop_track$description.y, pop_track$cluster.y))

pop_track$cluster.y <- factor(pop_track$cluster.y,levels = c("Cluster 3", "Cluster 2", "Cluster 4", "Cluster 1"))

pop_track_mids<-as.mids(pop_track,.imp="imputation",.id="userid")

## Clinical Data

## merge survey and clinical data
clinical_sample <- rename(left_join(clin_data,rename(clin_clinical,PNR=cpr),by="PNR"),bmi.self=bmi.x , bmi.clinical=bmi.y)
## merge with tracking data
clinical_sample <- left_join(clinical_sample,subject_tracking_clusters,by="userid")

## Relevelling rare categories
clinical_sample$educationW <- clinical_sample$education
clinical_sample$educationW[clinical_sample$education %in% c("Technical vocational education", "short cycle higher education")] <- "Technical vocational education or short cycle higher education"
clinical_sample$educationW[clinical_sample$education %in% c("Primary school", "Other")] <- "Primary school or other"

clinical_sample$occupationW <- clinical_sample$occupation
clinical_sample$occupationW[!(clinical_sample$occupation %in% c("employed","student"))] <- "Other"

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

#Confidence intervals for means: Calculated by integration over the model density. By (v) in Ferrari & Fumes we can get mean contrasts by integrating with the median contrasts.
lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[2,5],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[2,1],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[2,6],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[3,5],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[3,1],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[3,6],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[4,5],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[4,1],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseNight)[4,6],sigma=exp(pool_inf_baseNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

## estimates for 95%CI for mobileUseNight og continous BMI
confints_base_Night <- cbind(c(lowerCat2,lowerCat3,lowerCat4),c(estCat2,estCat3,estCat4),c(upperCat2,upperCat3,upperCat4))

# --------------------------------------------------------------------------- ##
## test for trend (night-time smartphone use and BMI continuous) in baseline Citizen Science sample (table 2 in paper)
# --------------------------------------------------------------------------- ##
#Because the skewness is constant in covariates, the shift between median and mean is constant in covariates.
#We can simply make a test for trend on the median to get a p value for trend on the mean, as median equality <=> mean equality (by (v) in Ferrari and Fumes).

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

#interval
lowerCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendNight)[2,5],sigma=exp(pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendNight)[2,1],sigma=exp(pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrendNight)[2,6],sigma=exp(pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrendNight$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

summary(pool_inf_baseTrendNight)$p[2]


# --------------------------------------------------------------------------- ##
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
## longitudinal analysis of risk scores of smartphone behavior and changes in BMI
# --------------------------------------------------------------------------- ##
#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup


## BMI continous

## Modelling numeric difference in bmi between baseline and followup
hist(bmi_followup_complete$difference,xlim=c(-10,10),breaks=100)
hist(bmi_followup_complete$weight.x-bmi_followup_complete$weight.y,xlim=c(-10,10),breaks=100)

plot(fitted(lm(difference~(mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),
     residuals(lm(difference~(mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))))
hist(residuals(lm(difference~(mobileUseNight.y:followup_time+followup_time+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),breaks=50)

## smartphone use during sleep period and changes in BMI:
m <- lm(bmi.fu~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","mobileUseNight.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_diff_Night <- summary(pool(with(bmi_followup_complete_mids,lm(bmi.fu~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))), conf.int = T)

cbind(model_summary_diff_Night$estimate[18:20],model_summary_diff_Night$`2.5 %`[18:20],model_summary_diff_Night$`97.5 %`[18:20])

## test for trend (mobileUseNight)
test_numNight <- with(bmi_followup_complete_mids,lm(bmi.fu~((as.numeric(mobileUseNight.y)):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))
test_TnumNight <- summary(pool(test_numNight), conf.int=T)


#Generally:
#Wald intervals with Robust=T are better for misspecified models.
#Profile likelihood intervals are better for models that are close to correct (mostly so for smaller sample sizes).

## changes in weight

## smartphone use during sleep period and changes in BMI:
m <- lm(weight.x~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+weight.y),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","mobileUseNight.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_weightdiff_Night <- summary(pool(with(bmi_followup_complete_mids,lm(weight.x~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+weight.y),weights=sample_weights))), conf.int = T)

cbind(model_summary_weightdiff_Night$estimate[18:20],model_summary_weightdiff_Night$`2.5 %`[18:20],model_summary_weightdiff_Night$`97.5 %`[18:20])

## test for trend (mobileUseNight)
test_numNightWeight <- with(bmi_followup_complete_mids,lm(weight.x~((as.numeric(mobileUseNight.y)):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+weight.y),weights=sample_weights))
test_TnumNightWeight <- summary(pool(test_numNightWeight), conf.int=T)


# -------------------------------------------------------------------------------------------------------------- ##
## cross-sectional analyses of self-reported and clusters of night-time smartphone use and BMI in population sample
# -------------------------------------------------------------------------------------------------------------- ##

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

#Confidence intervals for means: Calculated by integration over the model density. By (v) in Ferrari & Fumes we can get mean contrasts by integrating with the median contrasts.

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


#test for Trend (smartphone use during the sleep period and BMI continous in population sample:
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

#Confidence intervals for means: Calculated by integration over the model density. By (v) in Ferrari & Fumes we can get mean contrasts by integrating with the median contrasts.

lowerCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNightTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatNight <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNightTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTNightTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoTTrendNight <- rbind(c(lowerCatNight,estCatNight,upperCatNight) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTNightTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value)
confints_PopTrackNoTTrendNight <- cbind(confints_PopTrackNoTTrendNight,summary(pool_inf_PopTrackNoTNightTrend)[2,4])

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

#Confidence intervals for means: Calculated by integration over the model density. By (v) in Ferrari & Fumes we can get mean contrasts by integrating with the median contrasts.

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
## prediction:
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.

m <- gamlss(bmi~(cluster+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
# Erstatter estimerede v?rdier:
m$mu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 
# Laver pr?diktion:
predbmipop_sixmax <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("cluster","bmi","age","sex","education","occupation","sample_weights","imputation")])
# Finder MSE: 
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

#Confidence intervals for means: Calculated by integration over the model density. By (v) in Ferrari & Fumes we can get mean contrasts by integrating with the median contrasts.

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
## prediction:
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.

m <- gamlss(bmi~(cluster.y+age+sex+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster.y","bmi","age","sex","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
# Erstatter estimerede v?rdier:
m$mu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 
# Laver pr?diktion:
predbmipop_fourmax <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("cluster.y","bmi","age","sex","education","occupation","sample_weights","imputation")])
# Finder MSE: 
MSEbmipopfourmax <- mean((predbmipop_fourmax - pop_track$bmi[pop_track$imputation!=0])^2)


# --------------------------------------------------------------------------- ##
### Binary Outcomes for population sample
# --------------------------------------------------------------------------- ##

## BMI > 25

## Maximal posterior probability assignment

# Six clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No_mpSix <- summary(pool(Random25No), conf.int = T)
exp(cbind(modelRandom25No_mpSix$estimate[2:6],modelRandom25No_mpSix$`2.5 %`[2:6],modelRandom25No_mpSix$`97.5 %`[2:6]))
## Prediction:
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.
m <- glm((bmi>=25) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
# Erstatter estimerede v?rdier:
m$coefficients <- pool(Random25No)$pooled$estimate
# Laver pr?diktion:
predpopbin25_maxsix <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
# Finder MSE: 
MSEpopbin25_predmaxsix <- mean((expit(predpopbin25_maxsix)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)

# Four clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No_mpFour <- summary(pool(Random25No), conf.int = T)
exp(cbind(modelRandom25No_mpFour$estimate[2:4],modelRandom25No_mpFour$`2.5 %`[2:4],modelRandom25No_mpFour$`97.5 %`[2:4]))

## Prediction: (hvad bruger vi dette til? (12/09/2022))
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.
m <- glm((bmi>=25) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
# Erstatter estimerede v?rdier:
m$coefficients <- pool(Random25No)$pooled$estimate
# Laver pr?diktion:
predpopbin25_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
# Finder MSE: 
MSEpopbin25_predmaxfour <- mean((expit(predpopbin25_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)


## smartphone use during the sleep period and BMI>25 in population sample
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25NoTNight <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom25NoTNight <- summary(pool(Random25NoTNight), conf.int=T)
cbind(exp(modelRandom25NoTNight$estimate[2:4]),exp(modelRandom25NoTNight$`2.5 %`[2:4]),exp(modelRandom25NoTNight$`97.5 %`[2:4]))

#test for trend:
## smartphone use during the sleep period and BMI > 25
Random25NoTestNight <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseNight)+age+sex+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25NoTestNight), conf.int=T)

# --------------------------------------------------------------------------- ##
## sensitivity analyses (further adjusting for physical activity)

## self-reported night-time smartphone use

table(pop_track$physicalActivityDescription)
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+sex+education+occupation+physicalActivityDescription), weights=sample_weights,family=binomial))),conf.int=T)
Random25NoTNight <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+sex+education+occupation+physicalActivityDescription), weights=sample_weights,family=binomial))
modelRandom25NoTNight <- summary(pool(Random25NoTNight), conf.int=T)
cbind(exp(modelRandom25NoTNight$estimate[2:4]),exp(modelRandom25NoTNight$`2.5 %`[2:4]),exp(modelRandom25NoTNight$`97.5 %`[2:4]))

Random25NoTestNight <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseNight)+age+sex+education+occupation+physicalActivityDescription), weights=sample_weights,family=binomial))
summary(pool(Random25NoTestNight), conf.int=T)

# latent clusters (four clusters)
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster.y+age+sex+education+occupation+physicalActivityDescription), weights=sample_weights,family=binomial))
modelRandom25No_mpFour <- summary(pool(Random25No), conf.int = T)
exp(cbind(modelRandom25No_mpFour$estimate[2:4],modelRandom25No_mpFour$`2.5 %`[2:4],modelRandom25No_mpFour$`97.5 %`[2:4]))
## Prediction:
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.
m <- glm((bmi>=25) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
# Erstatter estimerede v?rdier:
m$coefficients <- pool(Random25No)$pooled$estimate
# Laver pr?diktion:
predpopbin25_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
# Finder MSE: 
MSEpopbin25_predmaxfour <- mean((expit(predpopbin25_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)


# --------------------------------------------------------------------------- ##
## BMI > 30

## Maximal posterior probability assignment

# Six clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No_mpSix <- summary(pool(Random30No), conf.int = T)
exp(cbind(modelRandom30No_mpSix$estimate[2:6],modelRandom30No_mpSix$`2.5 %`[2:6],modelRandom30No_mpSix$`97.5 %`[2:6]))
## Prediction:
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.
m <- glm((bmi>=30) ~ (cluster+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
# Erstatter estimerede v?rdier:
m$coefficients <- pool(Random30No)$pooled$estimate
# Laver pr?diktion:
predpopbin30_maxsix <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
# Finder MSE: 
MSEpopbin30_predmaxsix <- mean((expit(predpopbin30_maxsix)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)

# Four clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No_mpFour <- summary(pool(Random30No), conf.int = T)
exp(cbind(modelRandom30No_mpFour$estimate[2:4],modelRandom30No_mpFour$`2.5 %`[2:4],modelRandom30No_mpFour$`97.5 %`[2:4]))

## Prediction:
## vi fitter en enkelt model m med de samme kovariater og s? erstatter vi de fittede parametre i den enkelte model med dem fra vores poolede fit.
## Vi bruger s? objektet m til at lave pr?diktion.
m <- glm((bmi>=30) ~ (cluster.y+age+sex+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
# Erstatter estimerede v?rdier:
m$coefficients <- pool(Random30No)$pooled$estimate
# Laver pr?diktion:
predpopbin30_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
# Finder MSE: 
MSEpopbin30_predmaxfour <- mean((expit(predpopbin30_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)


## Smartphone use during the sleep period and BMI>30
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30NoTNight <- with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+age+sex+education+occupation), weights=sample_weights,family=binomial))
modelRandom30NoTNight <- summary(pool(Random30NoTNight), conf.int=T)
cbind(exp(modelRandom30NoTNight$estimate),
      exp(modelRandom30NoTNight$`2.5 %`),
      exp(modelRandom30NoTNight$`97.5 %`))[2:4,]

#test for trend 
## smartphone use during the sleep period and BMI >30
Random30NoTestNight <- with(pop_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseNight)+age+sex+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random30NoTestNight), conf.int=T)




###############################################################################
###############################################################################
###############################################################################

#Analysis of the clinical sample data - interest in biomarkers
#the clin_clinical data has also the last of the 245 subjects having only clinical information. (But we do not use this person, as there is no questionnaire or tracking)

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

## imputation
table(clinical_sample$imputation)

## age
publish(univariateTable(mobileUseNight ~ age.y,data=clinical_sample, column.percent=TRUE))
clinical_sample$age.x <- as.numeric(clinical_sample$age.x)

## educationW
publish(univariateTable(mobileUseNight ~ educationW,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseNight ~ occupationW,data=clinical_sample, column.percent=TRUE))



## BMI
publish(univariateTable(mobileUseNight ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

## night-time smartphone use
publish(univariateTable(mobileUseNight ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

#The subjects are scoring in the high end. Is this an issue or a characteristic of the data?

#Introducing interesting derived variables

clinical_sample$bmi <- as.numeric(clinical_sample$bmi.clinical)
clinical_sample$bmi25 <- as.numeric(clinical_sample$bmi.clinical>=25)
clinical_sample$bmi30 <- as.numeric(clinical_sample$bmi.clinical>=30)

## age at clinical examination - DOESN'T RUN!!!!!! BECAUSE OF NA's - SEE READING IN DATA COMMENT. USE age.x in NA cases?
table(clinical_sample$age.y, useNA="always")
clinical_sample$age <- coalesce(clinical_sample$age.y,as.character(clinical_sample$age.x))
clinical_sample$age<- as.numeric(str_c(substr(clinical_sample$age,1,1),substr(clinical_sample$age,2+(mod(nchar(clinical_sample$age),4)==1),2+(mod(nchar(clinical_sample$age),4)==1)),".",
                 substr(clinical_sample$age,3+(mod(nchar(clinical_sample$age),4)!=3),3+(mod(nchar(clinical_sample$age),4)!=3))))


# --------------------------------------------------------------------------- ##

## systolic blood pressure
clinical_sample$sbp<-rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T)
publish(univariateTable(mobileUseNight ~ sbp,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ sbp,data=clinical_sample, column.percent=TRUE))

# diastolic blood pressure
clinical_sample$dbp<-rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T)
publish(univariateTable(mobileUseNight ~ dbp,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ dbp,data=clinical_sample, column.percent=TRUE))

## hip waist ratio
clinical_sample$ratiowaisthip <- as.numeric(clinical_sample$ratiowaisthip)
publish(univariateTable(mobileUseNight ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))

## bmi clinical
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable(mobileUseNight ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))


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
publish(univariateTable(mobileUseNight ~ hdl,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ hdl,data=clinical_sample, column.percent=TRUE))

## LDL
clinical_sample$ldl <- as.numeric(clinical_sample$ldl)
publish(univariateTable(mobileUseNight~ ldl,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y~ ldl,data=clinical_sample, column.percent=TRUE))

## VLDL
clinical_sample$vldl <- as.numeric(clinical_sample$vldl)
publish(univariateTable(mobileUseNight~ vldl,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y~ vldl,data=clinical_sample, column.percent=TRUE))

## total cholesterol
clinical_sample$t_cholesterol <- as.numeric(clinical_sample$t_cholesterol)
publish(univariateTable(mobileUseNight ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))

## triglycerides
clinical_sample$triglycerids <- as.numeric(clinical_sample$triglycerids)
publish(univariateTable(mobileUseNight ~ triglycerids,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ triglycerids,data=clinical_sample, column.percent=TRUE))

## hba1c
clinical_sample$hba1c <- as.numeric(clinical_sample$hba1c)
publish(univariateTable(mobileUseNight ~ hba1c,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(cluster.y ~ hba1c,data=clinical_sample, column.percent=TRUE))

# --------------------------------------------------------------------------- ##

#Models - multiple testing issue if we are going to 'pick and choose' which responses we would like to look at.
table(clinical_sample$age.x)

#Transforming to mids for modelling and inference
table(clinical_sample$description.y, clinical_sample$cluster.y )
clinical_sample$cluster.y <- factor(clinical_sample$cluster.y,levels = c("Cluster 3", "Cluster 2", "Cluster 4", "Cluster 1"))

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")


## Maximal posterior probability assignment: Six clusters (we use four clusters for clinical sample!)

dbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(dbp) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

vldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

df_ints_mpSix <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_mpSix) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")


## Maximal posterior probability assignment: Four clusters

dbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(dbp) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

vldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster.y+age+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

## estimater for four clusters i clinical sample
df_ints_mpFour <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_mpFour) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")
dt_ints_mpFour


## night-time smartphone use 

## Analyses: mobileUseNight and biomarkers (age = age.x?? - brug age.y i stedet = kliniske)

dbp_intsNight <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
glu_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
hba1c_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
hdl_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
ldl_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
t_chol_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
sbp_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
tri_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
vldl_intsNight <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
wh_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]
bmi_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseNight+age.x+educationW+occupationW,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]


df_intsNight <- data.frame(rbind(hdl_intsNight[2,],ldl_intsNight[2,],vldl_intsNight[2,],t_chol_intsNight[2,],tri_intsNight[2,],hba1c_intsNight[2,],dbp_intsNight[2,],sbp_intsNight[2,],wh_intsNight[2,],glu_intsNight[2,],bmi_intsNight[2,]),
                      rbind(hdl_intsNight[3,],ldl_intsNight[3,],vldl_intsNight[3,],t_chol_intsNight[3,],tri_intsNight[3,],hba1c_intsNight[3,],dbp_intsNight[3,],sbp_intsNight[3,],wh_intsNight[3,],glu_intsNight[3,],bmi_intsNight[3,]),
                      rbind(hdl_intsNight[4,],ldl_intsNight[4,],vldl_intsNight[4,],t_chol_intsNight[4,],tri_intsNight[4,],hba1c_intsNight[4,],dbp_intsNight[4,],sbp_intsNight[4,],wh_intsNight[4,],glu_intsNight[4,],bmi_intsNight[4,]))

colnames(df_intsNight) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsNight) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsNight[c(9,11,8,7,4,1,2,3,5,6),]
