
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

estimate.pooler <- function(coef,sd){
  n_row <- nrow(coef)
  n_col <- ncol(coef)
  
  coefs <- rowMeans(coef)
  sds <- sqrt(rowMeans(sd^2)+sum((coef-coefs)^2)/n_col)
  df <- data.frame("estimate"=coefs,"sd"=sds,"lower.CI"=coefs-1.96*sds,"upper.CI"=coefs+1.96*sds)
  return(df)
}

#Reading in the data
#setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
#load("H:/SmartSleep backup IT issues/imputation/myImputationCSS-res-00001.RData")
#CSS0 <- imp_CSS$data
#CSS0$impnr=0
#CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
#CSS <- rbind(CSS0,CSS)

## load tracking data 
subject_tracking_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_clusters.csv")

#base_weights <- read.csv2("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/SmartSleepExpWeighted.csv")
#load("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/SmartSleep Experiment/full_imp_base.RData")
#base_data <- full_imp_Base
#base_data <- rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv"),imputation=imp_nr)
#base_data$zipCode<-as.numeric(base_data$zipCode)

## load baseline data
base_data <- rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv"),imputation=imp_nr)
#base_data$sample_weights <- coalesce(base_data$sample_weights,0)

## load followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")
#CSS$sample_weights <- coalesce(CSS$sample_weights,0)

## load population sample
pop_data <-rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv"),imputation=imp_nr)
#pop_data$sample_weights <- coalesce(pop_data$sample_weights,0)

## load clinical data (survey and clinical data)
clin_data <- rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample/imp_clinical.csv"),imputation=impnr)
clin_clinical <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/SmartSleep Clinical/Data/Rådata/SmartSleepClinicalData.csv")

# --------------------------------------------------------------------------- ##

#Baseline data with self-reports

## if no mobile phone = NA
base_data$pmpuScale[base_data$mobilephone=="No mobile phone"] <- NA

## risk profiles for baseline data 
base_data$selfScore <- (base_data$mobileUseBeforeSleep=="5-7 times per week")*4+(base_data$mobileUseBeforeSleep=="2-4 times per week")*3+(base_data$mobileUseBeforeSleep=="Once a week")*3+(base_data$mobileUseBeforeSleep=="Every month or less")*2+(base_data$mobileUseBeforeSleep=="Never")*1+
  (base_data$mobileUseNight=="Every night or almost every night")*4+(base_data$mobileUseNight=="A few times a week")*3+(base_data$mobileUseNight=="A few times a month or less")*2+(base_data$mobileUseNight=="Never")*1+
  (base_data$mobileCheck==">20 times per hour")*4+(base_data$mobileCheck=="11-20 times per hour")*4+(base_data$mobileCheck=="5-10 times per hour")*3+(base_data$mobileCheck=="1-4 times per hour")*2+(base_data$mobileCheck=="Every 2nd hour")*2+(base_data$mobileCheck=="Several times per day")*1+(base_data$mobileCheck=="Once a day")*1+
  (base_data$pmpuScale<14)*1+(base_data$pmpuScale>=14 & base_data$pmpuScale<17)*2+(base_data$pmpuScale>=17 & base_data$pmpuScale<19)*3+(base_data$pmpuScale>=19)*4
summary(base_data$selfScore[base_data$imputation!=0])
base_data$selfScoreCat<-NA
base_data$selfScoreCat[!is.na(base_data$selfScore)] <- "1"
base_data$selfScoreCat[base_data$selfScore>=8]="2"
base_data$selfScoreCat[base_data$selfScore>=10]="3"
base_data$selfScoreCat[base_data$selfScore>=12]="4"
table(base_data$selfScoreCat[base_data$imputation!=0],useNA="always") 

## bmi for baseline data 
base_data$bmi[base_data$height<=100]<-round((base_data$weight/(((base_data$height+100)/100)^2))[base_data$height<=100],1)
base_data$height[base_data$height<=100] <- base_data$height[base_data$height<=100]+100
base_data$height[base_data$CS_ID==586] <- 100
base_data$bmi[base_data$CS_ID==586] <- 48.0
base_data$bmi[base_data$height==base_data$weight]<-NA
base_data$bmi[base_data$bmi<14]<-NA
base_data$bmi[base_data$bmi>147]<-NA
base_data$bmi[base_data$bmi==0]<-NA
summary(base_data$bmi[base_data$imputation!=0])

# --------------------------------------------------------------------------- ##
#Followup sample - using quartile levels from baseline sample

## risk profiles for CSS
CSS$selfScore <- (CSS$mobileUseBeforeSleep=="5-7 times per week")*4+(CSS$mobileUseBeforeSleep=="2-4 times per week")*3+(CSS$mobileUseBeforeSleep=="Once a week")*3+(CSS$mobileUseBeforeSleep=="Every month or less")*2+(CSS$mobileUseBeforeSleep=="Never")*1+
  (CSS$mobileUseNight=="Every night or almost every night")*4+(CSS$mobileUseNight=="A few times a week")*3+(CSS$mobileUseNight=="A few times a month or less")*2+(CSS$mobileUseNight=="Never")*1+
  (CSS$mobileCheck==">20 times an hour")*4+(CSS$mobileCheck=="11-20 times an hour")*4+(CSS$mobileCheck=="5-10 times an hour")*3+(CSS$mobileCheck=="1-4 times an hour")*2+(CSS$mobileCheck=="Every 2nd hour")*2+(CSS$mobileCheck=="Several times a day")*1+(CSS$mobileCheck=="Once a day or less")*1+
  (CSS$pmpuScale<=14)*1+(CSS$pmpuScale>14 & CSS$pmpuScale<17)*2+(CSS$pmpuScale>=17 & CSS$pmpuScale<19)*3+(CSS$pmpuScale>=19)*4
summary(CSS$selfScore[CSS$impnr!=0])
CSS$selfScoreCat <- NA
CSS$selfScoreCat[!is.na(CSS$selfScore)]<-"1"
CSS$selfScoreCat[CSS$selfScore>=8]="2"
CSS$selfScoreCat[CSS$selfScore>=10]="3"
CSS$selfScoreCat[CSS$selfScore>=12]="4"
table(CSS$selfScoreCat[CSS$impnr!=0], useNA="always")

## bmi CSS
CSS$bmi[CSS$bmi==0] <- NA
CSS$bmi[CSS$height<100 & CSS$impnr!=0] <- (CSS$weight/(((CSS$height+100)/100)^2))[CSS$height<100  & CSS$impnr!=0]
CSS$height[CSS$height<100 & CSS$impnr!=0] <- CSS$height[CSS$height<100 & CSS$impnr!=0]+100 
CSS$bmi[CSS$height==CSS$weight]<-NA
summary(CSS$bmi[CSS$impnr!=0])

## merge survey and tracking data 
CSS_track <- inner_join(CSS,subject_tracking_clusters,by="userid")
CSS_track <- rename(CSS_track,imputation=impnr)
CSS_track$sample_weights <- as.numeric(CSS_track$sample_weights)

# --------------------------------------------------------------------------- ##
#Population sample

## bmi
pop_data$bmi[pop_data$bmi==0] <- NA
pop_data$bmi[pop_data$height<100 & pop_data$imputation!=0] <- (pop_data$weight/(((pop_data$height+100)/100)^2))[pop_data$height<100  & pop_data$imputation!=0]
pop_data$height[pop_data$height<100 & pop_data$imputation!=0] <- pop_data$height[pop_data$height<100 & pop_data$imputation!=0]+100 
pop_data$bmi[pop_data$height==pop_data$weight]<-NA
pop_data$bmi[pop_data$bmi<14]<-NA
summary(pop_data$bmi[pop_data$imputation!=0])

## risk profiles for population sample
pop_data$selfScore <- (pop_data$mobileUseBeforeSleep=="5-7 times per week")*4+(pop_data$mobileUseBeforeSleep=="2-4 times per week")*3+(pop_data$mobileUseBeforeSleep=="Once a week")*3+(pop_data$mobileUseBeforeSleep=="Every month or less")*2+(pop_data$mobileUseBeforeSleep=="Never")*1+
  (pop_data$mobileUseNight=="Every night or almost every night")*4+(pop_data$mobileUseNight=="A few times a week")*3+(pop_data$mobileUseNight=="A few times a month or less")*2+(pop_data$mobileUseNight=="Never")*1+
  (pop_data$mobileCheck==">20 times an hour")*4+(pop_data$mobileCheck=="11-20 times an hour")*4+(pop_data$mobileCheck=="5-10 times an hour")*3+(pop_data$mobileCheck=="1-4 times an hour")*2+(pop_data$mobileCheck=="Every 2nd hour")*2+(pop_data$mobileCheck=="Several times a day")*1+(pop_data$mobileCheck=="Once a day or less")*1+
  (pop_data$pmpuScale<=14)*1+(pop_data$pmpuScale>14 & pop_data$pmpuScale<17)*2+(pop_data$pmpuScale>=17 & pop_data$pmpuScale<19)*3+(pop_data$pmpuScale>=19)*4
summary(pop_data$selfScore[pop_data$imputation!=0])
pop_data$selfScoreCat <- NA
pop_data$selfScoreCat[!is.na(pop_data$selfScore)]<-"1"
pop_data$selfScoreCat[pop_data$selfScore>=8]="2"
pop_data$selfScoreCat[pop_data$selfScore>=10]="3"
pop_data$selfScoreCat[pop_data$selfScore>=12]="4"
table(pop_data$selfScoreCat[pop_data$imputation!=0])

## merge tracking and survey data for population sample
pop_track <- inner_join(pop_data,subject_tracking_clusters,by="userid")

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

#### BASE POPULATION

#Base BMI models (test til om modellen kører)

summary(glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0),family=binomial))
summary(glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0),family=binomial))
summary(lm(log(log(bmi))~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0)))

## test kontinuert bmi
plot(fitted(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0))),residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0))))
hist(residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0))),xlim=c(-20,20),breaks=200)
hist(simulate(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation!=0)))$sim_1,breaks=40) #The bell-shape is not that well suited


#Better alternative: Pretty good fit. A generalized family of models.

m <- gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==1)),family = BCCG)
#m_sum <- summary(m)
plot(m)
coef(m)
diag(vcov(m))
confint(m)

#m_simple <- gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==1)), weights=sample_weights,family = BCCG)
#plot(m_simple)

#Works with mids?

m_coefs <- m_sds <- matrix(nrow=length(coef(m)),ncol=20)

for (i in 1:20){
  m <- gamlss(bmi ~ (selfScoreCat+age+gender+education+occupation), sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)), family = BCCG)
  m_coefs[,i] <- coef(m)
  m_sds[,i] <- diag(vcov(m))[1:length(coef(m))] #or use rvcov for robust standard errors - if the variance model is thought to be incorrect
  
  #m_list[[i]] <- m
  #m_sum_list[[i]] <- m_sum
}

#Now we may simply combine results using Rubin's rules.

RRresult <- estimate.pooler(m_coefs,m_sds)
rownames(RRresult) <- names(coef(gamlss(bmi ~ (selfScoreCat+age+gender+education+occupation), sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==1)), family = BCCG)))
RRresult$p.value <- pnorm(q=0,mean=RRresult$estimate,sd=RRresult$sd)




## Using the mice package with mids objects
base_data_mids <- as.mids(base_data,.imp="imputation")
#<<<<<<< Updated upstream
mod30 <- with(base_data_mids,glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod25 <- with(base_data_mids,glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))

#=======
mod30 <- (glm.mids((bmi>=30)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=base_data_mids,family=binomial))
mod30_p <- (glm.mids((bmi>=30)~(selfScoreCat+age+gender+education+occupation)+sample_weights, weights=sample_weights, data=base_data_mids,family=binomial))
mod25 <- (glm.mids((bmi>=25)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=base_data_mids,family=binomial))
modnum <- (lm.mids(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=base_data_mids))
#>>>>>>> Stashed changes

## OR for BMI>30
model30 <- summary(pool(mod30),conf.int = T)
exp(model30$estimate)
exp(model30$`2.5 %`)
exp(model30$`97.5 %`)

D1(mod30,mod30_p) #Test function
summary(pool(mod25))
summary(pool(modnum))
D1(mod30,mod30_p)
summary(pool(mod25), conf.int = T)
summary(pool(modnum), conf.int = T)

## OR for BMI >25
model25 <- summary(pool(mod25), conf.int=T)
exp(model25$estimate)
exp(model25$`2.5 %`)
exp(model25$`97.5 %`)

#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

bmi_followup <- rename(inner_join(rename(CSS,imputation=impnr),base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

## difference mellem follow-up og baseline
bmi_followup$difference <- bmi_followup$bmi.fu-bmi_followup$bmi.base
mean(bmi_followup$difference[!is.na(bmi_followup$difference)])

bmi_followup$basebmi25=(bmi_followup$bmi.base>=25)
bmi_followup$basebmi30=(bmi_followup$bmi.base>=30)

hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))

## ændringer i bmi ja eller nej
bmi_followup$bmi25change <- as.numeric((bmi_followup$bmi.fu>=25)!=(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeUp <- as.numeric((bmi_followup$bmi.fu>=25)>(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeDown <- as.numeric((bmi_followup$bmi.fu>=25)<(bmi_followup$bmi.base>=25))

bmi_followup$bmi30change <- as.numeric((bmi_followup$bmi.fu>=30)!=(bmi_followup$bmi.base>=30))
bmi_followup$bmi30changeUp <- as.numeric((bmi_followup$bmi.fu>=30)>(bmi_followup$bmi.base>=30))
bmi_followup$bmi30changeDown <- as.numeric((bmi_followup$bmi.fu>=30)<(bmi_followup$bmi.base>=30))

#The differences are not skewed, but their distribution is more narrow than a normal distribution - is this critical?

#MUsing the mids object for simple lm. (for differencen)
bmi_followup$sample_weights <- bmi_followup$sample_weights.y
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

summary(pool(with(bmi_followup_mids,lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y))))

plot(fitted(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))),
     residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))))
hist(residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))),breaks=50)

#Alternative (better?) formulation

##ændringer i bmi over tid (men residual plottet siger at modellen er centreret omkring middelværdien)
summary(pool(with(bmi_followup_mids,lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y))), conf.int = T)

## residual plot
plot(fitted(lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))),
     residuals(lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))))
hist(residuals(lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))),breaks=50)

#Differences for indicators also (3 models per threshold: change, change from low to high group, and change from high to low group)

summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))
summary(glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))

summary(glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))
summary(glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))

summary(glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))
summary(glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))

#Using the mids object
## change bmi 25
#summary(pool(with(bmi_followup_mids,glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))))
#summary(pool(with(bmi_followup_mids,glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))))
summary(pool(with(bmi_followup_mids,glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)
summary(pool(with(bmi_followup_mids,glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

## change from low to high group
bmi_followup_mids <- as.mids(subset(bmi_followup,bmi.base<25),.imp="imputation")
summary(pool(with(bmi_followup_mids,glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)
bmi_followup_mids <- as.mids(subset(bmi_followup,bmi.base<30),.imp="imputation")
summary(pool(with(bmi_followup_mids,glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

#change from high to low group
bmi_followup_mids <- as.mids(subset(bmi_followup,bmi.base>=25),.imp="imputation")
summary(pool(with(bmi_followup_mids,glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)
bmi_followup_mids <- as.mids(subset(bmi_followup,bmi.base>=30),.imp="imputation")
summary(pool(with(bmi_followup_mids,glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

#Alternative (better?) formulation with more easily interpretable parameters
summary(glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))
summary(glm((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))

#And with the mids object class
 ## change 25 
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")
change25 <- with(bmi_followup_mids,glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))
modelchange25 <- summary(pool(change25), conf.int = T)
exp(modelchange25$estimate)
exp(modelchange25$`2.5 %`)
exp(modelchange25$`97.5 %`)
##summary(pool(with(bmi_followup_mids,glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))))

##change 30
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")
change30 <- with(bmi_followup_mids,glm((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))
modelchange30 <- summary(pool(change30), conf.int = T)
exp(modelchange30$estimate)
exp(modelchange30$`2.5 %`)
exp(modelchange30$`97.5 %`)

##OLD: summary(pool(with(bmi_followup_mids,glm((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))))

#New idea: Try to make long format where followup and baseline are at different time points, and then make an interaction effect with time with bmi (indicators) as response.

long_data <- data.frame("bmi"=c(bmi_followup$bmi.base,bmi_followup$bmi.fu),"userid"=bmi_followup$userid,"sample_weights"=bmi_followup$sample_weights.y,"gender"=bmi_followup$gender.y,"age"=bmi_followup$age.y,
                        education=bmi_followup$education.y,occupation=bmi_followup$occupation.y,selfScoreCat = bmi_followup$selfScoreCat.y,"time"=c(rep(0,length(bmi_followup$bmi.base)),rep(1,length(bmi_followup$bmi.fu))),"imputation"=bmi_followup$imputation)

long_data_mids <- as.mids(long_data,.imp="imputation")

#summary(pool(lm.mids(bmi~(selfScoreCat+age+gender+education+occupation)*time, weights=sample_weights, data=long_data_mids)))
summary(pool(with(long_data_mids,glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*time,family=binomial, weights=sample_weights))))
summary(pool(with(long_data_mids,glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*time,family=binomial, weights=sample_weights))))

#Now with the gamlss for numeric bmi

m <- gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,
            nu.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,weights=sample_weights, data=na.omit(subset(long_data,imputation==10)),family=BCCG,method=RS(100))

m_coefs_long <- m_sds_long <- matrix(nrow=length(coef(m)),ncol=20)

for (i in 1:20){
  m <- gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,
              nu.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,weights=sample_weights, data=na.omit(subset(long_data,imputation==i)), family = BCCG,method=RS(100))
  m_coefs_long[,i] <- coef(m)
  m_sds_long[,i] <- diag(vcov(m))[1:length(coef(m))] #or use rvcov for robust standard errors - if the variance model is thought to be incorrect
  
  #m_list[[i]] <- m
  #m_sum_list[[i]] <- m_sum
}

#Now we may simply combine results using Rubin's rules.

RRresult_long <- estimate.pooler(m_coefs_long,m_sds_long)
rownames(RRresult_long) <- names(coef(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,
                                        nu.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,weights=sample_weights, data=na.omit(subset(long_data,imputation==10)), family = BCCG)))

RRresult_long$p.value <- pnorm(q=0,mean=RRresult_long$estimate,sd=RRresult_long$sd)





#####Tracking data for the followup CSS sample

## Hvad er de forskellige clusters? cluster 1 = night-time user?, cluster 2 =morning user?, cluster = non-user?, cluster 4 = evening user?
CSS_track$bmi[CSS_track$bmi==0]<-NA

#Models
summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial))
hist(residuals(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial)))
plot(fitted(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial)),residuals(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial)))

summary(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial))
hist(residuals(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial)))
plot(fitted(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial)),residuals(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0),family=binomial)))

summary(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0)))
hist(residuals(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0))),breaks=40)
plot(fitted(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0))),residuals(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(CSS_track,imputation!=0))))


#Mice-based inference for three models:

CSS_track_mids<-as.mids(CSS_track,.imp="imputation",.id="userid")

summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))))
summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))))

m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
            nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100))
#robust=TRUE? Doesn't change much when the fit is this good.

m_coefs_CSS <- m_sds_CSS <- matrix(nrow=length(coef(m)),ncol=20)

for (i in 1:20){
  m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
              nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==i)), family = BCCG,method=RS(100))
  m_coefs_CSS[,i] <- coef(m)
  m_sds_CSS[,i] <- diag(vcov(m))[1:length(coef(m))] #or use rvcov for robust standard errors - if the variance model is thought to be incorrect
  
  #m_list[[i]] <- m
  #m_sum_list[[i]] <- m_sum
}

#Now we may simply combine results using Rubin's rules.

RRresult_CSStrack <- estimate.pooler(m_coefs_CSS,m_sds_CSS)
rownames(RRresult_CSStrack) <- names(coef(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
                                                 nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)), family = BCCG)))
RRresult_CSStrack$p.value <- pnorm(q=0,mean=RRresult_CSStrack$estimate,sd=RRresult_CSStrack$sd)








#Tracking data: Population sample (random sample) - same analysis

hist(pop_track$bmi,breaks=50,xlim=c(0,50))


#analyses

pop_track$sample_weights<-as.numeric(pop_track$sample_weights)
pop_track_mids<-as.mids(pop_track,.imp="imputation",.id="userid")

## regression analysis of clusters of night-time smartphone use and overweight/obesity #justeres for selfScoreCat?? ## hvad er de forskellige clusters?? ## fortolkning?? ## inkluderer imp_nr=0?
## bmi kontinuert 


m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
            nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100))

m_coefs_pop <- m_sds_pop <- matrix(nrow=length(coef(m)),ncol=20)

for (i in 1:20){
  m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
              nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==i)), family = BCCG,method=RS(100))
  m_coefs_pop[,i] <- coef(m)
  m_sds_pop[,i] <- diag(vcov(m))[1:length(coef(m))] #or use rvcov for robust standard errors - if the variance model is thought to be incorrect
  
  #m_list[[i]] <- m
  #m_sum_list[[i]] <- m_sum
}

#Now we may simply combine results using Rubin's rules.

RRresult_poptrack <- estimate.pooler(m_coefs_pop,m_sds_pop)
rownames(RRresult_poptrack) <- names(coef(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
                                                 nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)), family = BCCG)))
RRresult_poptrack$p.value <- pnorm(q=0,mean=RRresult_poptrack$estimate,sd=RRresult_poptrack$sd)


## BMI >=25
#summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25 <- with(pop_track_mids,glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25 <- summary(pool(Random25), conf.int=T)
exp(modelRandom25$estimate)
exp(modelRandom25$`2.5 %`)
exp(modelRandom25$`97.5 %`)

## no adjustment for selfScoreCat
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int=T)
exp(modelRandom25No$estimate)
exp(modelRandom25No$`2.5 %`)
exp(modelRandom25No$`97.5 %`)

## BMI >30
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30 <- with(pop_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30 <- summary(pool(Random30), conf.int = T)
exp(modelRandom30$estimate)
exp(modelRandom30$`2.5 %`)
exp(modelRandom30$`97.5 %`)

## no adjustment for selfScoreCat
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
exp(modelRandom30No$estimate)
exp(modelRandom30No$`2.5 %`)
exp(modelRandom30No$`97.5 %`)


####### BMI and risk profiles

## bmi >25
#summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25Risk <- with(pop_track_mids,glm((bmi>=25) ~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25Risk <- summary(pool(Random25Risk), conf.int = T)
exp(modelRandom25Risk$estimate)
exp(modelRandom25Risk$`2.5 %`)
exp(modelRandom25Risk$`97.5 %`)


## BMI >30
#summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30Risk <- with(pop_track_mids,glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30Risk <- summary(pool(Random30Risk), conf.int = T)
exp(modelRandom30Risk$estimate)
exp(modelRandom30Risk$`2.5 %`)
exp(modelRandom30Risk$`97.5 %`)

###############################################################################
###############################################################################
###############################################################################

#Analysis of the clinical sample data - interest in biomarkers

## merge survey and clinical data
clinical_sample <- rename(inner_join(clin_data,rename(clin_clinical,PNR=cpr),by="PNR"),bmi.self=bmi.x , bmi.clinical=bmi.y)
## merge with tracking data
clinical_sample <- inner_join(clinical_sample,subject_tracking_clusters,by="userid")

## no cleanng in bmi.clinical as it should be correct compared to self-reports
#clinical_sample$bmi.clinical[clinical_sample$bmi.clinical==0] <- NA
#clinical_sample$bmi.clinical[clinical_sample$height<100 & clinical_sample$imputation!=0] <- (clinical_sample$weight/(((clinical_sample$height+100)/100)^2))[clinical_sample$height<100  & clinical_sample$imputation!=0]
#clinical_sample$height[clinical_sample$height<100 & clinical_sample$imputation!=0] <- clinical_sample$height[clinical_sample$height<100 & clinical_sample$imputation!=0]+100 
#clinical_sample$bmi.clinical[clinical_sample$height==clinical_sample$weight]<-NA#
#clinical_sample$bmi.clinical[clinical_sample$bmi.clinical<14]<-NA

clinical_sample$selfScore <- (clinical_sample$mobileUseBeforeSleep=="5-7 times per week")*4+(clinical_sample$mobileUseBeforeSleep=="2-4 times per week")*3+(clinical_sample$mobileUseBeforeSleep=="Once a week")*3+(clinical_sample$mobileUseBeforeSleep=="Every month or less")*2+(clinical_sample$mobileUseBeforeSleep=="Never")*1+
  (clinical_sample$mobileUseNight=="Every night or almost every night")*4+(clinical_sample$mobileUseNight=="A few times a week")*3+(clinical_sample$mobileUseNight=="A few times a month or less")*2+(clinical_sample$mobileUseNight=="Never")*1+
  (clinical_sample$mobileCheck==">20 times an hour")*4+(clinical_sample$mobileCheck=="11-20 times an hour")*4+(clinical_sample$mobileCheck=="5-10 times an hour")*3+(clinical_sample$mobileCheck=="1-4 times an hour")*2+(clinical_sample$mobileCheck=="Every 2nd hour")*2+(clinical_sample$mobileCheck=="Several times a day")*1+(clinical_sample$mobileCheck=="Once a day or less")*1+
  (clinical_sample$pmpuScale<=14)*1+(clinical_sample$pmpuScale>14 & clinical_sample$pmpuScale<17)*2+(clinical_sample$pmpuScale>=17 & clinical_sample$pmpuScale<19)*3+(clinical_sample$pmpuScale>=19)*4
summary(clinical_sample$selfScore[clinical_sample$imputation!=0])
clinical_sample$selfScoreCat <- NA
clinical_sample$selfScoreCat[!is.na(clinical_sample$selfScore)]<-"1"
clinical_sample$selfScoreCat[clinical_sample$selfScore>=8]="2"
clinical_sample$selfScoreCat[clinical_sample$selfScore>=10]="3"
clinical_sample$selfScoreCat[clinical_sample$selfScore>=12]="4"
table(clinical_sample$selfScoreCat, useNA="always")

## descriptives of clinical sample
publish(univariateTable(selfScoreCat ~ age.x,data=clinical_sample, column.percent=TRUE))

publish(univariateTable(selfScoreCat ~ bmi,data=clinical_sample, column.percent=TRUE))

#The subjects are scoring in the high end. Is this an issue or a characteristic of the data?

#Introducing interesting derived variables

clinical_sample$bmi <- as.numeric(clinical_sample$bmi.clinical)
clinical_sample$bmi25 <- as.numeric(clinical_sample$bmi.clinical>=25)
clinical_sample$bmi30 <- as.numeric(clinical_sample$bmi.clinical>=30)

clinical_sample$age<- as.numeric(str_c(substr(clinical_sample$age.y,1,1),substr(clinical_sample$age.y,2+(mod(nchar(clinical_sample$age.y),4)==1),2+(mod(nchar(clinical_sample$age.y),4)==1)),".",
                 substr(clinical_sample$age.y,3+(mod(nchar(clinical_sample$age.y),4)!=3),3+(mod(nchar(clinical_sample$age.y),4)!=3))))

## systolic blood pressure
clinical_sample$sbp<-rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ sbp,data=clinical_sample, column.percent=TRUE))

## categorize systolic blood pressure into high (>=140) and normal (<140)
clinical_sample$sbpCat[clinical_sample$sbp>=140] <- "High"
clinical_sample$sbpCat[clinical_sample$sbp<140] <- "Normal"
publish(univariateTable( ~ sbpCat,data=clinical_sample, column.percent=TRUE))

# diastolic blood pressure
clinical_sample$dbp<-rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ dbp,data=clinical_sample, column.percent=TRUE))

clinical_sample$dbpCat[clinical_sample$dbp>=90] <- "High"
clinical_sample$dbpCat[clinical_sample$dbp<90] <- "Normal"
table(clinical_sample$dbpCat, useNA="always")

## hip waist ratio

publish(univariateTable(selfScoreCat ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))

clinical_sample$ratiowaisthip <- as.numeric(clinical_sample$ratiowaisthip)
clinical_sample$ratiowaisthipCat[clinical_sample$ratiowaisthip>=0.85] <- "High"
clinical_sample$ratiowaisthipCat[clinical_sample$ratiowaisthip<0.85] <- "Normal"
table(clinical_sample$ratiowaisthipCat)

#hdl, ldl, vldl, t_cholesterol, triglycerid, hba1c, (glucose), waist, hip, ratio waist hip, systolic bp og distolic bp 1-3: Ift. selvrapporteringer og tracking clusters

#Transforming to mids for modelling and inference

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")

#Looks into data

hist(as.numeric(clinical_sample$hdl),breaks=20)
hist(as.numeric(clinical_sample$ldl),breaks=20)
hist(as.numeric(clinical_sample$vldl),breaks=20)
hist(as.numeric(clinical_sample$triglycerids),breaks=40)
hist(as.numeric(clinical_sample$hba1c),breaks=40)
hist(as.numeric(clinical_sample$glucose),breaks=40)
hist(as.numeric(clinical_sample$waist),breaks=20)
hist(as.numeric(clinical_sample$hip),breaks=20)
hist(as.numeric(clinical_sample$ratiowaisthip),breaks=20)
hist(rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T),breaks=20)
hist(rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T),breaks=20)

clinical_sample$hdl <- as.numeric(clinical_sample$hdl)
publish(univariateTable(selfScoreCat ~ hdl,data=clinical_sample, column.percent=TRUE))

clinical_sample$ldl <- as.numeric(clinical_sample$ldl)
publish(univariateTable( selfScoreCat~ ldl,data=clinical_sample, column.percent=TRUE))

clinical_sample$triglycerids <- as.numeric(clinical_sample$triglycerids)
publish(univariateTable(selfScoreCat ~ triglycerids,data=clinical_sample, column.percent=TRUE))

clinical_sample$hba1c <- as.numeric(clinical_sample$hba1c)
publish(univariateTable(selfScoreCat ~ vldl,data=clinical_sample, column.percent=TRUE))

clinical_sample$vldl <- as.numeric(clinical_sample$vldl)
publish(univariateTable(selfScoreCat ~ hba1c,data=clinical_sample, column.percent=TRUE))

#Models - multiple testing issue if we are going to 'pick and collect' which responses we would like to look at.
table(clinical_sample$age.x)

#hdl
summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation))))
hist(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
lines(seq(from=min(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),na.rm=T),length.out=100),dnorm(x=seq(from=min(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),na.rm=T),length.out=100),mean=mean(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))),sd=sd(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))))
plot(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#ldl
summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))
hist(residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res=residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#vldl
summary(pool(gwith(data=clinical_mids, lm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,family=Gamma))))
hist(residuals(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#triglycerids
summary(pool(gwith(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,family=Gamma))))
hist(residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma)))
plot(fitted(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma)),residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),family=Gamma)))

#hba1c
summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))
hist(residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#glucose
summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))
hist(residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#ratiowaisthip
summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))
hist(residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#sbp
summary(pool(with(data=clinical_mids, lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))
hist(residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#dbp
summary(pool(with(data=clinical_mids, lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))
hist(residuals(lm(dpb ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(dpb ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(dpb ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(dpb ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(dpb ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

#Generally the residuals look reasonably centered, with a few positive outliers. The residual distributions on the first distribution actually look reasonably normal, save for the few (extreme) outliers.

#Seems that these models are appropriate, and that normal approximations of confidence interval will be reasonable too. We can however also just use the profile likelihood CI's.

