###### Analysis of associations between metabolic biomarkers and smartphone activity ######


#This script contains analyses of the multiply imputed data sets from the SmartSleep project.
#The multiple imputation results are combined using Rubin's rules.

library(dplyr)
library(matrixStats)
library(mice)
library(MASS)
library(miceadds)
library(mitml)
library(dataReporter)

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

## load population sample
pop_data <-rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv"),imputation=imp_nr)

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
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")

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

#Base BMI models

summary(glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0),family=binomial))
summary(glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0),family=binomial))
summary(lm(log(log(bmi))~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0)))

plot(fitted(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))))
hist(residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),xlim=c(-20,20),breaks=200)
hist(simulate(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0)))$sim_1,breaks=40) #The bell-shape is not that well suited
hist(simulate(lm(log(bmi)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0)))$sim_1,breaks=40) #The bell-shape is not that well suited here either

#box_cox_transformation - bmi values are heavily right-skewed so we transform them closer toward normality ## hvordan kan jeg f? mening af dette?
bc <- boxcox(bmi ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))
(lambda <- bc$x[which.max(bc$y)])
new_model <- lm(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))

hist(((base_data$bmi^lambda-1)/lambda),breaks=40,xlim=c(0.78,0.81))
hist(simulate(new_model)$sim_1,breaks=40,xlim=c(0.78,0.81))

summary(new_model, conf.int = T)


## Using the mice package with mids objects
base_data_mids <- as.mids(base_data,.imp="imputation")
mod30<-(glm.mids((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data_mids,family=binomial))
mod30_p<-(glm.mids((bmi>=30)~(selfScoreCat+age+gender+education+occupation)+sample_weights,data=base_data_mids,family=binomial))
mod25<-(glm.mids((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data_mids,family=binomial))
modnum<-(lm.mids(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data_mids))

summary(pool(mod30),conf.int = T)
D1(mod30,mod30_p) #Test function
summary(pool(mod25))
summary(pool(modnum))
D1(mod30,mod30_p)
summary(pool(mod25), conf.int = T)
summary(pool(modnum), conf.int = T)


#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

bmi_followup <- rename(inner_join(rename(CSS,imputation=impnr),base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

bmi_followup$difference <- bmi_followup$bmi.fu-bmi_followup$bmi.base

bmi_followup$basebmi25=(bmi_followup$bmi.base>=25)
bmi_followup$basebmi30=(bmi_followup$bmi.base>=30)

hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))
hist(simulate(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup))$sim_1,xlim=c(-10,10)) #Too wide?
summary(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0)))

plot(fitted(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)),residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)))
hist(residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)),breaks=200,xlim=c(-10,10))

bmi_followup$bmi25change <- as.numeric((bmi_followup$bmi.fu>=25)!=(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeUp <- as.numeric((bmi_followup$bmi.fu>=25)>(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeDown <- as.numeric((bmi_followup$bmi.fu>=25)<(bmi_followup$bmi.base>=25))

bmi_followup$bmi30change <- as.numeric((bmi_followup$bmi.fu>=30)!=(bmi_followup$bmi.base>=30))
bmi_followup$bmi30changeUp <- as.numeric((bmi_followup$bmi.fu>=30)>(bmi_followup$bmi.base>=30))
bmi_followup$bmi30changeDown <- as.numeric((bmi_followup$bmi.fu>=30)<(bmi_followup$bmi.base>=30))

#The differences are not skewed, but their distribution is more narrow than a normal distribution - is this critical?


#MUsing the mids object for simple lm.
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

summary(pool(lm.mids(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids)))

plot(fitted(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0))),
     residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0))))
hist(residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0))),breaks=50)

#Alternative (better?) formulation
summary(pool(lm.mids(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids)))
plot(fitted(lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0))),
     residuals(lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0))))
hist(residuals(lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0))),breaks=50)


#Differences for indicators also (3 models per threshold: change, change from low to high group, and change from high to low group)

summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))

summary(glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))

#Using the mids object
summary(pool(glm.mids(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)))
summary(pool(glm.mids(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)))

summary(pool(glm.mids(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)), conf.int = T)
summary(pool(glm.mids(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)), conf.int = T)

summary(pool(glm.mids(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)), conf.int = T)
summary(pool(glm.mids(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)), conf.int = T)

summary(pool(glm.mids(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)), conf.int = T)
summary(pool(glm.mids(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)), conf.int = T)

#Alternative (better?) formulation with more easily interpretable parameters
summary(glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))

#And with the mids object class
summary(pool(glm.mids((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)))
summary(pool(glm.mids((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup_mids,family=binomial)))


#New idea: Try to make long format where followup and baseline are at different time points, and then make an interaction effect with time with bmi (indicators) as response.

long_data <- data.frame("bmi"=c(bmi_followup$bmi.base,bmi_followup$bmi.fu),"userid"=bmi_followup$userid,"sample_weights"=bmi_followup$sample_weights.y,"gender"=bmi_followup$gender.y,"age"=bmi_followup$age.y,
                        education=bmi_followup$education.y,occupation=bmi_followup$occupation.y,selfScoreCat = bmi_followup$selfScoreCat,"time"=c(rep(0,length(bmi_followup$bmi.base)),rep(1,length(bmi_followup$bmi.fu))),"imputation"=bmi_followup$imputation)

long_data_mids <- as.mids(long_data,.imp="imputation")

summary(pool(lm.mids(bmi~(selfScoreCat+age+gender+education+occupation)*time*sample_weights,data=long_data_mids)))
summary(pool(glm.mids((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*time*sample_weights,family=binomial,data=long_data_mids)))
summary(pool(glm.mids((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*time*sample_weights,family=binomial,data=long_data_mids)))



#####Tracking data for the followup CSS sample

## Hvad er de forskellige clusters? cluster 1 = night-time user?, cluster 2 =morning user?, cluster = non-user?, cluster 4 = evening user?
CSS_track$bmi[CSS_track$bmi==0]<-NA

#Models
summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial))
hist(residuals(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial)))
plot(fitted(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial)),residuals(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial)))

summary(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial))
hist(residuals(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial)))
plot(fitted(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial)),residuals(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0),family=binomial)))

summary(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0)))
hist(residuals(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))),breaks=40)
plot(fitted(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))),residuals(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))))

#box_cox_transformation
bc <- boxcox(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))
(lambda2 <- bc$x[which.max(bc$y)])
#new_model <- lm(((bmi^lambda-1)/lambda) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))

hist(((CSS_track$bmi^lambda2-1)/lambda2))
hist(simulate(new_model)$sim_1,breaks=20,xlim=c(0.830,0.860))
plot(fitted(new_model),residuals(new_model))


summary(new_model)


#Mice-based inference for three models:

CSS_track_mids<-as.mids(CSS_track,.imp="impnr",.id="userid")

summary(pool(lm.mids(((bmi^lambda2-1)/lambda2) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track_mids)))
summary(pool(glm.mids((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track_mids,family=binomial)))
summary(pool(glm.mids((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track_mids,family=binomial)))



#Tracking data: Population sample (random sample) - same analysis

hist(pop_track$bmi,breaks=50,xlim=c(0,50))

#analyses

bc <- boxcox(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(pop_track,imputation!=0))
(lambda3 <- bc$x[which.max(bc$y)])
new_model <- lm(((bmi^lambda3-1)/lambda3) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(pop_track,imputation!=0))

hist(((pop_track$bmi^lambda3-1)/lambda3))
hist(simulate(new_model)$sim_1,breaks=20)
plot(fitted(new_model),residuals(new_model))

pop_track$sample_weights<-as.numeric(pop_track$sample_weights)
pop_track_mids<-as.mids(pop_track,.imp="imputation",.id="userid")

## regression analysis of clusters of night-time smartphone use and overweight/obesity #justeres for selfScoreCat?? ## hvad er de forskellige clusters?? ## fortolkning?? ## inkluderer imp_nr=0?
## bmi kontinuert 
summary(pool(lm.mids(((bmi^lambda3-1)/lambda3) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=pop_track_mids)),conf.int=T)

## BMI >=25
summary(pool(glm.mids((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=pop_track_mids,family=binomial)),conf.int=T)
## no adjustment for selfScoreCat
summary(pool(glm.mids((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation)*sample_weights-sample_weights,data=pop_track_mids,family=binomial)),conf.int=T)



#Analysis of the clinical sample data - interest in biomarkers

clinical_sample <- rename(inner_join(clin_data,rename(clin_clinical,PNR=cpr),by="PNR"),bmi.self=bmi.x , bmi.clinical=bmi.y)
clinical_sample <- inner_join(clin_data,subject_tracking_clusters,by="userid")

clinical_sample$bmi[clinical_sample$bmi==0] <- NA
clinical_sample$bmi[clinical_sample$height<100 & clinical_sample$imputation!=0] <- (clinical_sample$weight/(((clinical_sample$height+100)/100)^2))[clinical_sample$height<100  & clinical_sample$imputation!=0]
clinical_sample$height[clinical_sample$height<100 & clinical_sample$imputation!=0] <- clinical_sample$height[clinical_sample$height<100 & clinical_sample$imputation!=0]+100 
clinical_sample$bmi[clinical_sample$height==clinical_sample$weight]<-NA
clinical_sample$bmi[clinical_sample$bmi<14]<-NA

clinical_sample$selfScore <- (clinical_sample$mobileUseBeforeSleep=="5-7 times per week")*4+(clinical_sample$mobileUseBeforeSleep=="2-4 times per week")*3+(clinical_sample$mobileUseBeforeSleep=="Once a week")*3+(clinical_sample$mobileUseBeforeSleep=="Every month or less")*2+(clinical_sample$mobileUseBeforeSleep=="Never")*1+
  (clinical_sample$mobileUseNight=="Every night or almost every night")*4+(clinical_sample$mobileUseNight=="A few times a week")*3+(clinical_sample$mobileUseNight=="A few times a month or less")*2+(clinical_sample$mobileUseNight=="Never")*1+
  (clinical_sample$mobileCheck==">20 times an hour")*4+(clinical_sample$mobileCheck=="11-20 times an hour")*4+(clinical_sample$mobileCheck=="5-10 times an hour")*3+(clinical_sample$mobileCheck=="1-4 times an hour")*2+(clinical_sample$mobileCheck=="Every 2nd hour")*2+(clinical_sample$mobileCheck=="Several times a day")*1+(clinical_sample$mobileCheck=="Once a day or less")*1+
  (clinical_sample$pmpuScale<=14)*1+(clinical_sample$pmpuScale>14 & clinical_sample$pmpuScale<17)*2+(clinical_sample$pmpuScale>=17 & clinical_sample$pmpuScale<19)*3+(clinical_sample$pmpuScale>=19)*4
summary(clinical_sample$selfScore[clinical_sample$impnr!=0])
clinical_sample$selfScoreCat <- NA
clinical_sample$selfScoreCat[!is.na(clinical_sample$selfScore)]<-"1"
clinical_sample$selfScoreCat[clinical_sample$selfScore>=8]="2"
clinical_sample$selfScoreCat[clinical_sample$selfScore>=10]="3"
clinical_sample$selfScoreCat[clinical_sample$selfScore>=12]="4"

#The subjects are scoring in the high end. Is this an issue or a characteristic of the data?

#Introducing interesting derived variables

clinical_sample$bmi <- as.numeric(clinical_sample$bmi)
clinical_sample$bmi25 <- as.numeric(clinical_sample$bmi>=25)
clinical_sample$bmi30 <- as.numeric(clinical_sample$bmi>=30)

#hdl, ldl, vldl, t_cholesterol, triglycerid, hba1c, (glucose), waist, hip, ratio waist hip, systolic bp og distolic bp 1-3: Ift. selvrapporteringer og tracking clusters

#Transforming to mids for modelling and inference

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")

## BMI >30
summary(pool(glm.mids((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=pop_track_mids,family=binomial)),conf.int=T)
## no adjustment for selfScoreCat

