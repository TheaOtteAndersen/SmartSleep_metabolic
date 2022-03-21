
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
subject_tracking_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_clusters.csv")

## load baseline data
base_data <- rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv"),imputation=imp_nr)

## load followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")

## load population sample
pop_data <-rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv"),imputation=imp_nr)

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

# bar chart selfScoreCat 
ggplot(base_data, aes(x = factor(selfScoreCat))) +
  geom_bar()

## tjek risk profiles
publish(univariateTable( ~ mobileUseBeforeSleep,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ mobileCheck,data=base_data, column.percent=TRUE))
table((base_data$pmpuScale<14)*1+(base_data$pmpuScale>=14 & base_data$pmpuScale<17)*2+(base_data$pmpuScale>=17 & base_data$pmpuScale<19)*3+(base_data$pmpuScale>=19)*4)/20

publish(univariateTable( ~ selfScoreCat,data=base_data, column.percent=TRUE))

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

base_data_mids <- as.mids(base_data,.imp="imputation")

## BMI kategoriseringer ved baseline
table(base_data$bmi<25, base_data$selfScoreCat)/21
table(base_data$bmi>=25&base_data$bmi<30, base_data$selfScoreCat)/21
table(base_data$bmi>=30, base_data$selfScoreCat)/21

#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

bmi_followup <- rename(inner_join(rename(CSS,imputation=impnr),base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

## difference mellem follow-up og baseline
bmi_followup$difference <- bmi_followup$bmi.fu-bmi_followup$bmi.base
mean(bmi_followup$difference[!is.na(bmi_followup$difference)])

bmi_followup$basebmi25=(bmi_followup$bmi.base>=25)
bmi_followup$basebmi30=(bmi_followup$bmi.base>=30)

bmi_followup$sample_weights <- bmi_followup$sample_weights.y
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

#New idea: Try to make long format where followup and baseline are at different time points, and then make an interaction effect with time with bmi (indicators) as response.

long_data <- data.frame("bmi"=c(bmi_followup$bmi.base,bmi_followup$bmi.fu),"userid"=bmi_followup$userid,"sample_weights"=bmi_followup$sample_weights.y,"gender"=bmi_followup$gender.y,"age"=bmi_followup$age.y,
                        education=bmi_followup$education.y,occupation=bmi_followup$occupation.y,selfScoreCat = bmi_followup$selfScoreCat.y,"time"=c(rep(0,length(bmi_followup$bmi.base)),rep(1,length(bmi_followup$bmi.fu))),"imputation"=bmi_followup$imputation)

long_data_mids <- as.mids(long_data,.imp="imputation")



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

# bar chart selfScoreCat 
ggplot(CSS, aes(x = factor(selfScoreCat))) +
  geom_bar()

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

CSS_track_mids<-as.mids(CSS_track,.imp="imputation",.id="userid")

# --------------------------------------------------------------------------- ##
#Merging base and followup


#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

bmi_followup <- rename(inner_join(rename(CSS,imputation=impnr),base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

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
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

#New idea: Try to make long format where followup and baseline are at different time points, and then make an interaction effect with time with bmi (indicators) as response.

long_data <- data.frame("bmi"=c(bmi_followup$bmi.base,bmi_followup$bmi.fu),"userid"=bmi_followup$userid,"sample_weights"=bmi_followup$sample_weights.y,"gender"=bmi_followup$gender.y,"age"=bmi_followup$age.y,
                        education=bmi_followup$education.y,occupation=bmi_followup$occupation.y,selfScoreCat = bmi_followup$selfScoreCat.y,"time"=c(rep(0,length(bmi_followup$bmi.base)),rep(1,length(bmi_followup$bmi.fu))),"imputation"=bmi_followup$imputation)

long_data_mids <- as.mids(long_data,.imp="imputation")



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
pop_track$sample_weights<-as.numeric(pop_track$sample_weights)
pop_track_mids<-as.mids(pop_track,.imp="imputation",.id="userid")

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


#Checking the validity of Wald intervals

coef <- sds <- matrix(nrow=length(c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)),ncol=20)
for (i in 1:20){
  m <- gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  coef[,i] <- c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)
  sds[,i] <- sqrt(diag(vcov(m)))
}

ests <- estimate.pooler(coef,sds)[,1]

logL1 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==1)),family = BCCG))
logL2 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==2)),family = BCCG))
logL3 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==3)),family = BCCG))
logL4 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==4)),family = BCCG))
logL5 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==5)),family = BCCG))
logL6 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==6)),family = BCCG))
logL7 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==7)),family = BCCG))
logL8 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==8)),family = BCCG))
logL9 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==9)),family = BCCG))
logL10 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==10)),family = BCCG))
logL11 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==11)),family = BCCG))
logL12 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==12)),family = BCCG))
logL13 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==13)),family = BCCG))
logL14 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==14)),family = BCCG))
logL15 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==15)),family = BCCG))
logL16 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==16)),family = BCCG))
logL17 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==17)),family = BCCG))
logL18 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==18)),family = BCCG))
logL19 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==19)),family = BCCG))
logL20 <- gen.likelihood(gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~(selfScoreCat+age+gender+education+occupation), nu.formula =~ (selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==20)),family = BCCG))


logLhat <- (logL1(ests)+logL2(ests)+logL3(ests)+logL4(ests)+logL5(ests)+
  logL6(ests)+logL7(ests)+logL8(ests)+logL9(ests)+logL10(ests)+
  logL11(ests)+logL12(ests)+logL13(ests)+logL14(ests)+logL15(ests)+
  logL16(ests)+logL17(ests)+logL18(ests)+logL19(ests)+logL20(ests))/20

#Would ideally put the pooled estimates in here, and use them for generating profile likelihood intervals in the AVERAGE likelihood function across the imputations... Hence we would need to also generate that function by generating each of the likelihoods and making the average of the 20 likelihoods evaluated in given parameters.
#hatmucoefs <- c() #This should then contain the pooled estimates to put into the 20 likelihoods.

change_seq <- change_seq1 <- seq(from=-2.5,to=2.5,by=0.01)
change_seq2 <- seq(from=-0.1,to=0.1,by=0.0001)
out_seq <- numeric(0)
out_ints <- matrix(0,nrow=length(m$mu.coefficients),ncol=2)

for (k in 1:4){
  for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
    #if (k==1){
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)]))))) #c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}

out_seq <- numeric(0)
k=5
for (i in 1:(length(change_seq2))){ #*(k==1)+length(change_seq2)*(k>1)
  #if (k==1){
  out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL2(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL3(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL4(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL5(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL6(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL7(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL8(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL9(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL10(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL11(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL12(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL13(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL14(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL15(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL16(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL17(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL18(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL19(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL20(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])))))
  #}
  #if (k>1){
  #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
  #}
}
out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
plot(m$mu.coefficients[k]+change_seq2,out_seq)

out_seq <- numeric(0)
for (k in 6:(length(m$mu.coefficients))){
  for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
    #if (k==1){
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])))))
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}

#The resulting plots should be quadratically shaped around the MLE, if we would like to use the Wald approximation.
#The profile likelihood intervals seem quite narrow...



# --------------------------------------------------------------------------- ##
## cross-sectional associations between risk profiles and bmi (25, 30 & continous)
# --------------------------------------------------------------------------- ##


## Using the mice package with mids objects
mod30 <- with(base_data_mids,glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod25 <- with(base_data_mids,glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))

## test for trend
TEST <- with(base_data_mids,glm((bmi>=30)~((as.numeric(selfScoreCat))+age+gender+education+occupation), weights=sample_weights,family=binomial))
testT <- summary(pool(TEST), conf.int = T)

TEST2 <- with(base_data_mids,glm((bmi>=25)~((as.numeric(selfScoreCat))+age+gender+education+occupation), weights=sample_weights,family=binomial))
test2 <- summary(pool(TEST2), conf.int = T)


## OR for BMI>30
model30 <- summary(pool(mod30),conf.int = T)
exp(model30$estimate)
exp(model30$`2.5 %`)
exp(model30$`97.5 %`)

## OR for BMI >25
model25 <- summary(pool(mod25), conf.int=T)
exp(model25$estimate)
exp(model25$`2.5 %`)
exp(model25$`97.5 %`)

# --------------------------------------------------------------------------- ##
## longitudinal analysis of risk scores of smartphone behavior and changes in BMI
# --------------------------------------------------------------------------- ##
#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))

## ændringer i bmi ja eller nej
table(bmi_followup$bmi25change)
bmi_followup$bmi25change <- as.numeric((bmi_followup$bmi.fu>=25)!=(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeUp <- as.numeric((bmi_followup$bmi.fu>=25)>(bmi_followup$bmi.base>=25))
bmi_followup$bmi25changeDown <- as.numeric((bmi_followup$bmi.fu>=25)<(bmi_followup$bmi.base>=25))

#bmi_followup$bmi30change <- as.numeric((bmi_followup$bmi.fu>=30)!=(bmi_followup$bmi.base>=30))
#bmi_followup$bmi30changeUp <- as.numeric((bmi_followup$bmi.fu>=30)>(bmi_followup$bmi.base>=30))
#bmi_followup$bmi30changeDown <- as.numeric((bmi_followup$bmi.fu>=30)<(bmi_followup$bmi.base>=30))

#The differences are not skewed, but their distribution is more narrow than a normal distribution - is this critical?

#MUsing the mids object for simple lm. (for differencen)

#summary(pool(with(bmi_followup_mids,lm(difference~(selfScoreCat+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y))))

plot(fitted(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))),
     residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))))
hist(residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y, data=subset(bmi_followup,imputation!=0))),breaks=50)

#Alternative (better?) formulation

##ændringer i bmi over tid (men residual plottet siger at modellen er centreret omkring middelværdien)
#summary(pool(with(bmi_followup_mids,lm(bmi.fu~(bmi.base+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y))), conf.int = T)

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
#summary(pool(with(bmi_followup_mids,glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)
#summary(pool(with(bmi_followup_mids,glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

## change from low to high group (korrekte!)
#bmi_followup_mids_25risk <- as.mids(subset(bmi_followup,bmi.base<25),.imp="imputation")
#summary(pool(with(bmi_followup_mids_25risk,glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

#bmi_followup_mids_30risk <- as.mids(subset(bmi_followup,bmi.base<30),.imp="imputation")
#summary(pool(with(bmi_followup_mids,glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

#change from high to low group
#bmi_followup_mids_25case <- as.mids(subset(bmi_followup,bmi.base>=25),.imp="imputation")
#summary(pool(with(bmi_followup_mids_25case,glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)
#bmi_followup_mids_30case <- as.mids(subset(bmi_followup,bmi.base>=30),.imp="imputation")
#summary(pool(with(bmi_followup_mids_30case,glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

#Alternative (better?) formulation with more easily interpretable parameters
#summary(glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))
#summary(glm((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights, data=bmi_followup,family=binomial))

#And with the mids object class
 ## change 25 
#change25 <- with(bmi_followup_mids,glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))
#modelchange25 <- summary(pool(change25), conf.int = T)
#exp(modelchange25$estimate)
#exp(modelchange25$`2.5 %`)
#exp(modelchange25$`97.5 %`)
##summary(pool(with(bmi_followup_mids,glm((bmi.fu>=25) ~ (basebmi25+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))))

##change 30
#change30 <- with(bmi_followup_mids,glm((bmi.fu>=30) ~ (basebmi30+selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))
#modelchange30 <- summary(pool(change30), conf.int = T)
#exp(modelchange30$estimate)
#exp(modelchange30$`2.5 %`)
#exp(modelchange30$`97.5 %`)


## ----- ##
#Final change analyses (per Naja's requests)
## ----- ##

## change from low to high group for the subjects at risk
#bmi_followup_mids_25risk <- as.mids(subset(bmi_followup,bmi.base<25 & !is.na(sample_weights)),.imp="imputation")
#summary(pool(with(bmi_followup_mids_25risk,glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)
#bmi_followup_mids_30risk <- as.mids(subset(bmi_followup,bmi.base<30 & !is.na(sample_weights)),.imp="imputation")
#summary(pool(with(bmi_followup_mids,glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights,family=binomial))), conf.int = T)

## from below 25 to above 25
bmi_followup_mids_25risk <- as.mids(subset(bmi_followup,bmi.base<25),.imp="imputation")
summary(pool(with(bmi_followup_mids_25risk,glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y,family=binomial))), conf.int = T)

NewBmi25 <- with(bmi_followup_mids_25risk,glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y,family=binomial))
ModelNewBmi25 <- summary(pool(NewBmi25), conf.int = T)
exp(ModelNewBmi25$estimate)
exp(ModelNewBmi25$`2.5 %`)
exp(ModelNewBmi25$`97.5 %`)

## from below 30 to above 30
#bmi_followup_mids_30risk <- as.mids(subset(bmi_followup_30risk,bmi.base<30),.imp="imputation")
bmi_followup_mids_30risk <- as.mids(subset(bmi_followup,bmi.base<30),.imp="imputation")
#summary(pool(with(bmi_followup_mids,glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y,family=binomial))), conf.int = T)
newBmi30 <- with(bmi_followup_mids,glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y), weights=sample_weights.y,family=binomial))
modelnewBmi30 <- summary(pool(newBmi30), conf.int = T)
exp(modelnewBmi30$estimate)
exp(modelnewBmi30$`2.5 %`)
exp(modelnewBmi30$`97.5 %`)


## And then the long format for the numeric change:

#long_data <- data.frame("bmi"=c(bmi_followup$bmi.base,bmi_followup$bmi.fu),"userid"=bmi_followup$userid,"sample_weights"=bmi_followup$sample_weights.y,"gender"=bmi_followup$gender.y,"age"=bmi_followup$age.y,
#                        education=bmi_followup$education.y,occupation=bmi_followup$occupation.y,selfScoreCat = bmi_followup$selfScoreCat.y,"time"=c(rep(0,length(bmi_followup$bmi.base)),rep(1,length(bmi_followup$bmi.fu))),"imputation"=bmi_followup$imputation)

#long_data_mids <- as.mids(long_data[!is.na(long_data$sample_weights),],.imp="imputation")

m <- gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time,
            nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==10)),family=BCCG,method=RS(100),robust=T) #(selfScoreCat+age+gender+education+occupation)*time
logL <- gen.likelihood(m)
conf.res <- confint(m,what=c("mu")) #Wald type

prof.dev(m,which="mu",min=1,max=25) #Doesn't work
prof.term(model=m,criterion="GD",min=-35,max=25,step=1)$CI #Doesn't seem to work

#Manual comparison of likelihood based confidence intervals and wald type intervals


coef <- sds <- matrix(nrow=length(c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)),ncol=20)
for (i in 1:20){
  m <- gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time,
              nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==i)),family=BCCG,method=RS(100),robust=T) #(selfScoreCat+age+gender+education+occupation)*time
  
  coef[,i] <- c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)
  sds[,i] <- sqrt(diag(vcov(m)))
}

ests <- estimate.pooler(coef,sds)[,1]

logL1 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==1)),family=BCCG,method=RS(100)))
logL2 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==2)),family=BCCG,method=RS(100)))
logL3 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==3)),family=BCCG,method=RS(100)))
logL4 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==4)),family=BCCG,method=RS(100)))
logL5 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==5)),family=BCCG,method=RS(100)))
logL6 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==6)),family=BCCG,method=RS(100)))
logL7 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==7)),family=BCCG,method=RS(100)))
logL8 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==8)),family=BCCG,method=RS(100)))
logL9 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==9)),family=BCCG,method=RS(100)))
logL10 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==10)),family=BCCG,method=RS(100)))
logL11 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==11)),family=BCCG,method=RS(100)))
logL12 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==12)),family=BCCG,method=RS(100)))
logL13 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==13)),family=BCCG,method=RS(100)))
logL14 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==14)),family=BCCG,method=RS(100)))
logL15 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==15)),family=BCCG,method=RS(100)))
logL16 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==16)),family=BCCG,method=RS(100)))
logL17 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==17)),family=BCCG,method=RS(100)))
logL18 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==18)),family=BCCG,method=RS(100)))
logL19 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==19)),family=BCCG,method=RS(100)))
logL20 <- gen.likelihood(gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ time, nu.formula = ~ time,weights=sample_weights, data=na.omit(subset(long_data,imputation==20)),family=BCCG,method=RS(100)))


logLhat <- (logL1(ests)+logL2(ests)+logL3(ests)+logL4(ests)+logL5(ests)+
              logL6(ests)+logL7(ests)+logL8(ests)+logL9(ests)+logL10(ests)+
              logL11(ests)+logL12(ests)+logL13(ests)+logL14(ests)+logL15(ests)+
              logL16(ests)+logL17(ests)+logL18(ests)+logL19(ests)+logL20(ests))/20

#Would ideally put the pooled estimates in here, and use them for generating profile likelihood intervals in the AVERAGE likelihood function across the imputations... Hence we would need to also generate that function by generating each of the likelihoods and making the average of the 20 likelihoods evaluated in given parameters.
#hatmucoefs <- c() #This should then contain the pooled estimates to put into the 20 likelihoods.

change_seq <- change_seq1 <- seq(from=-2.5,to=2.5,by=0.01)
change_seq2 <- seq(from=-0.1,to=0.1,by=0.0001)
out_seq <- numeric(0)
out_ints <- matrix(0,nrow=length(m$mu.coefficients),ncol=2)


for (k in 1:4){
for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
#if (k==1){
  out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                     logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])))))
  #}
#if (k>1){
#out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
#}
}
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}


out_seq <- numeric(0)
k=5
for (i in 1:(length(change_seq2))){ #*(k==1)+length(change_seq2)*(k>1)
  out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL2(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL3(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL4(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL5(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL6(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL7(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL8(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL9(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL10(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL11(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL12(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL13(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL14(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL15(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL16(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL17(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL18(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL19(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL20(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])))))
  #}
  #if (k>1){
  #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq22[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
  #}
}
out_ints[k,] <- c(ests[k]+min(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]))
plot(m$mu.coefficients[k]+change_seq2,out_seq)


out_seq <- numeric(0)
for (k in 6:(length(m$mu.coefficients))){
  for (i in 1:(length(change_seq2))){ #*(k==1)+length(change_seq2)*(k>1)
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL2(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL3(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL4(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL5(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL6(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL7(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL8(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL9(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL10(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL11(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL12(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL13(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL14(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL15(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL16(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL17(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL18(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL19(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                       logL20(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])))))
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq22[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq2,out_seq)
}


#Comparison in intervals
cbind(out_ints,conf.res,coef(m))

#Conclusion in this case: The likelihood based intervals are very narrow... Is something wrong with the implementation?
#If the implementation can be trusted the Wald type confidence intervals are very conservative.

#The plots are quite nicely quadratic in their shapes. Hence the Wald approximation is good and there must be an issue with the implementation above.

#Generally:
#Wald intervals with Robust=T are better for misspecified models.
#Profile likelihood intervals are better for models that are close to correct (mostly so for smaller sample sizes).

#Here the important estimates are the time x category interaction term estimates - but the pooler only provides a single imputation fit...
summary(pool(with(long_data_mids,gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,
                                        nu.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,weights=sample_weights,family=BCCG,method=RS(100),robust=T))))
                                        
## 95%CI
confint(pool(long_data_mids,gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,
                                   nu.formula = ~ (selfScoreCat+age+gender+education+occupation)*time,weights=sample_weights,family=BCCG,method=RS(100))))


confint(pool(with(long_data_mids,gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ 1,
                                        nu.formula = ~ 1,weights=sample_weights,family=BCCG,method=RS(100)))),type="likelihood") #It doesn't want to give me what I want...

logLpool <- gen.likelihood(pool(with(long_data_mids,gamlss(bmi~(selfScoreCat+age+gender+education+occupation)*time, sigma.formula = ~ 1,
                                                           nu.formula = ~ 1,weights=sample_weights,family=BCCG,method=RS(100)))))

#Bonus
summary(pool(with(long_data_mids,glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*time,family=binomial, weights=sample_weights))))
summary(pool(with(long_data_mids,glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*time,family=binomial, weights=sample_weights))))



# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#####Tracking data for the followup CSS sample

## Hvad er de forskellige clusters? cluster 1 = night-time user?, cluster 2 =morning user?, cluster = non-user?, cluster 4 = evening user?

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

summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))))

m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
            nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100))
#robust=TRUE? Doesn't change much when the fit is this good.

#m.pool <- pool(with(CSS_track_mids,gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
#                                          nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights,family=BCCG,method=RS(100))))
summary(pool(with(CSS_track_mids,gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
                           nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights,family=BCCG,method=RS(100)))))

#Checking Wald validity

logL <- gen.likelihood(m)

logLhat <- logL()


coef <- sds <- matrix(nrow=length(c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)),ncol=20)
for (i in 1:20){
  m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
              nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==i)),family=BCCG,method=RS(100))
  
  coef[,i] <- c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)
  sds[,i] <- sqrt(diag(vcov(m)))
}

ests <- estimate.pooler(coef,sds)[,1]

logL1 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==1)),family=BCCG,method=RS(100)))
logL2 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==2)),family=BCCG,method=RS(100)))
logL3 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==3)),family=BCCG,method=RS(100)))
logL4 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==4)),family=BCCG,method=RS(100)))
logL5 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==5)),family=BCCG,method=RS(100)))
logL6 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==6)),family=BCCG,method=RS(100)))
logL7 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==7)),family=BCCG,method=RS(100)))
logL8 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==8)),family=BCCG,method=RS(100)))
logL9 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==9)),family=BCCG,method=RS(100)))
logL10 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100)))
logL11 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==11)),family=BCCG,method=RS(100)))
logL12 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==12)),family=BCCG,method=RS(100)))
logL13 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==13)),family=BCCG,method=RS(100)))
logL14 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==14)),family=BCCG,method=RS(100)))
logL15 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==15)),family=BCCG,method=RS(100)))
logL16 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==16)),family=BCCG,method=RS(100)))
logL17 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==17)),family=BCCG,method=RS(100)))
logL18 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==18)),family=BCCG,method=RS(100)))
logL19 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==19)),family=BCCG,method=RS(100)))
logL20 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==20)),family=BCCG,method=RS(100)))


logLhat <- (logL1(ests)+logL2(ests)+logL3(ests)+logL4(ests)+logL5(ests)+
              logL6(ests)+logL7(ests)+logL8(ests)+logL9(ests)+logL10(ests)+
              logL11(ests)+logL12(ests)+logL13(ests)+logL14(ests)+logL15(ests)+
              logL16(ests)+logL17(ests)+logL18(ests)+logL19(ests)+logL20(ests))/20

change_seq <- change_seq1 <- seq(from=-2.5,to=2.5,by=0.01)
change_seq2 <- seq(from=-0.1,to=0.1,by=0.0001)
out_seq <- numeric(0)

out_ints <- matrix(0,nrow=length(m$mu.coefficients),ncol=2)

for (k in 1:7){
  for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
    #if (k==1){
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])))))
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}

out_seq <- numeric(0)
k=8
for (i in 1:(length(change_seq2))){ #*(k==1)+length(change_seq2)*(k>1)
  #if (k==1){
  out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL2(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL3(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL4(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL5(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL6(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL7(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL8(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL9(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL10(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL11(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL12(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL13(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL14(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL15(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL16(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL17(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL18(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL19(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL20(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])))))
  #}
  #if (k>1){
  #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq22[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
  #}
}
out_ints[k,] <- c(ests[k]+min(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]))
plot(m$mu.coefficients[k]+change_seq2,out_seq)

out_seq <- numeric(0)
for (k in 9:(length(m$mu.coefficients))){
  for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
    #if (k==1){
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])))))
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}

#The resulting plots should be quadratically shaped around the MLE, if we would like to use the Wald approximation.


# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#Tracking data: Population sample (random sample) - same analysis

hist(pop_track$bmi,breaks=50,xlim=c(0,50))

ggplot(pop_track, aes(x = factor(selfScoreCat))) +
  geom_bar()
ggplot(pop_track, aes(x = factor(cluster))) +
  geom_bar()

#analyses

## regression analysis of clusters of night-time smartphone use and overweight/obesity #justeres for selfScoreCat?? ## hvad er de forskellige clusters?? ## fortolkning?? ## inkluderer imp_nr=0?
## bmi kontinuert 

m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
            nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100))

## no adjustment for risk profiles
summary(pool(with(pop_track_mids,gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation),
                                        nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation),weights=sample_weights,family=BCCG,method=RS(100)))))

## no adjustment for clusters
summary(pool(with(pop_track_mids,gamlss(bmi~(selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (selfScoreCat+age+gender+education+occupation),
                                        nu.formula = ~ (selfScoreCat+age+gender+education+occupation),weights=sample_weights,family=BCCG,method=RS(100)))))


## mutually adjusting for clusters/risk profiles
summary(pool(with(pop_track_mids,gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
                                        nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights,family=BCCG,method=RS(100)))))

#Checking validity of Wald yet again

logL <- gen.likelihood(m)



coef <- sds <- matrix(nrow=length(c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)),ncol=20)
for (i in 1:20){
  m <- gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),
              nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==i)),family=BCCG,method=RS(100))
  
  coef[,i] <- c(m$mu.coefficients,m$sigma.coefficients,m$nu.coefficients)
  sds[,i] <- sqrt(diag(vcov(m)))
}

ests <- estimate.pooler(coef,sds)[,1]

logL1 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==1)),family=BCCG,method=RS(100)))
logL2 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==2)),family=BCCG,method=RS(100)))
logL3 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==3)),family=BCCG,method=RS(100)))
logL4 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==4)),family=BCCG,method=RS(100)))
logL5 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==5)),family=BCCG,method=RS(100)))
logL6 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==6)),family=BCCG,method=RS(100)))
logL7 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==7)),family=BCCG,method=RS(100)))
logL8 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==8)),family=BCCG,method=RS(100)))
logL9 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==9)),family=BCCG,method=RS(100)))
logL10 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100)))
logL11 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==11)),family=BCCG,method=RS(100)))
logL12 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==12)),family=BCCG,method=RS(100)))
logL13 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==13)),family=BCCG,method=RS(100)))
logL14 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==14)),family=BCCG,method=RS(100)))
logL15 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==15)),family=BCCG,method=RS(100)))
logL16 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==16)),family=BCCG,method=RS(100)))
logL17 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==17)),family=BCCG,method=RS(100)))
logL18 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==18)),family=BCCG,method=RS(100)))
logL19 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==19)),family=BCCG,method=RS(100)))
logL20 <- gen.likelihood(gamlss(bmi~(cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),nu.formula = ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation),weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==20)),family=BCCG,method=RS(100)))


logLhat <- (logL1(ests)+logL2(ests)+logL3(ests)+logL4(ests)+logL5(ests)+
              logL6(ests)+logL7(ests)+logL8(ests)+logL9(ests)+logL10(ests)+
              logL11(ests)+logL12(ests)+logL13(ests)+logL14(ests)+logL15(ests)+
              logL16(ests)+logL17(ests)+logL18(ests)+logL19(ests)+logL20(ests))/20

logLhat <- logL()
change_seq <- change_seq1 <- seq(from=-2.5,to=2.5,by=0.01)
change_seq2 <- seq(from=-0.1,to=0.1,by=0.0001)
out_seq <- numeric(0)

out_ints <- matrix(0,nrow=length(m$mu.coefficients),ncol=2)

for (k in 1:7){
  for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
    #if (k==1){
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])))))
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}

out_seq <- numeric(0)
k=8
for (i in 1:(length(change_seq2))){ #*(k==1)+length(change_seq2)*(k>1)
  #if (k==1){
  out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL2(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL3(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL4(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL5(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL6(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL7(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL8(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL9(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL10(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL11(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL12(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL13(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL14(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL15(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL16(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL17(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL18(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL19(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])),
                                     logL20(c(ests[0:(k-1)],ests[k]+change_seq2[i],ests[(k+1):length(ests)])))))
  #}
  #if (k>1){
  #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq22[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
  #}
}
out_ints[k,] <- c(ests[k]+min(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq2[which(out_seq<=qchisq(p=0.95,df=1))]))
plot(m$mu.coefficients[k]+change_seq2,out_seq)

out_seq <- numeric(0)
for (k in 9:(length(m$mu.coefficients))){
  for (i in 1:(length(change_seq))){ #*(k==1)+length(change_seq2)*(k>1)
    #if (k==1){
    out_seq[i] <- -2*(logLhat - mean(c(logL1(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL2(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL3(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL4(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL5(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL6(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL7(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL8(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL9(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL10(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL11(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL12(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL13(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL14(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL15(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL16(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL17(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL18(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL19(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])),
                                       logL20(c(ests[0:(k-1)],ests[k]+change_seq[i],ests[(k+1):length(ests)])))))
    #}
    #if (k>1){
    #out_seq[i] <- abs(logLhat - logL(c(m$mu.coefficients[0:(k-1)],m$mu.coefficients[k]+change_seq2[i],m$mu.coefficients[(k+1):length(m$mu.coefficients)],m$sigma.coefficients,m$nu.coefficients)))
    #}
  }
  out_ints[k,] <- c(ests[k]+min(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]),ests[k]+max(change_seq[which(out_seq<=qchisq(p=0.95,df=1))]))
  plot(m$mu.coefficients[k]+change_seq,out_seq)
}

#The resulting plots should be quadratically shaped around the MLE, if we would like to use the Wald approximation. 



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
#summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
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


####### BMI and self-reported risk profiles

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

## test for trend (selfScoreCat numeric)
Random30Risk <- with(pop_track_mids,glm((bmi>=30) ~ ((as.numeric(selfScoreCat))+age+gender+education+occupation), weights=sample_weights,family=binomial))


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

## BMI clinical
table(clinical_sample$bmi.clinical)
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable( ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

## descriptives of clinical sample
publish(univariateTable(selfScoreCat ~ age.x,data=clinical_sample, column.percent=TRUE))

publish(univariateTable(selfScoreCat ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

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
clinical_sample$sbpCat[clinical_sample$sbp>=140] <- 1#"High"
clinical_sample$sbpCat[clinical_sample$sbp<140] <- 0#"Normal"
publish(univariateTable( ~ sbpCat,data=clinical_sample, column.percent=TRUE))

# diastolic blood pressure
clinical_sample$dbp<-rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ dbp,data=clinical_sample, column.percent=TRUE))

clinical_sample$dbpCat[clinical_sample$dbp>=90] <- 1#"High"
clinical_sample$dbpCat[clinical_sample$dbp<90] <- 0#"Normal"
table(clinical_sample$dbpCat, useNA="always")

## hip waist ratio

clinical_sample$ratiowaisthip <- as.numeric(clinical_sample$ratiowaisthip)
publish(univariateTable(selfScoreCat ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))

clinical_sample$ratiowaisthipCat[clinical_sample$ratiowaisthip>=0.85] <- 1#"High"
clinical_sample$ratiowaisthipCat[clinical_sample$ratiowaisthip<0.85] <- 0#"Normal"
table(clinical_sample$ratiowaisthipCat)

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
hist(rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T),breaks=20)
hist(rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T),breaks=20)

## HDL
clinical_sample$hdl <- as.numeric(clinical_sample$hdl)
publish(univariateTable(selfScoreCat ~ hdl,data=clinical_sample, column.percent=TRUE))

## categorize
clinical_sample$hdlCat[clinical_sample$hdl<=1.2] <- 1#"Bad"
clinical_sample$hdlCat[clinical_sample$hdl>1.2] <- 0#"Good"
table(clinical_sample$hdlCat)

publish(univariateTable(selfScoreCat ~ hdlCat,data=clinical_sample, column.percent=TRUE))

## LDL
clinical_sample$ldl <- as.numeric(clinical_sample$ldl)
publish(univariateTable( selfScoreCat~ ldl,data=clinical_sample, column.percent=TRUE))

## categorize LDL
clinical_sample$ldlCat[clinical_sample$ldl>=4.3] <- 1#"Bad"
clinical_sample$ldlCat[clinical_sample$ldl<4.3] <- 0#"Good"
publish(univariateTable(selfScoreCat ~ ldlCat,data=clinical_sample, column.percent=TRUE))

## triglycerides

clinical_sample$triglycerids <- as.numeric(clinical_sample$triglycerids)
publish(univariateTable(selfScoreCat ~ triglycerids,data=clinical_sample, column.percent=TRUE))

clinical_sample$triglyceridsCat[clinical_sample$triglycerids>=1.2] <- 1#"Bad"
clinical_sample$triglyceridsCat[clinical_sample$triglycerids<1.2] <- 0#"Good"
publish(univariateTable(selfScoreCat ~ triglyceridsCat,data=clinical_sample, column.percent=TRUE))


## hba1c
clinical_sample$hba1c <- as.numeric(clinical_sample$hba1c)
publish(univariateTable(selfScoreCat ~ hba1c,data=clinical_sample, column.percent=TRUE))

## categorize hba1c
clinical_sample$hba1cCat[clinical_sample$hba1c>=44] <- 1#"Bad"
clinical_sample$hba1cCat[clinical_sample$hba1c<44] <- 0#"Good"
publish(univariateTable(selfScoreCat ~ hba1cCat,data=clinical_sample, column.percent=TRUE))


## total cholesterol
clinical_sample$t_cholesterol <- as.numeric(clinical_sample$t_cholesterol)
publish(univariateTable(selfScoreCat ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))

clinical_sample$t_cholesterolCat[clinical_sample$t_cholesterol>=5] <- 1#"Bad"
clinical_sample$t_cholesterolCat[clinical_sample$t_cholesterol<5] <- 0#"Good"
publish(univariateTable(selfScoreCat ~ t_cholesterolCat,data=clinical_sample, column.percent=TRUE))


#Models - multiple testing issue if we are going to 'pick and choose' which responses we would like to look at.
table(clinical_sample$age.x)

#Transforming to mids for modelling and inference

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")

#hdl
#hdl_sum<-cbind(summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
lines(seq(from=min(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),length.out=100),dnorm(x=seq(from=min(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),length.out=100),mean=mean(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))),sd=sd(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))))
plot(residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))


#ldl 
#ldl_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res=residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))
plot(fitted(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

cbind(confint(glm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#vldl
#vldl_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),
      cbind(coef(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))-1.96*summary(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2],
            coef(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))+1.96*summary(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2]))


#t_cholesterol
#t_cholesterol_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#                summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(lm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(lm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#triglycerids
#tri_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, glm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(glm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),
      cbind(coef(glm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))-1.96*summary(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2],
            coef(glm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))+1.96*summary(glm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2]))

#hba1c
#hba1c_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res <- residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

#glucose
#glu_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res <- residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))
plot(fitted(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

cbind(confint(glm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#ratiowaisthip
#wh_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation))))$std.error[2:4])
hist(residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#sbp
#sbp_sum<-cbind(summary(pool(with(data=clinical_mids, lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(sbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(sbp) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(sbp) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#dbp
#dbp_sum<-cbind(summary(pool(with(data=clinical_mids, lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(dbp) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(dbp) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#Generally the residuals look reasonably centered, with a few positive outliers. The residual distributions on the first imputation actually look reasonably normal, save for the few (extreme) outliers.

#Seems that these models are appropriate, and that normal approximations of confidence interval will be reasonable too. We can however also just use the profile likelihood CI's.

#Wald confidence intervals

dbp_int <- summary(pool(with(data=clinical_mids, lm(dbp ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_int <- summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_int <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_ints <- data.frame(rbind(dbp_int[2,],sbp_int[2,],hba1c_int[2,],hdl_int[2,],ldl_int[2,],vldl_int[2,],t_chol_int [2,],tri_int[2,],glu_int[2,],wh_int[2,]),
                      rbind(dbp_int[3,],sbp_int[3,],hba1c_int[3,],hdl_int[3,],ldl_int[3,],vldl_int[3,],t_chol_int [3,],tri_int[3,],glu_int[3,],wh_int[3,]),
                      rbind(dbp_int[4,],sbp_int[4,],hba1c_int[4,],hdl_int[4,],ldl_int[4,],vldl_int[4,],t_chol_int [4,],tri_int[4,],glu_int[4,],wh_int[4,]))

colnames(df_ints) <- c("type.1.estimate","type.1.lower","type.1.upper","type.1.pvalue","type.2.estimate","type.2.lower","type.2.upper","type.2.pvalue","type.4.estimate","type.4.lower","type.4.upper","type.4.pvalue")
rownames(df_ints) <- c("dbp","sbp","hba1c","hdl","ldl","vldl","total cholesterol","triglycerides","glucose","w/h ratio")
df_ints

# ----- #

#Now looking instead at classifications - where there are quite few poeple with values above the thresholds...

#hdl_class_sum<-cbind(summary(pool(with(clinical_mids,glm(as.numeric(hdlCat) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,family=binomial))))$estimate[2:4],
#               summary(pool(with(clinical_mids,glm(as.numeric(hdlCat) ~ cluster1prob+cluster2prob+cluster4prob+age.x+education+occupation,na.action=na.omit,family=binomial))))$std.error[2:4])

#ldl_class_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(ldlCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$estimate[2:4],
#               summary(pool(with(data=clinical_mids, glm(as.numeric(ldlCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$std.error[2:4])

#t_chol_class_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterolCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$estimate[2:4],
#                summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterolCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$std.error[2:4])

#tri_class_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(triglyceridsCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$estimate[2:4],
#               summary(pool(with(data=clinical_mids, glm(as.numeric(triglyceridsCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$std.error[2:4])

#hba1c_class_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(hba1cCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$estimate[2:4],
#                 summary(pool(with(data=clinical_mids, glm(as.numeric(hba1cCat) ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))))$std.error[2:4])

hdl_class_int <- summary(pool(with(data=clinical_mids, glm(hdlCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_class_sum[,1]-1.96*hdl_class_sum[,2],hdl_class_sum[,1]+1.96*hdl_class_sum[,2])

ldl_class_int <- summary(pool(with(data=clinical_mids, glm(ldlCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_class_sum[,1]-1.96*ldl_class_sum[,2],ldl_class_sum[,1]+1.96*ldl_class_sum[,2])

t_chol_class_int <- summary(pool(with(data=clinical_mids, glm(t_cholesterolCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_class_sum[,1]-1.96*t_chol_class_sum[,2],t_chol_class_sum[,1]+1.96*t_chol_class_sum[,2])

tri_class_int <- summary(pool(with(data=clinical_mids, glm(triglyceridsCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_class_sum[,1]-1.96*tri_class_sum[,2],tri_class_sum[,1]+1.96*tri_class_sum[,2])

hba1c_class_int <- summary(pool(with(data=clinical_mids, glm(hba1cCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_class_sum[,1]-1.96*hba1c_class_sum[,2],hba1c_class_sum[,1]+1.96*hba1c_class_sum[,2])

dbp_class_int <- summary(pool(with(data=clinical_mids, glm(dbpCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]

sbp_class_int <- summary(pool(with(data=clinical_mids, glm(sbpCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]

rwh_class_int <- summary(pool(with(data=clinical_mids, glm(ratiowaisthipCat ~ cluster1prob+cluster2prob+cluster4prob+age+gender+education+occupation,na.action=na.omit,family=binomial))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]


df_class_ints <- data.frame(rbind(hdl_class_int[2,],ldl_class_int[2,],t_chol_class_int[2,],tri_class_int[2,],hba1c_int[2,],dbp_class_int[2,],sbp_class_int[2,],rwh_class_int[2,]),
                            rbind(hdl_class_int[3,],ldl_class_int[3,],t_chol_class_int[3,],tri_class_int[3,],hba1c_int[3,],dbp_class_int[3,],sbp_class_int[3,],rwh_class_int[3,]),
                            rbind(hdl_class_int[4,],ldl_class_int[4,],t_chol_class_int[4,],tri_class_int[4,],hba1c_int[4,],dbp_class_int[4,],sbp_class_int[4,],rwh_class_int[4,]))

colnames(df_class_ints) <- c("type.1.estimate","type.1.lower","type.1.upper","type.1.pvalue","type.2.estimate","type.2.lower","type.2.upper","type.2.pvalue","type.4.estimate","type.4.lower","type.4.upper","type.4.pvalue")
rownames(df_class_ints) <- c("hdl","ldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio")
df_class_ints
