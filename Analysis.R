###### Analysis of associations between metabolic biomarkers and smartphone activity ######

#IMPORTANT: DATA WILL CONTAIN WEIGHTS AND IMP NUMBERS FROM THE BEGINNING. SO THERE WILL BE NO NEED TO COLLECT THESE. HENCE, DELETE THE LINES THAT ADD THESE TO THE DATA.


#This script contains analyses of the multiply imputed data sets from the SmartSleep project.
#The multiple imputation results are combined using Rubin's rules.

library(dplyr)
library(matrixStats)
library(mice)
library(MASS)

#Reading in the data
#setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
#load("H:/SmartSleep backup IT issues/imputation/myImputationCSS-res-00001.RData")
#CSS0 <- imp_CSS$data
#CSS0$impnr=0
#CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
#CSS <- rbind(CSS0,CSS)

subject_tracking_clusters <- read.csv2("H:/SmartSleep backup IT issues/subject_tracking_clusters.csv") #Change to S-drive location

#base_weights <- read.csv2("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/SmartSleepExpWeighted.csv")
#load("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/SmartSleep Experiment/full_imp_base.RData")
#base_data <- full_imp_Base
#base_data <- rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv"),imputation=imp_nr)
#base_data$zipCode<-as.numeric(base_data$zipCode)

base_data <- rename(read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv"),imputation=imp_nr)

#Baseline data with self-reports

base_data$selfScore <- (base_data$mobileUseBeforeSleep=="5-7 times per week")*4+(base_data$mobileUseBeforeSleep=="2-4 times per week")*3+(base_data$mobileUseBeforeSleep=="Once a week")*3+(base_data$mobileUseBeforeSleep=="Every month or less")*2+(base_data$mobileUseBeforeSleep=="Never")*1+
  (base_data$mobileUseNight=="Every night or almost every night")*4+(base_data$mobileUseNight=="A few times a week")*3+(base_data$mobileUseNight=="A few times a month or less")*2+(base_data$mobileUseNight=="Never")*1+
  (base_data$mobileCheck==">20 times an hour")*4+(base_data$mobileCheck=="11-20 times an hour")*4+(base_data$mobileCheck=="5-10 times an hour")*3+(base_data$mobileCheck=="1-4 times an hour")*2+(base_data$mobileCheck=="Every 2nd hour")*2+(base_data$mobileCheck=="Several times a day")*1+(base_data$mobileCheck=="Once a day or less")*1+
  (base_data$pmpuScale<=14)*1+(base_data$pmpuScale>14 & base_data$pmpuScale<=16.76)*2+(base_data$pmpuScale>16.76 & base_data$pmpuScale<=19)*3+(base_data$pmpuScale>19)*4
base_data$selfScoreCat <- "1"
base_data$selfScoreCat[base_data$selfScore>=8]="2"
base_data$selfScoreCat[base_data$selfScore>=9.517]="3"
base_data$selfScoreCat[base_data$selfScore>=12]="4"


#Followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")
CSS_track <- inner_join(CSS,subject_tracking_clusters,by="userid")

CSS_track$selfScore <- (CSS_track$mobileUseBeforeSleep=="5-7 times per week")*4+(CSS_track$mobileUseBeforeSleep=="2-4 times per week")*3+(CSS_track$mobileUseBeforeSleep=="Once a week")*3+(CSS_track$mobileUseBeforeSleep=="Every month or less")*2+(CSS_track$mobileUseBeforeSleep=="Never")*1+
  (CSS_track$mobileUseNight=="Every night or almost every night")*4+(CSS_track$mobileUseNight=="A few times a week")*3+(CSS_track$mobileUseNight=="A few times a month or less")*2+(CSS_track$mobileUseNight=="Never")*1+
  (CSS_track$mobileCheck==">20 times an hour")*4+(CSS_track$mobileCheck=="11-20 times an hour")*4+(CSS_track$mobileCheck=="5-10 times an hour")*3+(CSS_track$mobileCheck=="1-4 times an hour")*2+(CSS_track$mobileCheck=="Every 2nd hour")*2+(CSS_track$mobileCheck=="Several times a day")*1+(CSS_track$mobileCheck=="Once a day or less")*1+
  (CSS_track$pmpuScale<=14)*1+(CSS_track$pmpuScale>14 & CSS_track$pmpuScale<=16.76)*2+(CSS_track$pmpuScale>16.76 & CSS_track$pmpuScale<=19)*3+(CSS_track$pmpuScale>19)*4
CSS_track$selfScoreCat <- "1"
CSS_track$selfScoreCat[CSS_track$selfScore>=8]="2"
CSS_track$selfScoreCat[CSS_track$selfScore>=9.517]="3"
CSS_track$selfScoreCat[CSS_track$selfScore>=12]="4"


#Looking at general patterns

#Cakculated bmi from imputed height and weight vs directly imputed bmi - which one to use? At the moment we use the directly imputed bmi.
bmi_ind <- which(abs(CSS$bmi-CSS$weight/((CSS$height/100)^2))>1)
cor(CSS$bmi[bmi_ind],(CSS$weight/((CSS$height/100)^2))[bmi_ind])
plot(CSS$bmi[bmi_ind],(CSS$weight/((CSS$height/100)^2))[bmi_ind])
CSS$bmi[CSS$bmi==0]<-NA

hist(CSS$bmi,breaks=40)
hist(log(CSS$bmi),breaks=40)
hist(log(log(CSS$bmi)),breaks=40)
hist(sqrt(CSS$bmi),breaks=40)
hist(CSS$age,breaks=40)
hist(CSS$height,breaks=40)
hist(CSS$weight,breaks=40)


bmi_ind <- which(abs(base_data$bmi-base_data$weight/((base_data$height/100)^2))>1)
cor(base_data$bmi[bmi_ind],(base_data$weight/((base_data$height/100)^2))[bmi_ind])
plot(base_data$bmi[bmi_ind],(base_data$weight/((base_data$height/100)^2))[bmi_ind])
base_data$bmi[base_data$bmi==0]<-NA

hist(base_data$bmi,breaks=40)
hist(log(base_data$bmi),breaks=40)
hist(log(log(base_data$bmi)),breaks=40)
hist(sqrt(base_data$bmi),breaks=40)
hist(base_data$age,breaks=40)
hist(base_data$height,breaks=40)
hist(base_data$weight,breaks=40)



#Regression analyses


#Base BMI models

summary(glm((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0),family=binomial))
summary(glm((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0),family=binomial))
summary(lm(log(log(bmi))~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0)))

plot(fitted(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))))
hist(residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),xlim=c(-20,20),breaks=200)
hist(simulate(lm(bmi~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),breaks=40) #The bell-shape is not that well suited
hist(simulate(lm(log(bmi)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),breaks=40) #The bell-shape is not that well suited here either

#box_cox_transformation - bmi values are heavily right-skewed so we transform them closer toward normality
bc <- boxcox(bmi ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))
(lambda <- bc$x[which.max(bc$y)])
new_model <- lm(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))

hist(((base_data$bmi^lambda-1)/lambda),breaks=40,xlim=c(0.91,0.96))
hist(simulate(new_model),breaks=10,xlim=c(0.91,0.96))

summary(new_model)



#Method 1: Using the mice package with mids objects
base_data_mids <- as.mids(base_data,.imp="imputation")
mod30<-(glm.mids((bmi>=30)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data_mids,family=binomial))
mod25<-(glm.mids((bmi>=25)~(selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data_mids,family=binomial))
modnum<-(lm.mids(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data_mids))

summary(pool(mod30))
summary(pool(mod25))
summary(pool(modnum))

#Method 2: Estimates adjusted with Rubin's rules (raw computations)

sd25<-sd30<-sdnum<-
  matrix(nrow=length(summary(glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data))$coefficients[,2]),ncol=20)

est25<-est30<-estnum<-
  matrix(nrow=length(summary(glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=base_data))$coefficients[,1]),ncol=20)


for (i in 1:20){
  sd25[,i]<-summary(glm((bmi>=25) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation==i),family=binomial))$coefficients[,2]
  sd30[,i]<-summary(glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation==i),family=binomial))$coefficients[,2]
  sdnum[,i]<-summary(lm(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation==i)))$coefficients[,2]
  
  est25[,i]<-summary(glm((bmi>=25) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation==i),family=binomial))$coefficients[,1]
  est30[,i]<-summary(glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation==i),family=binomial))$coefficients[,1]
  estnum[,i]<-summary(lm(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation==i)))$coefficients[,1]
}

sd25Value<-sqrt(rowMeans(sd25^2)+rowSums((est25-rowMeans(est25))^2)/19) # sqrt(mean estimator variance + estimated imputation estimator variance): SE=sqrt(V_total).
sd30Value<-sqrt(rowMeans(sd30^2)+rowSums((est30-rowMeans(est30))^2)/19)
sdnumValue<-sqrt(rowMeans(sdnum^2)+rowSums((estnum-rowMeans(estnum))^2)/19)


#Summary tables:
bmi25_indicator_table <- data.frame("ests"<-coef(glm((bmi>=25) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0),family=binomial)),
                                       "sds"<- sd25Value)
bmi30_indicator_table <- data.frame("ests"<-coef(glm((bmi>=30) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0),family=binomial)),
                                         "sds"<- sd30Value)
bmi_num_table <- data.frame("ests"<-coef(lm(((bmi^lambda-1)/lambda) ~ (selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(base_data,imputation!=0))),
                                           "sds"<- sdnumValue)


#BMI followup difference - match with emailAddress

bmi_followup <- inner_join(rename(CSS_track,imputation=impnr),base_data,by=c("emailAdress","imputation"))
bmi_followup$difference <- bmi_followup$bmi.x-bmi_followup$bmi.y

hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,1500))
hist(simulate(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)),xlim=c(-10,10)) #Too wide?
summary(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation!=0)))

plot(fitted(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)),residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)))
hist(residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup)),breaks=200,xlim=c(-10,10))

#The differences are not skewed.

#Method 1: Using the mids object
bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

summay(lm.mids(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids))


#Method 2: Rubin's rules for variance in raw computations

bmi_followup_ests <- matrix(nrow=length(coef(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup))),ncol=20)
bmi_followup_sds <- matrix(nrow=length(coef(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup))),ncol=20)

for (i in 1:20){
bmi_followup_ests[,i] <- summary(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i)))$coefficients[,1]
bmi_followup_sds[,i] <- summary(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i)))$coefficients[,2]
}

#Summary table
bmi_followup_results <- data.frame("est"=rowMeans(bmi_followup_ests),
                                   "sd"=sqrt(rowMeans(bmi_followup_sds^2)+rowSums((bmi_followup_ests-rowMeans(bmi_followup_ests))^2)/19))

bmi_followup_results$sd-summary(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==1)))$coefficients[,2]


#Differences for indicators also (3 models per threshold: change, change from low to high group, and change from high to low group)

bmi_followup$bmi25change <- as.numeric((bmi_followup$bmi.x>=25)!=(bmi_followup$bmi.y>=25))
bmi_followup$bmi25changeUp <- as.numeric((bmi_followup$bmi.x>=25)>(bmi_followup$bmi.y>=25))
bmi_followup$bmi25changeDown <- as.numeric((bmi_followup$bmi.x>=25)<(bmi_followup$bmi.y>=25))

summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))

bmi_followup$bmi30change <- as.numeric((bmi_followup$bmi.x>=30)!=(bmi_followup$bmi.y>=30))
bmi_followup$bmi30changeUp <- as.numeric((bmi_followup$bmi.x>=30)>(bmi_followup$bmi.y>=30))
bmi_followup$bmi30changeDown <- as.numeric((bmi_followup$bmi.x>=30)<(bmi_followup$bmi.y>=30))

summary(glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))
summary(glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))


#Method 1: Using the mids object

summary(pool(glm.mids(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids,family=binomial)))
summary(pool(glm.mids(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids,family=binomial)))

summary(pool(glm.mids(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids,family=binomial)))
summary(pool(glm.mids(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids,family=binomial)))

summary(pool(glm.mids(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids,family=binomial)))
summary(pool(glm.mids(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=base_data_mids,family=binomial)))


##Method 2: Rubin's rule for variances - taking into account the imputations - raw computations

sd25<-sd25up<-sd25down<-sd30<-sd30up<-sd30down<-
  matrix(nrow=length(summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup))$coefficients[,2]),ncol=20)

est25<-est25up<-est25down<-est30<-est30up<-est30down<-
  matrix(nrow=length(summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup))$coefficients[,1]),ncol=20)

for (i in 1:20){
  sd25[,i]<-summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,2]
  sd25up[,i]<-summary(glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,2]
  sd25down[,i]<-summary(glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,2]
  
  sd30[,i]<-summary(glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,2]
  sd30up[,i]<-summary(glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,2]
  sd30down[,i]<-summary(glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,2]
  
  est25[,i]<-summary(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,1]
  est25up[,i]<-summary(glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,1]
  est25down[,i]<-summary(glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,1]
  
  est30[,i]<-summary(glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,1]
  est30up[,i]<-summary(glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,1]
  est30down[,i]<-summary(glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=subset(bmi_followup,imputation==i),family=binomial))$coefficients[,1]
}

sd25Value<-sqrt(rowMeans(sd25^2)+rowSums((est25-rowMeans(est25))^2)/19) # sqrt(mean estimator variance + estimated imputation estimator variance): SE=sqrt(V_total).
sd25upValue<-sqrt(rowMeans(sd25up^2)+rowSums((est25up-rowMeans(est25up))^2)/19)
sd25downValue<-sqrt(rowMeans(sd25down^2)+rowSums((est25down-rowMeans(est25down))^2)/19)

sd30Value<-sqrt(rowMeans(sd30^2)+rowSums((est30-rowMeans(est30))^2)/19)
sd30upValue<-sqrt(rowMeans(sd30up^2)+rowSums((est30up-rowMeans(est30up))^2)/19)
sd30downValue<-sqrt(rowMeans(sd30down^2)+rowSums((est30down-rowMeans(est30down))^2)/19)


#Summary tables
change25_indicator_table <- data.frame("ests"<-coef(glm(bmi25change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial)),
                                      "sds"<- sd25Value)
change25up_indicator_table <- data.frame("ests"<-coef(glm(bmi25changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial)),
                                       "sds"<- sd25upValue)
change25down_indicator_table <- data.frame("ests"<-coef(glm(bmi25changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial)),
                                       "sds"<- sd25downValue)

change30_indicator_table <- data.frame("ests"<-coef(glm(bmi30change ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial)),
                                       "sds"<- sd30Value)
change30up_indicator_table <- data.frame("ests"<-coef(glm(bmi30changeUp ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial)),
                                         "sds"<- sd30upValue)
change30down_indicator_table <- data.frame("ests"<-c(coef(glm(bmi30changeDown ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*sample_weights.y-sample_weights.y,data=bmi_followup,family=binomial))),
                                           "sds"<- sd30downValue)




#####Tracking data for the followup CSS sample

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
(lambda <- bc$x[which.max(bc$y)])
new_model <- lm(((bmi^lambda-1)/lambda) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))

hist(((CSS_track$bmi^lambda-1)/lambda))
hist(simulate(new_model)$sim_1,breaks=20,xlim=c(0.830,0.860))
plot(fitted(new_model),residuals(new_model))


summary(new_model)


#Method 1: Mice-based inference for three models:

CSS_track_mids<-as.mids(CSS_track,.imp="impnr",.id="userid")

summary(pool(lm.mids(((bmi^lambda-1)/lambda) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track_mids)))
summary(pool(glm.mids((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track_mids,family=binomial)))
summary(pool(glm.mids((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track_mids,family=binomial)))


#Method 2: Adjusting for multiple imputations in the three kinds of models using Rubin's rules (raw calculations)

sd25<-sd30<-sdnum<-
  matrix(nrow=length(summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track))$coefficients[,2]),ncol=20)

est25<-est30<-estnum<-
  matrix(nrow=length(summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=CSS_track))$coefficients[,1]),ncol=20)


for (i in 1:20){
  sd25[,i]<-summary(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr==i)),family=binomial)$coefficients[,2]
  sd30[,i]<-summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr==i)),family=binomial)$coefficients[,2]
  sdnum[,i]<-summary(lm(((bmi^lambda-1)/lambda) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr==i)))$coefficients[,2]
  
  est25[,i]<-summary(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr==i)),family=binomial)$coefficients[,1]
  est30[,i]<-summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr==i)),family=binomial)$coefficients[,1]
  estnum[,i]<-summary(lm(((bmi^lambda-1)/lambda) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr==i)))$coefficients[,1]
}

sd25Value<-sqrt(rowMeans(sd25^2)+rowSums((est25-rowMeans(est25))^2)/19) # sqrt(mean estimator variance + estimated impnr estimator variance): SE=sqrt(V_total).
sd30Value<-sqrt(rowMeans(sd30^2)+rowSums((est30-rowMeans(est30))^2)/19)
sdnumValue<-sqrt(rowMeans(sdnum^2)+rowSums((estnum-rowMeans(estnum))^2)/19)

#Summary tables
bmi25_indicator_table_track <- data.frame("ests"<-coef(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))),
                                    "sds"<- sd25Value)
bmi30_indicator_table_track <- data.frame("ests"<-coef(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))),
                                    "sds"<- sd30Value)
bmi_num_table_track <- data.frame("ests"<-coef(glm(((bmi^lambda-1)/lambda) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weights-sample_weights,data=subset(CSS_track,impnr!=0))),
                            "sds"<- sdnumValue)

