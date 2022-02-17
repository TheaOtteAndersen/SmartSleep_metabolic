###### Analysis of associations between metabolic biomarkers and smartphone activity ######

#This script contains analyses of the multiply imputed data sets from the SmartSleep project.
#The multiple imputation results are combined using Rubin's rules.

library(dplyr)

#Reading in the data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
CSS <- subset(CSS,impnr!=0)

subject_tracking_clusters <- read.csv2("H:/SmartSleep backup IT issues/subject_tracking_clusters.csv")

#Looking at general patterns

#Cakculated bmi from imputed height and weight vs directly imputed bmi - which one to use?
bmi_ind <- which(abs(CSS$bmi-CSS$weight/((CSS$height/100)^2))>1)
cor(CSS$bmi[bmi_ind],(CSS$weight/((CSS$height/100)^2))[bmi_ind])
plot(CSS$bmi[bmi_ind],(CSS$weight/((CSS$height/100)^2))[bmi_ind])

hist(CSS$age,breaks=40)
hist(CSS$height,breaks=40)
hist(CSS$weight,breaks=40)


#Regression analyses


#Baseline data with self-reports
base_weights <- read.csv2("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/SmartSleepExpWeighted.csv")
load("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/SmartSleep Experiment/full_imp_base.RData")
base_data <- full_imp_Base
base_data$zipCode<-as.numeric(base_data$zipCode)

base_data$selfScore <- (base_data$mobileUseBeforeSleep=="5-7 times per week")*4+(base_data$mobileUseBeforeSleep=="2-4 times per week")*3+(base_data$mobileUseBeforeSleep=="Once a week")*3+(base_data$mobileUseBeforeSleep=="Every month or less")*2+(base_data$mobileUseBeforeSleep=="Never")*1+
  (base_data$mobileUseNight=="Every night or almost every night")*4+(base_data$mobileUseNight=="A few times a week")*3+(base_data$mobileUseNight=="A few times a month or less")*2+(base_data$mobileUseNight=="Never")*1+
  (base_data$mobileCheck==">20 times an hour")*4+(base_data$mobileCheck=="11-20 times an hour")*4+(base_data$mobileCheck=="5-10 times an hour")*3+(base_data$mobileCheck=="1-4 times an hour")*2+(base_data$mobileCheck=="Every 2nd hour")*2+(base_data$mobileCheck=="Several times a day")*1+(base_data$mobileCheck=="Once a day or less")*1+
  (base_data$pmpuScale<=14)*1+(base_data$pmpuScale>14 & base_data$pmpuScale<=16.76)*2+(base_data$pmpuScale>16.76 & base_data$pmpuScale<=19)*3+(base_data$pmpuScale>19)*4
base_data$selfScoreCat <- "1"
base_data$selfScoreCat[base_data$selfScore>=8]="2"
base_data$selfScoreCat[base_data$selfScore>=9.517]="3"
base_data$selfScoreCat[base_data$selfScore>=12]="4"

base_data$sample_weights <- as.numeric(unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")))$Weights) #use respondDate and variables that are relevant for weight sizes. Then take unique of the data.

#check for matching stuff
sum(base_data$responseDate[1:(nrow(base_data)/21)]!=unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")))$responseDate)
sum(base_data$age[1:(nrow(base_data)/21)]!=unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")))$age)
sum(base_data$gender[1:(nrow(base_data)/21)]!=unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")))$gender)
sum(base_data$zipCode[1:(nrow(base_data)/21)]!=unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")))$zipCode)
#So we have correct weights

base_data$CSSid <- unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","riseTime","sleepTime","education","occupation","mobilephone","pmpuScale")))$id

#We have as many non-na weights as there are rows in the base_weights data. This is a success.

#Base BMI - match with emailAddress



#BMI followup difference


#Tracking data

#Population sample
pop_weights <- read.csv2("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/ Population Sample/PopulationWeighted.csv")


#Followup sample
CSS_track <- inner_join(rename(CSS,userid=RespondKey),subject_tracking_clusters,by="userid")

CSS_track$selfScore <- (CSS_track$mobileUseBeforeSleep=="5-7 times per week")*4+(CSS_track$mobileUseBeforeSleep=="2-4 times per week")*3+(CSS_track$mobileUseBeforeSleep=="Once a week")*3+(CSS_track$mobileUseBeforeSleep=="Every month or less")*2+(CSS_track$mobileUseBeforeSleep=="Never")*1+
  (CSS_track$mobileUseNight=="Every night or almost every night")*4+(CSS_track$mobileUseNight=="A few times a week")*3+(CSS_track$mobileUseNight=="A few times a month or less")*2+(CSS_track$mobileUseNight=="Never")*1+
  (CSS_track$mobileCheck==">20 times an hour")*4+(CSS_track$mobileCheck=="11-20 times an hour")*4+(CSS_track$mobileCheck=="5-10 times an hour")*3+(CSS_track$mobileCheck=="1-4 times an hour")*2+(CSS_track$mobileCheck=="Every 2nd hour")*2+(CSS_track$mobileCheck=="Several times a day")*1+(CSS_track$mobileCheck=="Once a day or less")*1+
  (CSS_track$pmpuScale<=14)*1+(CSS_track$pmpuScale>14 & CSS_track$pmpuScale<=16.76)*2+(CSS_track$pmpuScale>16.76 & CSS_track$pmpuScale<=19)*3+(CSS_track$pmpuScale>19)*4
CSS_track$selfScoreCat <- "1"
CSS_track$selfScoreCat[CSS_track$selfScore>=8]="2"
CSS_track$selfScoreCat[CSS_track$selfScore>=9.517]="3"
CSS_track$selfScoreCat[CSS_track$selfScore>=12]="4"

CSS_weights <- read.csv2("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/FollowUpWeighted.csv")
CSS_track$sample_weight<-as.numeric(left_join(CSS_track,rename(CSS_weights,userid=RespondKey),by="userid")$Weights)

#Models
summary(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial))
hist(residuals(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial)))
plot(fitted(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial)),residuals(glm((bmi>=30) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial)))

summary(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial))
hist(residuals(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial)))
plot(fitted(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial)),residuals(glm((bmi>=25) ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track,family=binomial)))

summary(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track))
hist(residuals(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track)),breaks=40)
plot(fitted(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track)),residuals(lm(bmi ~ (cluster1prob+cluster2prob+cluster4prob+selfScoreCat+age+gender+education+occupation)*sample_weight,data=CSS_track)))
