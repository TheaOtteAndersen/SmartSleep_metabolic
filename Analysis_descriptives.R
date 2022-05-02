#----------------------------------------------------------------------
## author: Thea Otte Andersen
## created: 17 February 2022
## Project: Night-time smartphone use and obesity and metabolic biomarkers - 2nd PhD paper 
#----------------------------------------------------------------------

## =============================================================================
## load libraries, import data, --------------------------------------------- ##

## Load libraries ----------------------------------------------------------- ##
library(Publish)
library(ggplot2)
library(readxl)
library(extrafont)
library(dplyr)
library(mice)

## =============================================================================
# Set path ------------------------------------------------------------------ ##
## Baseline Citizen Science Sample
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment")
base_data <- read.table("imp_Experiment.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
table(base_data$imputation)
base_data <- subset(base_data,imputation!=0)

#popu$sample_weights <- as.numeric(unique(left_join(popu[1:(nrow(popu)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","education","occupation","mobilephone","pmpuScale")))$Weights) #use respondDate and variables that are relevant for weight sizes. Then take unique of the data.
table(base_data$mobileUseNight, useNA = "always")

## exclude participants with no mobile phone N=1220
base_data <- subset(base_data, mobilephone!="No mobile phone")
table(base_data$mobilephone, useNA = "always")

## self-reported smartphone profiles
publish(univariateTable( ~ mobileUseBeforeSleep,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ mobileCheck,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ pmpuScale,data=base_data, column.percent=TRUE))
summary(base_data$pmpuScale)
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))

base_data$selfScore <- ((base_data$mobileUseBeforeSleep=="5-7 times per week")*4+(base_data$mobileUseBeforeSleep=="2-4 times per week")*3+(base_data$mobileUseBeforeSleep=="Once a week")*3+(base_data$mobileUseBeforeSleep=="Every month or less")*2+(base_data$mobileUseBeforeSleep=="Never")*1+
  (base_data$mobileUseNight=="Every night or almost every night")*4+(base_data$mobileUseNight=="A few times a week")*3+(base_data$mobileUseNight=="A few times a month or less")*2+(base_data$mobileUseNight=="Never")*1+
  (base_data$mobileCheck==">20 times per hour")*4+(base_data$mobileCheck=="11-20 times per hour")*4+(base_data$mobileCheck=="5-10 times per hour")*3+(base_data$mobileCheck=="1-4 times per hour")*2+(base_data$mobileCheck=="Every 2nd hour")*2+(base_data$mobileCheck=="Several times per day")*1+(base_data$mobileCheck=="Once a day")*1+
  (base_data$pmpuScale<14)*1+(base_data$pmpuScale>=14 & base_data$pmpuScale<17)*2+(base_data$pmpuScale>=17 & base_data$pmpuScale<19)*3+(base_data$pmpuScale>=19)*4)-4
summary(base_data$selfScore)
base_data$selfScore <- as.numeric(base_data$selfScore)

table(base_data$selfScore, useNA="always")
base_data$selfScoreCat<-NA
base_data$selfScoreCat[!is.na(base_data$selfScore)] <- "1"
base_data$selfScoreCat[base_data$selfScore>=4]="2"
base_data$selfScoreCat[base_data$selfScore>=6]="3"
base_data$selfScoreCat[base_data$selfScore>=8]="4"
table(base_data$selfScoreCat,useNA="always")/25
base_data$selfScoreCat <- as.factor(base_data$selfScoreCat)
prop.table(table(base_data$selfScoreCat,useNA="always")/25 )

## variables in table 1 (total)
## age
publish(univariateTable( ~ age,data=base_data, column.percent=TRUE))
#table(base_data$age)

## SD for pooled mean age

table(base_data$imputation)
sqrt( mean(aggregate(subset(base_data,imputation!=0)$age,by=list(subset(base_data,imputation!=0)$imputation),FUN=var)$x)+sum((aggregate(subset(base_data,imputation!=0)$age,by=list(subset(base_data,imputation!=0)$imputation),FUN=mean)$x-mean(subset(base_data,imputation!=0)$age))^2)/19
)


## gender
table(base_data$gender, useNA="always")/25
prop.table(table(base_data$gender)/25)
base_data$gender[base_data$gender=="Do not want to answer"] <- "Other"

## education
table(base_data$education, useNA="always")/25
prop.table(table(base_data$education))

## occupation
table(base_data$occupation, useNA="always")/25
prop.table(table(base_data$occupation)/25)

## table 1 according to subjective smartphone behavior profiles
table(base_data$selfScoreCat, useNA="always")/25
prop.table(table(base_data$selfScoreCat)/25)

## age ## how to calculate SD for each category??
publish(univariateTable(selfScoreCat ~ age,data=base_data, column.percent=TRUE))

## SD for risk profiles 1
sqrt( 
  mean(aggregate(subset(base_data,imputation!=0 & selfScoreCat==1)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==1)$imputation),FUN=var)$x)+sum((aggregate(subset(base_data,imputation!=0 & selfScoreCat==1)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==1)$imputation),FUN=mean)$x-mean(subset(base_data,imputation!=0 & selfScoreCat==1)$age))^2)/19
)

## SD for risk profile 2
sqrt( 
  mean(aggregate(subset(base_data,imputation!=0 & selfScoreCat==2)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==2)$imputation),FUN=var)$x)+sum((aggregate(subset(base_data,imputation!=0 & selfScoreCat==2)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==2)$imputation),FUN=mean)$x-mean(subset(base_data,imputation!=0 & selfScoreCat==2)$age))^2)/19
)

## SD for risk profile 3
sqrt( 
  mean(aggregate(subset(base_data,imputation!=0 & selfScoreCat==3)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==3)$imputation),FUN=var)$x)+sum((aggregate(subset(base_data,imputation!=0 & selfScoreCat==3)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==3)$imputation),FUN=mean)$x-mean(subset(base_data,imputation!=0 & selfScoreCat==3)$age))^2)/19
)
## SD for risk profile 4
sqrt( 
  mean(aggregate(subset(base_data,imputation!=0 & selfScoreCat==4)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==4)$imputation),FUN=var)$x)+sum((aggregate(subset(base_data,imputation!=0 & selfScoreCat==4)$age,by=list(subset(base_data,imputation!=0 & selfScoreCat==4)$imputation),FUN=mean)$x-mean(subset(base_data,imputation!=0 & selfScoreCat==4)$age))^2)/19
)

## gender
publish(univariateTable(selfScoreCat ~ gender,data=base_data, column.percent=TRUE))

## education
publish(univariateTable(selfScoreCat ~ education,data=base_data, column.percent=TRUE))

## occupation
publish(univariateTable(selfScoreCat ~ occupation,data=base_data, column.percent=TRUE))

## BMI
publish(univariateTable(selfScoreCat ~ bmi,data=base_data, column.percent=TRUE))
aggregate(base_data$bmi, by=list(base_data$selfScoreCat), FUN=median, na.rm=T)
aggregate(base_data$bmi, by=list(base_data$selfScoreCat), FUN=quantile, probs=0.75, na.rm=T)

# --------------------------------------------------------------------------- ##

## Population Sample
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample")
popu <- read.table("imp_population.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)

table(popu$imputation, useNA = "always")
popu <- subset(popu,imputation!=0)

## load tracking data subjects
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data")
#tracking <- read.table("subject_tracking_clusters.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
tracking_popu <- subset(subject_tracking_clusters, random==1)

popuAlt <- merge(tracking_popu, popu, by="userid")
popuAlt <- subset(popuAlt,imputation!=0)

## distribution in only survey data
## age distribution
publish(univariateTable( ~ agePNR,data=popu, column.percent=TRUE))

## SD for pooled mean age
sqrt( mean(aggregate(subset(popu,imputation!=0)$agePNR,by=list(subset(popu,imputation!=0)$imputation),FUN=var)$x)+sum((aggregate(subset(popu,imputation!=0)$agePNR,by=list(subset(popu,imputation!=0)$imputation),FUN=mean)$x-mean(subset(popu,imputation!=0)$agePNR))^2)/19
)

## gender
table(popu$Gender, useNA="always")/20
prop.table(table(popu$Gender))

## education
table(popu$education, useNA="always")/20
prop.table(table(popu$education))

## occupation
table(popu$occupation, useNA="always")/20
prop.table(table(popu$occupation))

## bmi
table(popu$bmi)

popu$bmiTEST <- as.numeric(popu$bmi)
table(popu$bmiTEST)
mean(popu$bmiTEST)
popu$bmi[popu$bmi==0] <- NA
popu$bmi[popu$height<100] <- (popu$weight/(((popu$height+100)/100)^2))[popu$height<100]
popu$height[popu$height<100 & popu$imputation!=0] <- popu$height[popu$height<100 & popu$imputation!=0]+100 
popu$bmi[popu$height==popu$weight]<-NA
popu$bmi[popu$bmi<14]<-NA
table(popu$bmi)

## self-reported smartphone behavior
table(popu$selfScore)


## self-reported smartphone profiles
summary(base_data$mobileCheck)
popu$selfScore <- (popu$mobileUseBeforeSleep=="5-7 times per week")*4+(popu$mobileUseBeforeSleep=="2-4 times per week")*3+(popu$mobileUseBeforeSleep=="Once a week")*3+(popu$mobileUseBeforeSleep=="Every month or less")*2+(popu$mobileUseBeforeSleep=="Never")*1+
  (popu$mobileUseNight=="Every night or almost every night")*4+(popu$mobileUseNight=="A few times a week")*3+(popu$mobileUseNight=="A few times a month or less")*2+(popu$mobileUseNight=="Never")*1+
  (popu$mobileCheck==">20 times an hour")*4+(popu$mobileCheck=="11-20 times an hour")*4+(popu$mobileCheck=="5-10 times an hour")*3+(popu$mobileCheck=="1-4 times an hour")*2+(popu$mobileCheck=="Every 2nd hour")*2+(popu$mobileCheck=="Several times a day")*1+(popu$mobileCheck=="Once a day or less")*1+
  (popu$pmpuScale<14)*1+(popu$pmpuScale>=14 & popu$pmpuScale<17)*2+(popu$pmpuScale>=17 & popu$pmpuScale<19)*3+(popu$pmpuScale>=19)*4
summary(popu$selfScore)

table(popu$selfScore, useNA="always")
popu$selfScoreCat<-NA
popu$selfScoreCat[!is.na(popu$selfScore)] <- "1"
popu$selfScoreCat[popu$selfScore>=8]="2"
popu$selfScoreCat[popu$selfScore>=10]="3"
popu$selfScoreCat[popu$selfScore>=12]="4"
table(popu$selfScoreCat,useNA="always")/20 
popu$selfScoreCat <- as.factor(popu$selfScoreCat)
prop.table(table(popu$selfScoreCat,useNA="always"))

# --------------------------------------------------------------------------- ##
  
#Reading in the data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
CSS <- subset(CSS,impnr!=0)
CSS$userid <- CSS$RespondKey

tracking_CSS <- subset(subject_tracking_clusters, followup==1)
table(subject_tracking_clusters$followup)
  
CSSAlt <- merge(tracking_CSS, CSS, by="userid")


## distributin only in survey data

## age
publish(univariateTable( ~ age,data=CSS, column.percent=TRUE))

## SD for pooled mean age
age_mean <- with(CSS, expr=c("SD_age"=stats::sd(age)))
withPool_MI(age_mean)

sqrt( mean(aggregate(subset(CSS,impnr!=0)$age,by=list(subset(CSS,impnr!=0)$impnr),FUN=var)$x)+sum((aggregate(subset(CSS,impnr!=0)$age,by=list(subset(CSS,impnr!=0)$impnr),FUN=mean)$x-mean(subset(CSS,impnr!=0)$age))^2)/19
)


## gender
table(CSS$gender, useNA="always")/20
prop.table(table(CSS$gender))

## education
table(CSS$education, useNA="always")/20
prop.table(table(CSS$education))

## occupation
table(CSS$occupation, useNA="always")
prop.table(table(CSS$education))

## risk profiles of self-reported smartphone behavior
## self-reported smartphone profiles
table(CSS$mobileCheck)
CSS$selfScore <- (CSS$mobileUseBeforeSleep=="5-7 times per week")*4+(CSS$mobileUseBeforeSleep=="2-4 times per week")*3+(CSS$mobileUseBeforeSleep=="Once a week")*3+(CSS$mobileUseBeforeSleep=="Every month or less")*2+(CSS$mobileUseBeforeSleep=="Never")*1+
  (CSS$mobileUseNight=="Every night or almost every night")*4+(CSS$mobileUseNight=="A few times a week")*3+(CSS$mobileUseNight=="A few times a month or less")*2+(CSS$mobileUseNight=="Never")*1+
  (CSS$mobileCheck==">20 times an hour")*4+(CSS$mobileCheck=="11-20 times an hour")*4+(CSS$mobileCheck=="5-10 times an hour")*3+(CSS$mobileCheck=="1-4 times an hour")*2+(CSS$mobileCheck=="Every 2nd hour")*2+(CSS$mobileCheck=="Several times a day")*1+(CSS$mobileCheck=="Once a day or less")*1+
  (CSS$pmpuScale<14)*1+(CSS$pmpuScale>=14 & CSS$pmpuScale<17)*2+(CSS$pmpuScale>=17 & CSS$pmpuScale<19)*3+(CSS$pmpuScale>=19)*4
summary(CSS$selfScore)

table(CSS$selfScore, useNA="always")
CSS$selfScoreCat<-NA
CSS$selfScoreCat[!is.na(CSS$selfScore)] <- "1"
CSS$selfScoreCat[CSS$selfScore>=8]="2"
CSS$selfScoreCat[CSS$selfScore>=10]="3"
CSS$selfScoreCat[CSS$selfScore>=12]="4"
table(CSS$selfScoreCat,useNA="always")/20 
CSS$selfScoreCat <- as.factor(CSS$selfScoreCat)
prop.table(table(CSS$selfScoreCat,useNA="always"))


# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##
## sammenligning med den danske population

## baseline citizen Science sample
## load baseline data and only keep imputation = 0
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment")
base <- read.table("imp_Experiment.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
base <- subset(base,imputation==0)

## age groups
base$ageCat[base$age>=16 & base$age<20] <- "16-19"
base$ageCat[base$age>=20 & base$age<25] <- "20-24"
base$ageCat[base$age>=25 & base$age<30] <- "25-29"
base$ageCat[base$age>=30 & base$age<35] <- "30-34"
base$ageCat[base$age>=35 & base$age<40] <- "35-39"
base$ageCat[base$age>=40 & base$age<45] <- "40-44"
base$ageCat[base$age>=45 & base$age<50] <- "45-49"
base$ageCat[base$age>=50 & base$age<55] <- "50-54"
base$ageCat[base$age>=55 & base$age<60] <- "55-59"
base$ageCat[base$age>=60 & base$age<65] <- "60-64"
base$ageCat[base$age>=65 & base$age<70] <- "65-69"
base$ageCat[base$age>=70 & base$age<75] <- "70-74"
base$ageCat[base$age>=75 & base$age<80] <- "75-79"
base$ageCat[base$age>=80] <- "80+"

publish(univariateTable( ~ ageCat,data=base, column.percent=TRUE))

## region
publish(univariateTable( ~ regionD,data=base, column.percent=TRUE))

## gender
publish(univariateTable( ~ gender,data=base, column.percent=TRUE))


## Population sample

## Population Sample
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample")
po <- read.table("imp_population.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)

table(po$imputation, useNA = "always")
po <- subset(po,imputation==0)

## age groups
po$ageCat[po$agePNR>=18 & po$agePNR<23] <- "18-22"
po$ageCat[po$agePNR>=23 & po$agePNR<28] <- "23-27"
po$ageCat[po$agePNR>=28 & po$agePNR<33] <- "28-32"
po$ageCat[po$agePNR>=33 & po$agePNR<38] <- "33-37"
po$ageCat[po$agePNR>=38 & po$agePNR<43] <- "38-42"
po$ageCat[po$agePNR>=43 & po$agePNR<48] <- "43-47"
po$ageCat[po$agePNR>=48] <- "48+"

publish(univariateTable( ~ ageCat,data=po, column.percent=TRUE))

publish(univariateTable( ~ Gender,data=po, column.percent=TRUE))

## load weighted datasæt where regionD is
setwd("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Population Sample/")
poW <- read.table("PopulationWeighted.csv", header=TRUE, fill=TRUE, quote = "", sep=";", stringsAsFactors = TRUE)

publish(univariateTable( ~regionD ,data=poW, column.percent=TRUE))

## Citizen science follow-up sample

setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
CS <- subset(CS,impnr==0)

## age 

CS$ageCat[CS$age>=16 & CS$age<21] <- "16-20"
CS$ageCat[CS$age>=21 & CS$age<26] <- "21-25"
CS$ageCat[CS$age>=26 & CS$age<31] <- "26-30"
CS$ageCat[CS$age>=31 & CS$age<36] <- "31-35"
CS$ageCat[CS$age>=36 & CS$age<41] <- "36-40"
CS$ageCat[CS$age>=41 & CS$age<46] <- "41-45"
CS$ageCat[CS$age>=46 & CS$age<51] <- "46-50"
CS$ageCat[CS$age>=51 & CS$age<56] <- "51-55"
CS$ageCat[CS$age>=56 & CS$age<61] <- "56-60"
CS$ageCat[CS$age>=61 & CS$age<66] <- "61-65"
CS$ageCat[CS$age>=66] <- "+65"
publish(univariateTable( ~ ageCat,data=CS, column.percent=TRUE))

publish(univariateTable( ~ gender,data=CS, column.percent=TRUE))

## load weighted datasæt where regionD is
setwd("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/")
CSW <- read.table("FollowUpWeighted.csv", header=TRUE, fill=TRUE, quote = "", sep=";", stringsAsFactors = TRUE)


publish(univariateTable( ~ regionD,data=CSW, column.percent=TRUE))

