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
table(base_data$selfScoreCat, base_data$gender)/25

## education
table(base_data$selfScoreCat, base_data$education)/25
publish(univariateTable(selfScoreCat ~ education,data=base_data, column.percent=TRUE))

## occupation
publish(univariateTable(selfScoreCat ~ occupation,data=base_data, column.percent=TRUE))

## BMI 
table(base_data$bmi)
base_data$bmi<- gsub(",", ".", base_data$bmi)
base_data$bmi <- as.numeric(base_data$bmi)

publish(univariateTable(selfScoreCat ~ bmi,data=base_data, column.percent=TRUE))
aggregate(base_data$bmi, by=list(base_data$selfScoreCat), FUN=median, na.rm=T)
aggregate(base_data$bmi, by=list(base_data$selfScoreCat), FUN=quantile, probs=0.25, na.rm=T)
aggregate(base_data$bmi, by=list(base_data$selfScoreCat), FUN=quantile, probs=0.75, na.rm=T)

## bmi >=25
base_data$bmiCat[base_data$bmi<25] <- "<25"
base_data$bmiCat[base_data$bmi>=25 & base_data$bmi<30] <- "25-30"
base_data$bmiCat[base_data$bmi>=30] <- ">=30"
table(base_data$bmiCat)/25
prop.table(table(base_data$bmiCat))
table(base_data$selfScoreCat, base_data$bmiCat)/25
publish(univariateTable(selfScoreCat ~ bmiCat,data=base_data, column.percent=TRUE))
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
table(popu$Gender, useNA="always")/25
prop.table(table(popu$Gender))

## education
table(popu$education, useNA="always")/25
prop.table(table(popu$education))

## occupation
table(popu$occupation, useNA="always")/25
prop.table(table(popu$occupation))

## bmi
table(popu$bmi)
popu$bmi<- gsub(",", ".", popu$bmi)
popu$bmi <- as.numeric(popu$bmi)

popu$bmiCat[popu$bmi<25] <- "<25"
popu$bmiCat[popu$bmi>=25 & popu$bmi<30] <- "25-30"
popu$bmiCat[popu$bmi>=30] <- ">=30"
table(popu$bmiCat, useNA="always")/25
prop.table(table(popu$bmiCat, useNA="always")/25)

## self-reported smartphone profiles
summary(base_data$mobileCheck)
popu$selfScore <- ((popu$mobileUseBeforeSleep=="5-7 times per week")*4+(popu$mobileUseBeforeSleep=="2-4 times per week")*3+(popu$mobileUseBeforeSleep=="Once a week")*3+(popu$mobileUseBeforeSleep=="Every month or less")*2+(popu$mobileUseBeforeSleep=="Never")*1+
  (popu$mobileUseNight=="Every night or almost every night")*4+(popu$mobileUseNight=="A few times a week")*3+(popu$mobileUseNight=="A few times a month or less")*2+(popu$mobileUseNight=="Never")*1+
  (popu$mobileCheck==">20 times an hour")*4+(popu$mobileCheck=="11-20 times an hour")*4+(popu$mobileCheck=="5-10 times an hour")*3+(popu$mobileCheck=="1-4 times an hour")*2+(popu$mobileCheck=="Every 2nd hour")*2+(popu$mobileCheck=="Several times a day")*1+(popu$mobileCheck=="Once a day or less")*1+
  (popu$pmpuScale<14)*1+(popu$pmpuScale>=14 & popu$pmpuScale<17)*2+(popu$pmpuScale>=17 & popu$pmpuScale<19)*3+(popu$pmpuScale>=19)*4)-4
summary(popu$selfScore)

table(popu$selfScore, useNA="always")
popu$selfScoreCat[!is.na(popu$selfScore)] <- "1"
popu$selfScoreCat[popu$selfScore>=4]="2"
popu$selfScoreCat[popu$selfScore>=6]="3"
popu$selfScoreCat[popu$selfScore>=8]="4"
table(popu$selfScoreCat,useNA="always")/25 
popu$selfScoreCat <- as.factor(popu$selfScoreCat)
prop.table(table(popu$selfScoreCat,useNA="always"))

# --------------------------------------------------------------------------- ##
  
#Reading in the data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
table(CSS$imputation)
CSS <- subset(CSS,imputation!=0)

tracking_CSS <- subset(subject_tracking_clusters, followup==1)
table(subject_tracking_clusters$followup)
  
CSSAlt <- merge(tracking_CSS, CSS, by="userid")


## distributin only in survey data

## age
publish(univariateTable( ~ age,data=CSS, column.percent=TRUE))

## SD for pooled mean age
age_mean <- with(CSS, expr=c("SD_age"=stats::sd(age)))
#withPool_MI(age_mean)

sqrt( mean(aggregate(subset(CSS,imputation!=0)$age,by=list(subset(CSS,imputation!=0)$imputation),FUN=var)$x)+sum((aggregate(subset(CSS,imputation!=0)$age,by=list(subset(CSS,imputation!=0)$imputation),FUN=mean)$x-mean(subset(CSS,imputation!=0)$age))^2)/19
)

## gender
table(CSS$gender, useNA="always")/25
prop.table(table(CSS$gender))

## education
table(CSS$education, useNA="always")/25
prop.table(table(CSS$education))

## occupation
table(CSS$occupation, useNA="always")/25
prop.table(table(CSS$occupation))

## risk profiles of self-reported smartphone behavior
## self-reported smartphone profiles
table(CSS$mobileCheck)
CSS$selfScore <- ((CSS$mobileUseBeforeSleep=="5-7 times per week")*4+(CSS$mobileUseBeforeSleep=="2-4 times per week")*3+(CSS$mobileUseBeforeSleep=="Once a week")*3+(CSS$mobileUseBeforeSleep=="Every month or less")*2+(CSS$mobileUseBeforeSleep=="Never")*1+
  (CSS$mobileUseNight=="Every night or almost every night")*4+(CSS$mobileUseNight=="A few times a week")*3+(CSS$mobileUseNight=="A few times a month or less")*2+(CSS$mobileUseNight=="Never")*1+
  (CSS$mobileCheck==">20 times an hour")*4+(CSS$mobileCheck=="11-20 times an hour")*4+(CSS$mobileCheck=="5-10 times an hour")*3+(CSS$mobileCheck=="1-4 times an hour")*2+(CSS$mobileCheck=="Every 2nd hour")*2+(CSS$mobileCheck=="Several times a day")*1+(CSS$mobileCheck=="Once a day or less")*1+
  (CSS$pmpuScale<14)*1+(CSS$pmpuScale>=14 & CSS$pmpuScale<17)*2+(CSS$pmpuScale>=17 & CSS$pmpuScale<19)*3+(CSS$pmpuScale>=19)*4)-4
summary(CSS$selfScore)

table(CSS$selfScore, useNA="always")
CSS$selfScoreCat[!is.na(CSS$selfScore)] <- "1"
CSS$selfScoreCat[CSS$selfScore>=4]="2"
CSS$selfScoreCat[CSS$selfScore>=6]="3"
CSS$selfScoreCat[CSS$selfScore>=8]="4"
table(CSS$selfScoreCat,useNA="always")/25 
prop.table(table(CSS$selfScoreCat,useNA="always"))

## bmi
table(CSS$bmi, useNA="always")
CSS$bmiCat[CSS$bmi<25] <- "<25"
CSS$bmiCat[CSS$bmi>=25 & popu$bmi<30] <- "25-30"
CSS$bmiCat[CSS$bmi>=30] <- ">=30"
table(CSS$bmiCat, useNA="always")/25
prop.table(table(CSS$bmiCat, useNA="always")/25)

# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##
## sammenligning med den danske population

## baseline citizen Science sample
## load baseline data and only keep imputation = 0
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment")
base <- read.table("imp_Experiment.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
base <- subset(base,imputation==0)

prop.table(table(base$ageCat))
prop.table(table(base$ageCat, useNA="always"))

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
prop.table(table(po$ageCat, useNA="always"))

## sex
publish(univariateTable( ~ Gender,data=po, column.percent=TRUE))

prop.table(table(popu$regionD))


## Citizen science follow-up sample

setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
table(CS$imputation)
CS <- subset(CS,imputation==0)

## age 
table(CS$age)
prop.table(table(CS$ageCat, useNA="always"))

publish(univariateTable( ~ gender,data=CS, column.percent=TRUE))
prop.table(table(CS$sex, useNA="always"))
prop.table(table(CS$regionD, useNA="always"))

