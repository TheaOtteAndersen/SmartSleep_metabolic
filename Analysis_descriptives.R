#----------------------------------------------------------------------
## author: Thea Otte Andersen
## created: 17 February 2022
## Project: Night-time smartphone use and obesity and metabolic biomarkers - 2nd PhD paper 
#----------------------------------------------------------------------

## =============================================================================
## load libraries, import data, --------------------------------------------- ##

## Load libraries ----------------------------------------------------------- ##
library("Publish")
library(ggplot2)
library(readxl)
library(extrafont)
library(dplyr)

## =============================================================================
# Set path ------------------------------------------------------------------ ##

## Baseline Citizen Science Sample
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment")
base_data <- read.table("imp_Experiment.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
table(base_data$imp_nr)
base_data <- subset(base_data,imp_nr!=0)

## include weights ## OBS virker ikke
base_weights <- read.csv2("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/SmartSleepExpWeighted.csv")

#base_data$zipCode<-as.numeric(base_data$zipCode)
#base_data$riseTime<-as.numeric(base_data$riseTime)
#base_data$sleepTime<-as.numeric(base_data$sleepTime)
#base_weights$riseTime <- as.numeric(base_weights$riseTime)
#base_weights$sleepTime <- as.numeric(base_weights$sleepTime)
#table(base_weights$riseTime)

#base_data$sample_weights <- as.numeric(unique(left_join(base_data[1:(nrow(base_data)/21),],base_weights[,c("id","responseDate","Weights","age","gender","zipCode","education","occupation","mobilephone","pmpuScale")],by=c("responseDate","age","gender","zipCode","education","occupation","mobilephone","pmpuScale")))$Weights) #use respondDate and variables that are relevant for weight sizes. Then take unique of the data.
table(base_data$mobileUseNight, useNA = "always")

## exclude participants with no mobile phone
base_data <- subset(base_data, mobilephone!="No mobile phone")
  
table(base_data$mobileCheck, useNA = "always")
## self-reported smartphone profiles
base_data$selfScore <- (base_data$mobileUseBeforeSleep=="5-7 times per week")*4+(base_data$mobileUseBeforeSleep=="2-4 times per week")*3+(base_data$mobileUseBeforeSleep=="Once a week")*3+(base_data$mobileUseBeforeSleep=="Every month or less")*2+(base_data$mobileUseBeforeSleep=="Never")*1+
  (base_data$mobileUseNight=="Every night or almost every night")*4+(base_data$mobileUseNight=="A few times a week")*3+(base_data$mobileUseNight=="A few times a month or less")*2+(base_data$mobileUseNight=="Never")*1+
  (base_data$mobileCheck==">20 times per hour")*4+(base_data$mobileCheck=="11-20 times per hour")*4+(base_data$mobileCheck=="5-10 times per hour")*3+(base_data$mobileCheck=="1-4 times per hour")*2+(base_data$mobileCheck=="Every 2nd hour")*2+(base_data$mobileCheck=="Several times per day")*1+(base_data$mobileCheck=="Once a day")*1+
  (base_data$pmpuScale<=14)*1+(base_data$pmpuScale>14 & base_data$pmpuScale<=17)*2+(base_data$pmpuScale>17 & base_data$pmpuScale<=19)*3+(base_data$pmpuScale>19)*4
summary(base_data$selfScore)

table(base_data$selfScore, useNA="always")
base_data$selfScoreCat[base_data$selfScore<=7] <- "1"
base_data$selfScoreCat[base_data$selfScore>7 & base_data$selfScore<=9]="2"
base_data$selfScoreCat[base_data$selfScore>9 & base_data$selfScore<=12]="3"
base_data$selfScoreCat[base_data$selfScore>12]="4"
base_data$selfScoreCat <- as.factor(base_data$selfScoreCat)
table(base_data$selfScoreCat)

## variables in table 1 (total)
## age
publish(univariateTable( ~ age,data=base_data, column.percent=TRUE))
table(base_data$age)

## gender
publish(univariateTable( ~ gender,data=base_data, column.percent=TRUE))
base_data$gender[base_data$gender=="Do not want to answer"] <- "Other"

## education
publish(univariateTable( ~ education,data=base_data, column.percent=TRUE))

## occupation
publish(univariateTable( ~ occupation,data=base_data, column.percent=TRUE))

## table 1 according to subjective smartpone behavior profiles
publish(univariateTable( ~ selfScoreCat,data=base_data, column.percent=TRUE))

publish(univariateTable(selfScoreCat ~ age,data=base_data, column.percent=TRUE))

publish(univariateTable(selfScoreCat ~ gender,data=base_data, column.percent=TRUE))

publish(univariateTable(selfScoreCat ~ education,data=base_data, column.percent=TRUE))

publish(univariateTable(selfScoreCat ~ occupation,data=base_data, column.percent=TRUE))

# --------------------------------------------------------------------------- ##

## Population Sample
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample")
popu <- read.table("imp_population_w0.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)

table(popu$impnr, useNA = "always")
popu <- subset(popu,impnr!=0)

popu$userid <- popu$RespondKey
table(popu$userid)

setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data")
tracking <- read.table("subject_tracking_clusters.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
tracking_popu <- subset(tracking, random==1)
table(tracking$random)

popuAlt <- merge(tracking_popu, popu, by="userid")

# --------------------------------------------------------------------------- ##
  
#Reading in the data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
CSS <- subset(CSS,impnr!=0)
CSS$userid <- CSS$RespondKey

tracking_CSS <- subset(tracking, followup==1)
table(tracking$followup)
  
CSSAlt <- merge(tracking_CSS, CSS, by="userid")
