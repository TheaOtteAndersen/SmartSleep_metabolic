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


## self-reported smartphone profiles
base_data$selfScore <- (base_data$mobileUseBeforeSleep=="5-7 times per week")*4+(base_data$mobileUseBeforeSleep=="2-4 times per week")*3+(base_data$mobileUseBeforeSleep=="Once a week")*3+(base_data$mobileUseBeforeSleep=="Every month or less")*2+(base_data$mobileUseBeforeSleep=="Never")*1+
  (base_data$mobileUseNight=="Every night or almost every night")*4+(base_data$mobileUseNight=="A few times a week")*3+(base_data$mobileUseNight=="A few times a month or less")*2+(base_data$mobileUseNight=="Never")*1+
  (base_data$mobileCheck==">20 times an hour")*4+(base_data$mobileCheck=="11-20 times an hour")*4+(base_data$mobileCheck=="5-10 times an hour")*3+(base_data$mobileCheck=="1-4 times an hour")*2+(base_data$mobileCheck=="Every 2nd hour")*2+(base_data$mobileCheck=="Several times a day")*1+(base_data$mobileCheck=="Once a day or less")*1+
  (base_data$pmpuScale<=14)*1+(base_data$pmpuScale>14 & base_data$pmpuScale<=16.76)*2+(base_data$pmpuScale>16.76 & base_data$pmpuScale<=19)*3+(base_data$pmpuScale>19)*4
base_data$selfScoreCat <- "1"
base_data$selfScoreCat[base_data$selfScore>=8]="2"
base_data$selfScoreCat[base_data$selfScore>=9.517]="3"
base_data$selfScoreCat[base_data$selfScore>=12]="4"
table(base_data$selfScoreCat)
