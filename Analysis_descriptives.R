

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

## rubins rules til mean(sd)
estimate.pooler.unitary <- function(coef,sd){
  n_col <- length(coef)
  
  coefs <- mean(coef)
  sds <- sqrt(mean(sd^2)+(1+1/n_col)*sum((coef-coefs)^2)/(n_col-1))
  df <- data.frame("estimate"=coefs,"sd"=sds,"lower.CI"=coefs-1.96*sds,"upper.CI"=coefs+1.96*sds)
  return(df)
}

## self-reported night-time smartphone use and smartphone use before sleep onset
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))
base_data$mobileUseNight <- factor(base_data$mobileUseNight, levels = c("Never", "A few times a month or less", "A few times a week", "Every night or almost every night"))

## variables in table 1 (total)
## age
publish(univariateTable( ~ age,data=base_data, column.percent=TRUE))
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i)$age)
  sdev[i] <- sd(subset(base_data,imputation==i)$age) 
}

estimate.pooler.unitary(m,sdev)

## SD for pooled mean age
## SD for Never
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="Never")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="Never")$age) 
}

estimate.pooler.unitary(m,sdev)

## SD for a few times a month
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="A few times a month or less")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="A few times a month or less")$age) 
}

estimate.pooler.unitary(m,sdev)

## SD for "a few times a week"
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="A few times a week")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="A few times a week")$age) 
}

estimate.pooler.unitary(m,sdev)

## every night or almost every night
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="Every night or almost every night")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="Every night or almost every night")$age) 
}

estimate.pooler.unitary(m,sdev)

##

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

## table 1 according to night-time smartphone use
table(base_data$mobileUseNight, useNA="always")/25
prop.table(table(base_data$mobileUseNight)/25)

## age ## how to calculate SD for each category??
publish(univariateTable(mobileUseNight ~ age,data=base_data, column.percent=TRUE))

## SD for Never
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="Never")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="Never")$age) 
}

estimate.pooler.unitary(m,sdev)

## age and SD for every month
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="A few times a month or less")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="A few times a month or less")$age) 
}

estimate.pooler.unitary(m,sdev)

## every week
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="A few times a week")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="A few times a week")$age) 
}

estimate.pooler.unitary(m,sdev)

## every night or almost every night
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="Every night or almost every night")$age)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="Every night or almost every night")$age) 
}

estimate.pooler.unitary(m,sdev)

## gender
publish(univariateTable(mobileUseNight ~ gender,data=base_data, column.percent=TRUE))
table(base_data$mobileUseNight, base_data$gender)/25

## education
table(base_data$mobileUseNight, base_data$education)/25
publish(univariateTable(mobileUseNight ~ education,data=base_data, column.percent=TRUE))

## occupation
table(base_data$mobileUseNight, base_data$occupation)/25
publish(univariateTable(mobileUseNight ~ occupation,data=base_data, column.percent=TRUE))

## BMI 
table(base_data$bmi, useNA="always")
base_data$bmi<- gsub(",", ".", base_data$bmi)
base_data$bmi <- as.numeric(base_data$bmi)

publish(univariateTable(mobileUseNight ~ bmi,data=base_data, column.percent=TRUE))

## SD for Never (bmi)
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="Never")$bmi)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="Never")$bmi) 
}

estimate.pooler.unitary(m,sdev)

## SD for every month (bmi)
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="A few times a month or less")$bmi)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="A few times a month or less")$bmi) 
}

estimate.pooler.unitary(m,sdev)

## SD for every week (bmi)
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="A few times a week")$bmi)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="A few times a week")$bmi) 
}

estimate.pooler.unitary(m,sdev)

## SD for every night or almost every night (bmi)
m <- numeric(0)
sdev <- numeric(0)
for (i in 1:25){
  m[i] <- mean(subset(base_data,imputation==i & mobileUseNight=="Every night or almost every night")$bmi)
  sdev[i] <- sd(subset(base_data,imputation==i & mobileUseNight=="Every night or almost every night")$bmi) 
}

estimate.pooler.unitary(m,sdev)


## bmi >=25
base_data$bmiCat[base_data$bmi<25] <- "<25"
base_data$bmiCat[base_data$bmi>=25 & base_data$bmi<30] <- "25-30"
base_data$bmiCat[base_data$bmi>=30] <- ">=30"
table(base_data$bmiCat)/25
prop.table(table(base_data$bmiCat))
table(base_data$mobileUseNight, base_data$bmiCat)/25
publish(univariateTable(mobileUseNight ~ bmiCat,data=base_data, column.percent=TRUE))

publish(univariateTable(mobileUseNight ~ bmi,data=base_data, column.percent=TRUE))

# --------------------------------------------------------------------------- ##


## Population Sample
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample")
popu <- read.table("imp_population.csv", header=TRUE, fill=TRUE, stringsAsFactors = TRUE, strip.white = TRUE, sep=";")

table(popu$imputation, useNA = "always")
popu <- subset(popu,imputation!=0)

## load tracking data subjects
#setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data")
#tracking <- read.table("subject_tracking_clusters.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)
#tracking_popu <- subset(subject_tracking_clusters, random==1)

#popuAlt <- merge(tracking_popu, popu, by="userid")
#popuAlt <- subset(popuAlt,imputation!=0)

## distribution in only survey data
## age distribution
publish(univariateTable( ~ agePNR,data=popu, column.percent=TRUE))

## SD for pooled mean age
sqrt( mean(aggregate(subset(popu,imputation!=0)$agePNR,by=list(subset(popu,imputation!=0)$imputation),FUN=var)$x)+sum((aggregate(subset(popu,imputation!=0)$agePNR,by=list(subset(popu,imputation!=0)$imputation),FUN=mean)$x-mean(subset(popu,imputation!=0)$agePNR))^2)/19
)


## sex
table(popu$sex, useNA="always")/25
prop.table(table(popu$sex))

## education
table(popu$education, useNA="always")/25
prop.table(table(popu$education))

## occupation
table(popu$occupation, useNA="always")/25
prop.table(table(popu$occupation))

## night-time smartphone use
table(popu$mobileUseNight, useNA="always")/25
publish(univariateTable( ~ mobileUseNight,data=popu, column.percent=TRUE))

## bmi
popu$bmi <- as.numeric(popu$bmi)
table(popu$bmi)
publish(univariateTable(mobileUseNight ~ bmi,data=popu, column.percent=TRUE))

##---------------------------------------------------------------------------##
## confounders according to night-time smartphone use

## age
publish(univariateTable(mobileUseNight ~ agePNR,data=popu, column.percent=TRUE))

## SD for mobileUseNight=Never
sqrt( 
  mean(aggregate(subset(popu,imputation!=0 & mobileUseNight=="Never")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="Never")$imputation),FUN=var)$x)+sum((aggregate(subset(popu,imputation!=0 & mobileUseNight=="Never")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="Never")$imputation),FUN=mean)$x-mean(subset(popu,imputation!=0 & mobileUseNight=="Never")$agePNR))^2)/25
)

## SD formobileUseNight=A few times a month or less
sqrt( 
  mean(aggregate(subset(popu,imputation!=0 & mobileUseNight=="A few times a month or less")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="A few times a month or less")$imputation),FUN=var)$x)+sum((aggregate(subset(popu,imputation!=0 & mobileUseNight=="A few times a month or less")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="A few times a month or less")$imputation),FUN=mean)$x-mean(subset(popu,imputation!=0 & mobileUseNight=="A few times a month or less")$agePNR))^2)/25
)

## SD for mobileUseNight=A few times a week
sqrt( 
  mean(aggregate(subset(popu,imputation!=0 & mobileUseNight=="A few times a week")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="A few times a week")$imputation),FUN=var)$x)+sum((aggregate(subset(popu,imputation!=0 & mobileUseNight=="A few times a week")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="A few times a week")$imputation),FUN=mean)$x-mean(subset(popu,imputation!=0 & mobileUseNight=="A few times a week")$agePNR))^2)/25
)

## SD for mobileUseNight=Every night or almost every night
sqrt( 
  mean(aggregate(subset(popu,imputation!=0 & mobileUseNight=="Every night or almost every night")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="Every night or almost every night")$imputation),FUN=var)$x)+sum((aggregate(subset(popu,imputation!=0 & mobileUseNight=="Every night or almost every night")$agePNR,by=list(subset(popu,imputation!=0 & mobileUseNight=="Every night or almost every night")$imputation),FUN=mean)$x-mean(subset(popu,imputation!=0 & mobileUseNight=="Every night or almost every night")$agePNR))^2)/25
)

## sex
publish(univariateTable(mobileUseNight ~ sex,data=popu, column.percent=TRUE))
table(popu$mobileUseNight, popu$sex)/25

## education
table(popu$mobileUseNight, popu$education)/25
publish(univariateTable(mobileUseNight ~ education,data=popu, column.percent=TRUE))

## occupation
table(popu$mobileUseNight, popu$occupation)/25
publish(univariateTable(mobileUseNight ~ occupation,data=popu, column.percent=TRUE))

# --------------------------------------------------------------------------- ##
## clinical sample

setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample")
clinical <- read.table("imp_clinical.csv", header=TRUE, fill=TRUE, sep=";", stringsAsFactors = TRUE, strip.white = TRUE)

## imputation
table(clinical$imputation, useNA = "always")
clinical <- subset(clinical,imputation!=0)

#names(clinical)
## distribution in only survey data
## age distribution
publish(univariateTable( ~ age,data=clinical, column.percent=TRUE))

## SD for pooled mean age
sqrt( mean(aggregate(subset(clinical,imputation!=0)$age,by=list(subset(clinical,imputation!=0)$imputation),FUN=var)$x)+sum((aggregate(subset(clinical,imputation!=0)$age,by=list(subset(clinical,imputation!=0)$imputation),FUN=mean)$x-mean(subset(clinical,imputation!=0)$age))^2)/19
)

## education
table(clinical$education, useNA="always")/25
prop.table(table(clinical$education))

#recategorise education
table(clinical$education, useNA="always")
clinical$educationCat[clinical$education=="long cycle higher education"] <- "long cycle higher education"
clinical$educationCat[clinical$education=="medium cycle higher education"] <- "medium cycle higher education"
clinical$educationCat[clinical$education=="short cycle higher education"] <- "short cycle higher education/Technical vocational education"
clinical$educationCat[clinical$education=="Technical vocational education"] <- "short cycle higher education/Technical vocational education"
clinical$educationCat[clinical$education=="Upper secondary education"] <- "Upper secondary education"
clinical$educationCat[clinical$education=="Other"] <- "Primary school or other"
clinical$educationCat[clinical$education=="Primary school"] <- "Primary school or other"
table(clinical$educationCat, useNA="always")/25
prop.table(table(clinical$educationCat, useNA="always"))

## occupation
table(clinical$occupation, useNA="always")/25
prop.table(table(clinical$occupation))

## recategorize occupation
clinical$occupationCat[clinical$occupation=="unemployed"] <- "Other"
clinical$occupationCat[clinical$occupation=="longterm sick leave"] <- "Other"
clinical$occupationCat[clinical$occupation=="Other"] <- "Other"
clinical$occupationCat[clinical$occupation=="employed"] <- "Employed"
clinical$occupationCat[clinical$occupation=="outside labor market"] <- "Other"
clinical$occupationCat[clinical$occupation=="student"] <- "student"
table(clinical$occupationCat, useNA="always")/25
prop.table(table(clinical$occupationCat, useNA="always"))


## night-time smartphone use
table(clinical$mobileUseNight, useNA="always")/25
publish(univariateTable( ~ mobileUseNight,data=clinical, column.percent=TRUE))

##----------------------------------------------------------------------------##
## confounders according to night-time smartphone use
## age
publish(univariateTable(mobileUseNight ~ age,data=clinical, column.percent=TRUE))

## SD for mobileUseNight=Never
sqrt( 
  mean(aggregate(subset(clinical,imputation!=0 & mobileUseNight=="Never")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="Never")$imputation),FUN=var)$x)+sum((aggregate(subset(clinical,imputation!=0 & mobileUseNight=="Never")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="Never")$imputation),FUN=mean)$x-mean(subset(clinical,imputation!=0 & mobileUseNight=="Never")$age))^2)/25
)

## SD formobileUseNight=A few times a month or less
sqrt( 
  mean(aggregate(subset(clinical,imputation!=0 & mobileUseNight=="A few times a month or less")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="A few times a month or less")$imputation),FUN=var)$x)+sum((aggregate(subset(clinical,imputation!=0 & mobileUseNight=="A few times a month or less")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="A few times a month or less")$imputation),FUN=mean)$x-mean(subset(clinical,imputation!=0 & mobileUseNight=="A few times a month or less")$age))^2)/25
)

## SD for mobileUseNight=A few times a week
sqrt( 
  mean(aggregate(subset(clinical,imputation!=0 & mobileUseNight=="A few times a week")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="A few times a week")$imputation),FUN=var)$x)+sum((aggregate(subset(clinical,imputation!=0 & mobileUseNight=="A few times a week")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="A few times a week")$imputation),FUN=mean)$x-mean(subset(clinical,imputation!=0 & mobileUseNight=="A few times a week")$age))^2)/25
)

## SD for mobileUseNight=Every night or almost every night
sqrt( 
  mean(aggregate(subset(clinical,imputation!=0 & mobileUseNight=="Every night or almost every night")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="Every night or almost every night")$imputation),FUN=var)$x)+sum((aggregate(subset(clinical,imputation!=0 & mobileUseNight=="Every night or almost every night")$age,by=list(subset(clinical,imputation!=0 & mobileUseNight=="Every night or almost every night")$imputation),FUN=mean)$x-mean(subset(clinical,imputation!=0 & mobileUseNight=="Every night or almost every night")$age))^2)/25
)

## education
table(clinical$mobileUseNight, clinical$education)/25
publish(univariateTable(mobileUseNight ~ education,data=clinical, column.percent=TRUE))

# recategorized education
table(clinical$mobileUseNight, clinical$educationCat)/25
publish(univariateTable(mobileUseNight ~ educationCat,data=clinical, column.percent=TRUE))


## occupation
table(clinical$mobileUseNight, clinical$occupation)/25
publish(univariateTable(mobileUseNight ~ occupation,data=clinical, column.percent=TRUE))

# recategorized occupation
table(clinical$mobileUseNight, clinical$occupationCat)/25
publish(univariateTable(mobileUseNight ~ occupationCat,data=clinical, column.percent=TRUE))





# --------------------------------------------------------------------------- ##
  
#Reading in the data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
CSS <- read.csv2("Citizen Science Sample/imp_citizenScience.csv")
table(CSS$imputation)
CSS <- subset(CSS,imputation!=0)

#tracking_CSS <- subset(subject_tracking_clusters, followup==1)
#table(subject_tracking_clusters$followup)
  
CSSAlt <- merge(tracking_CSS, CSS, by="userid")

## distribution only in survey data

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


## night-time smartphone use
table(CSS$mobileUseNight, useNA="always")/25
publish(univariateTable( ~ mobileUseNight,data=CSS, column.percent=TRUE))

## bmi
table(CSS$bmi, useNA="always")
CSS$bmiCat[CSS$bmi<25] <- "<25"
CSS$bmiCat[CSS$bmi>=25 & CSS$bmi<30] <- "25-30"
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

