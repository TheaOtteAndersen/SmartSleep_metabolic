#----------------------------------------------------------------------
## authors: Thea Otte Andersen
## created: 17 November 2021
## Version: 
## last-updated: 
##           By: 
##     Update #:
#----------------------------------------------------------------------
## 

### Commentary: Calculation of weights used for the SmartSleep Population Sample (N=4522)

## ========================================================================== ##
## load libraries, import data, --------------------------------------------- ##
library(survey)
library(readxl)

## ========================================================================== ##
## distribution in Denmark
zip <- read_excel("S:/SUND-IFSV-SmartSleep/Thea/weights/Denmark_ZIP_clean.xls")
## ========================================================================== ##

## ========================================================================== ##
## Population sample N=4522
## ========================================================================== ##

## load data set
setwd("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Population Sample/")
conv <- read.table("Population.csv", header=TRUE, fill=TRUE, quote = "", sep=";", stringsAsFactors = TRUE)


#renaming zip file

zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Hovedstaden"] <- "Capital region"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Sjælland"] <- "Region zealand"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Syddanmark"] <- "Region of southern denmark"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Midtjylland"] <- "Central denmark region"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Nordjylland"] <- "North denmark region"
unique(zip$ADRESSERINGSNAVN)

#recoding postalcode (NA = 137)
table(conv$postnr, useNA="always")
length(unique(conv$postnr))
conv <- subset(conv, is.na(conv$postnr)==FALSE)

conv$regionD <- c("")

for (i in (unique(conv$postnr))) {
  if(i %in% zip$POSTNR) {
    conv$regionD[conv$postnr == i] <- zip$ADRESSERINGSNAVN[zip$POSTNR == i]
  } else {
    conv$regionD[conv$postnr == i] <- NA
  }
}

table(conv$regionD, useNA="always")
conv$postnr[is.na(conv$regionD)]

#Age
table(conv$agePNR, useNA="always")

conv$ageCat <- c("")

## Age categories from 18-50 years
conv$ageCat[conv$agePNR < 15 ] <- NA
conv$ageCat[conv$agePNR >= 18 & conv$agePNR <= 22] <- "18-22"
conv$ageCat[conv$agePNR >= 23 & conv$agePNR <= 27] <- "23-27"
conv$ageCat[conv$agePNR >= 28 & conv$agePNR <= 32] <- "28-32"
conv$ageCat[conv$agePNR >= 33 & conv$agePNR <= 37] <- "33-37"
conv$ageCat[conv$agePNR >= 38 & conv$agePNR <= 42] <- "38-42"
conv$ageCat[conv$agePNR >= 43 & conv$agePNR <= 47] <- "43-47"
conv$ageCat[conv$agePNR >= 48 & conv$agePNR <= 51] <- "48-50"
conv$ageCat[conv$ageCat == ""] <- NA
table(conv$ageCat, useNA="always")

#Gender 
table(conv$Gender, useNA="always")

conv$sex <- c("")
conv$sex[conv$Gender== "Female" ] <- "Woman"
conv$sex[conv$Gender == "Male" ] <- "Man"
conv$sex[conv$Gender == "Other" ] <- NA
conv$sex[conv$sex == "" ] <- NA
table(conv$sex, useNA="always")

#create dataset without missing data on age, gender and region. (137 participants are excluded)
conv <- subset(conv, is.na(conv$sex)==FALSE & is.na(conv$ageCat)==FALSE & is.na(conv$regionD)==FALSE)

##Create unweighted survey object
conv.svy.unweighted <- svydesign(ids=~1, data=conv)

#Setting up population dataframes. Data comes from register FOLK1A, 1st quarter 2020 (2020 K3): https://statistikbanken.dk/statbank5a/selectvarval/define.asp?PLanguage=0&MainTable=FOLK1A&TabStrip=Select
## distrubuion on region and gender is restricted to agegroup 18-50

gender.dist <- data.frame(sex = c("Man", "Woman"),
                          Freq = nrow(conv) * c(0.50, 0.50))

age.dist <- data.frame(ageCat = c("18-22", "23-27", "28-32", "33-37", "38-42", "43-47", "48-50"),
                       Freq = nrow(conv) * c(0.15, 0.17, 0.15, 0.14, 0.14, 0.16, 0.09))

region.dist <- data.frame(regionD = c("Capital region", "Region zealand", "Region of southern denmark", "Central denmark region","North denmark region"),
                          Freq = nrow(conv) * c(0.35, 0.13, 0.19, 0.23, 0.10))

#Calculating weights
conv.svy.rake <- rake(design = conv.svy.unweighted,
                      sample.margins = list(~sex, ~ageCat, ~regionD),
                      population.margins = list(gender.dist, age.dist,region.dist))

summary(weights(conv.svy.rake))

#Trim weights if to small ()smaller than 0.3 or too large (larger than 3)
conv.svy.rake.trim <- trimWeights(conv.svy.rake, lower=0.3, upper=3,
                                  strict=TRUE) 
summary(weights(conv.svy.rake.trim))

#Append weights to dataset
convW <- cbind(conv, Weights = weights(conv.svy.rake.trim))

## ========================================================================== ##
## save data (convW)
## ========================================================================== ##
#N=
save(convW, file="S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Population Sample/PopulationWeighted.RData")
write.table(convW, file="S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Population Sample/PopulationWeighted.csv",sep=";",row.names=FALSE, quote = FALSE)
