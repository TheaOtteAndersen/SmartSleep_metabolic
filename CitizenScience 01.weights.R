#----------------------------------------------------------------------
## authors: Thea Otte Andersen
## created: 17 November 2021
## Version: 
## last-updated: 
##           By: 
##     Update #:
#----------------------------------------------------------------------
## 

### Commentary: Calculation of weights used for the SmartSleep Baseline citizen science sample (N=25135) and follow-up (N=1885)

## ========================================================================== ##
## load libraries, import data, --------------------------------------------- ##
library(survey)
library(readxl)

## ========================================================================== ##
## distribution in Denmark
zip <- read_excel("S:/SUND-IFSV-SmartSleep/Thea/weights/Denmark_ZIP_clean.xls")
## ========================================================================== ##

## ========================================================================== ##
## baseline citizen science sample N=25,135
## ========================================================================== ##

## load data set
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/SmartSleep Experimentet/Data/Renset data")
conv <- read.table("SmartSleepExp.csv", header=TRUE, fill=TRUE, quote = "", sep=";", stringsAsFactors = TRUE, strip.white = TRUE)

#renaming zip file
unique(zip$ADRESSERINGSNAVN)
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Hovedstaden"] <- "Capital region"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Sjælland"] <- "Region zealand"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Syddanmark"] <- "Region of southern denmark"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Midtjylland"] <- "Central denmark region"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Nordjylland"] <- "North denmark region"
unique(zip$ADRESSERINGSNAVN)

#recoding postalcode
table(conv$zipCode)
#sum(is.na(conv$zipCode))
length(unique(conv$zipCode))
conv$regionD <- c("")

for (i in (unique(conv$zipCode))) {
  if(i %in% zip$POSTNR) {
    conv$regionD[conv$zipCode == i] <- zip$ADRESSERINGSNAVN[zip$POSTNR == i]
  } else {
    conv$regionD[conv$zipCode == i] <- NA
  }
}

table(conv$regionD, useNA="always")
conv$zipCode[is.na(conv$regionD)]

#Age
table(conv$age, useNA="always")
conv$ageCat <- c("")

## Age categories from 16-99 years
conv$ageCat[conv$age < 15 ] <- NA
conv$ageCat[conv$age >= 15 & conv$age <= 19] <- "15-19"
conv$ageCat[conv$age >= 20 & conv$age <= 24] <- "20-24"
conv$ageCat[conv$age >= 25 & conv$age <= 29] <- "25-29"
conv$ageCat[conv$age >= 30 & conv$age <= 34] <- "30-34"
conv$ageCat[conv$age >= 35 & conv$age <= 39] <- "35-39"
conv$ageCat[conv$age >= 40 & conv$age <= 44] <- "40-44"
conv$ageCat[conv$age >= 45 & conv$age <= 49] <- "45-49"
conv$ageCat[conv$age >= 50 & conv$age <= 54] <- "50-54"
conv$ageCat[conv$age >= 55 & conv$age <= 59] <- "55-59"
conv$ageCat[conv$age >= 60 & conv$age <= 64] <- "60-64"
conv$ageCat[conv$age >= 65 & conv$age <= 69] <- "65-69"
conv$ageCat[conv$age >= 70 & conv$age <= 74] <- "70-74"
conv$ageCat[conv$age >= 75 & conv$age <= 79] <- "75-79"
conv$ageCat[conv$age >=80] <- "80+"
table(conv$ageCat, useNA="always")

#Gender 
conv$sex <- c("")
table(conv$gender, useNA="always")
conv$sex[conv$gender == "Female" ] <- "Woman"
conv$sex[conv$gender == "Male" ] <- "Man"
conv$sex[conv$gender == "Do not want to answer" ] <- NA
conv$sex[conv$gender == "Other" ] <- NA
table(conv$sex, useNA="always")

#create dataset without missing data on age, gender and region. (210 participants are excluded)
conv <- subset(conv, is.na(conv$sex)==FALSE & is.na(conv$ageCat)==FALSE & is.na(conv$regionD)==FALSE)

##Create unweighted survey object
conv.svy.unweighted <- svydesign(ids=~1, data=conv)

#Setting up population dataframes. Data comes from register FOLK1A, 4th quarter 2018 (2018 K4): https://statistikbanken.dk/statbank5a/selectvarval/define.asp?PLanguage=0&MainTable=FOLK1A&TabStrip=Select
## distrubuion on region and gender is for all men and women aged 16-81 years

gender.dist <- data.frame(sex = c("Man", "Woman"),
                          Freq = nrow(conv) * c(0.50, 0.50))

age.dist <- data.frame(ageCat = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64","65-69", "70-74","75-79", "80+"),
                       Freq = nrow(conv) * c(0.06, 0.08, 0.08, 0.07, 0.07, 0.08, 0.08, 0.09, 0.08, 0.07, 0.07, 0.07, 0.05, 0.05))

region.dist <- data.frame(regionD = c("Capital region", "Region zealand", "Region of southern denmark", "Central denmark region","North denmark region"),
                          Freq = nrow(conv) * c(0.32, 0.14, 0.21, 0.23, 0.10))


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

save(convW, file="S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/SmartSleepExpWeighted.RData")
write.table(convW, file="S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/SmartSleepExpWeighted.csv",sep=";",row.names=FALSE, quote = FALSE)

## ========================================================================== ##
## Follow-up citizen science data
## ========================================================================== ##

## ========================================================================== ##
## load libraries, import data, --------------------------------------------- ##
library(survey)
library(readxl)

## load data set
## questionnaire data
setwd("S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/")
conv <- read.table("FollowUp.csv", header=TRUE, fill=TRUE, quote = "", sep=";", stringsAsFactors = TRUE, strip.white = TRUE)

zip <- read_excel("S:/SUND-IFSV-SmartSleep/Thea/weights/Denmark_ZIP_clean.xls")

#renaming zip file
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Hovedstaden"] <- "Capital region"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Sjælland"] <- "Region zealand"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Syddanmark"] <- "Region of southern denmark"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Midtjylland"] <- "Central denmark region"
zip$ADRESSERINGSNAVN[zip$ADRESSERINGSNAVN == "Region Nordjylland"] <- "North denmark region"
unique(zip$ADRESSERINGSNAVN)

#str(conv, list.len=425)

#conv <- subset(conv, select=-c(age.x, gender.x,liveAlone.x, education.x, educationOther.x, occupation.x, occupationOther.x, mobilephone,
                              mobilephoneOther, mobileCheck.x ,pmpu1.x, pmpu2.x, pmpu3.x, pmpu4.x, pmpu5.x, pmpu6.x, pmpu7.x, mobileAddiction.x,
                              soMeFace.x, soMeInsta.x, soMeSnap.x, soMeTwit.x, soMeLinkd.x, soMePint.x, soMeYout, soMeRed.x, soMeWhats.x, soMeMess.x, soMeOther.x, soMeNoUse.x, soMeFrequency.x, 
                              sleepTime.x, riseTime.x, sleepUncalm.x, sleepFallingASleep.x, sleepMorning.x, sleepAwakeNight.x, mobileUseBeforeSleep.x, disturbance.x, mobileUseNight.x, 
                              otherDevicesNight, typeOfDevice, numberMobileCheck.x, numberDeviceCheck, timeMobileUse.x, timeDeviceUse, checkFirstHour.x, checkNight.x, checkMorning.x, checkDeviceFirstHour, 
                              checkDeviceNight, checkDeviceMorning, reasonMobileUseWait.x, reasonMobileUseHabit.x,reasonMobileUseHelp.x, reasonMobileUseFomo.x, reasonMobileUseClock.x, reasonMobileUseOther.x,
                              reasonMobileUseSpecify, reasonAwakeFallAsleep.x, reasonAwakeMyself.x, reasonAwakeNoiseMobile.x, reasonAwakeNoise.x, reasonAwakeCheckMobile.x, reasonAwakeOther.x, reasonAwakeSpecify, 
                              reasonDeviceUseWait, reasonDeviceUseHabit , reasonDeviceUseFallAsleep, reasonDeviceUseFomo, reasonDeviceUseClock, reasonDeviceUseOther , reasonDeviceUseSpecify , reasonDeviceAwakeFallAsleep, 
                              reasonDeviceAwakeNoiseDevice, reasonDeviceAwakeOtherNoise, reasonDeviceAwakeMyself, reasonDeviceAwakeCheck, reasonDeviceAwakeOther, reasonDeviceAwakeSpecify, 
                              platformStreaming.x, platformSoMe.x, platformMusic.x, platformDeviceSoMe, platformDeviceMusic, platformDeviceNews, platformDeviceShopping, 
                              platformDeviceSport, platformDeviceGames, platformDeviceEbook, platformDeviceMail, platformDeviceCalls, platformDeviceAlarm, platformDeviceOther, platformDeviceSpecify, 
                              someAnswers, accomplished.x, dropOut, id, ageD, genderD, regionD, 
                              occupationD.x, educationD.x, educationS, mobileCheckCat, pmpu1Re.x, pmpu1ReNum.x,  pmpu2Re.x, pmpu2ReNum.x, pmpu3Re.x, pmpu3ReNum.x, pmpu4Re.x, pmpu4ReNum.x, pmpu5Re.x, pmpu5ReNum.x, 
                              pmpu6Re.x, pmpu6ReNum.x, pmpu7Re.x, pmpu7ReNum.x, pmpuScale.x, pmpuScaleCat.x, sleepUncalmRe.x , sleepUncalmReNum.x, sleepFallingASleepRe.x, sleepFallingASleepReNum.x, sleepMorningRe.x, 
                              sleepMorningReNum.x, sleepAwakeNightRe.x, sleepAwakeNightReNum.x, sleepQualityScale.x, sleepTimeD, riseTimeD, sleepDur, sleepDurDCat, mobileUseBeforeSleepCat, mobileUseNightCat, heightD, 
                              bmi.x, bmiCat, stressControlRe, stressControlReNum, stressProblemsRe, stressProblemsReNum, stressGoodLifeRe, stressGoodLifeReNum, stressCopingRe, stressCopingReNum, stressScale.x, 
                              platformNews.x, platformShopping.x, platformSport.x, platformGames.x, platformEbook.x, platformMail.x,
                              platformCalls.x, platformAlarm.x, platformOther.x, platformSpecify, platformDeviceStreaming, 
                              locationMobileNight.x, partnerMobileUse.x, precautionNo.x, precautionSilent.x, precautionTurnOff.x, precautionOutsideBedroom, precautionSpecify, health.x, height.x, 
                              weight.x, stress, stressControl, stressProblems, stressGoodLife, stressCoping, headache.x))

#recoding postalcode
table(conv$zipCode, useNA = "always")
length(unique(conv$zipCode))
conv$regionD <- c("")

for (i in (unique(conv$zipCode))) {
  if(i %in% zip$POSTNR) {
    conv$regionD[conv$zipCode == i] <- zip$ADRESSERINGSNAVN[zip$POSTNR == i]
  } else {
    conv$regionD[conv$zipCode == i] <- NA
  }
}

table(conv$regionD, useNA="always")
conv$zipCode[is.na(conv$regionD)]

#Age (16-81)
table(conv$age.y, useNA="always")
conv$ageCat <- c("")

## Age categories from 16-81 years

conv$ageCat[conv$age.y >= 16 & conv$age.y <= 20] <- "16-20"
conv$ageCat[conv$age.y >= 21 & conv$age.y <= 25] <- "21-25"
conv$ageCat[conv$age.y >= 26 & conv$age.y <= 30] <- "26-30"
conv$ageCat[conv$age.y >= 31 & conv$age.y <= 35] <- "31-35"
conv$ageCat[conv$age.y >= 36 & conv$age.y <= 40] <- "36-40"
conv$ageCat[conv$age.y >= 41 & conv$age.y <= 45] <- "41-45"
conv$ageCat[conv$age.y >= 46 & conv$age.y <= 50] <- "46-50"
conv$ageCat[conv$age.y >= 51 & conv$age.y <= 55] <- "51-55"
conv$ageCat[conv$age.y >= 56 & conv$age.y <= 60] <- "56-60"
conv$ageCat[conv$age.y >= 61 & conv$age.y <= 65] <- "61-65"
conv$ageCat[conv$age.y >= 66 ] <- "66+"
conv$ageCat[conv$ageCat == "" ] <- NA
table(conv$ageCat, useNA="always")

#Gender 
conv$sex <- c("")
table(conv$gender.y, useNA="always")
conv$sex[conv$gender.y == "Female" ] <- "Woman"
conv$sex[conv$gender.y == "Male" ] <- "Man"
conv$sex[conv$gender.y == "Other" ] <- NA
conv$sex[conv$sex == "" ] <- NA
table(conv$sex, useNA="always")

#create data set without missing data on age, gender and region. (53 participants are excluded)
conv <- subset(conv, is.na(conv$sex)==FALSE & is.na(conv$ageCat)==FALSE & is.na(conv$regionD)==FALSE)

##Create unweighted survey object
conv.svy.unweighted <- svydesign(ids=~1, data=conv)

#Setting up population dataframes. Data comes from register FOLK1A, 1st quarter 2020 (2020 K1): https://statistikbanken.dk/statbank5a/selectvarval/define.asp?PLanguage=0&MainTable=FOLK1A&TabStrip=Select
## distrubuion on region and gender is for all men and women aged 16-81 years

gender.dist <- data.frame(sex = c("Man", "Woman"),
                          Freq = nrow(conv) * c(0.50, 0.50))

age.dist <- data.frame(ageCat = c("16-20", "21-25", "26-30", "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", "61-65","66+"),
                       Freq = nrow(conv) * c(0.08, 0.09, 0.09, 0.08, 0.07, 0.08, 0.08, 0.09, 0.08, 0.07, 0.19))

region.dist <- data.frame(regionD = c("Capital region", "Region zealand", "Region of southern denmark", "Central denmark region","North denmark region"),
                          Freq = nrow(conv) * c(0.32, 0.14, 0.21, 0.23, 0.10))


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
#N=1832
save(convW, file="S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/FollowUpWeighted.RData")
write.table(convW, file="S:/SUND-IFSV-SmartSleep/Thea/Clusters, obesity and metabolic biomarkers/Data/Citizen Science Sample/FollowUpWeighted.csv",sep=";",row.names=FALSE, quote = FALSE)


