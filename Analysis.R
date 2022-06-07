
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
library(ggtern)

expit = function(x) exp(x)/(1+exp(x))

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
subject_tracking_six_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_clusters.csv")
subject_tracking_four_clusters <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Tracking data/subject_tracking_four_clusters.csv")

## Mapping compositions to R^(D-1)
SBP6 <- matrix(c(1,-1,-1,-1,-1,-1,0,-1,-1,1,1,1,0,-1,1,0,0,0,0,0,0,1,1,-1,0,0,0,1,-1,0),byrow=T,nrow=5,ncol=6)
## Coefficients are: Intercept, No activity vs rest (1 vs 2,3,4,5,6), little activity vs more activity, little offset vs little onset, moderate activity vs much activity, moderate offset vs moderate onset
SBP4 <- matrix(c(1,1,-1,1,1,-1,0,1,1,0,0,-1),byrow=T,nrow=3,ncol=4)
## Coefficients are: Intercept, No activity vs rest (1 vs 2,3,4), little activity vs more activity, offset with activity vs onset with activity


X6 <- subject_tracking_six_clusters[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob")]
X4 <- subject_tracking_four_clusters[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob")]

Phi.generator <- function(X) { #Orthonormal basis
  K <- nrow(X)
  L <- ncol(X)
  A <- matrix(rep(NA,K*L),nrow = K, ncol = L)
  for (i in 1:K){
    r <- sum(X[i,]>0)
    s <- sum(X[i,]<0)
    A[i,] <- X[i,]*sqrt((r*s)/(r+s))*((X[i,]>0)*1/r + (X[i,]<0)*1/s)
  }
  return(A)
}

Orthogonal.coordinate.generator <- function(X,l) { #Inputs data and base index, outputs coordinates in a system with orthogonal basis elements
  D <- ncol(X)
  M <- nrow(X)
  A <- matrix(rep(NA,(D-1)*M),ncol = D-1, nrow = M)
  for (m in 1:M){
  for (i in 1:(D-1)){
    A[m,i] <- log2(unlist(c(X[m,l],X[m,-l]))/(prod(unlist(c(X[m,l],X[m,-l]))[(i+1):D])^(1/(D-i))))[i]
  }
  }
  return(A)
}

Phi6 <- Phi.generator(SBP6)
Phi4 <- Phi.generator(SBP4)

ilrX6 <- as.matrix(log(X6)) %*% t(Phi6)
ilrX4 <- as.matrix(log(X4)) %*% t(Phi4)

clrX6 <- ilrX6 %*% Phi6
clrX4 <- ilrX4 %*% Phi4

## The ilr coordinates can be readily used as covariates in regression models. They are coefficients with respect to a basis of the simplex space (in this case the SBP basis).
## The clr coordinates should not be used as covariates in regression models as they will lead to singular covariance matrices, as they are constrained to have component sums zero and thus are not coefficients with respect to a basis of the simplex space.

## We may afterwards back-transform fitted coefficients to coefficients in the original simplex space, thereby getting coefficients associated with each component.

subject_tracking_six_clusters$ilr1 <- ilrX6[,1]
subject_tracking_six_clusters$ilr2 <- ilrX6[,2]
subject_tracking_six_clusters$ilr3 <- ilrX6[,3]
subject_tracking_six_clusters$ilr4 <- ilrX6[,4]
subject_tracking_six_clusters$ilr5 <- ilrX6[,5]

subject_tracking_four_clusters$ilr1 <- ilrX4[,1]
subject_tracking_four_clusters$ilr2 <- ilrX4[,2]
subject_tracking_four_clusters$ilr3 <- ilrX4[,3]

## Alternative basis 

ilrX6.orthogonal1 <- Orthogonal.coordinate.generator(X6,1)
ilrX6.orthogonal2 <- Orthogonal.coordinate.generator(X6,2)
ilrX6.orthogonal3 <- Orthogonal.coordinate.generator(X6,3)
ilrX6.orthogonal4 <- Orthogonal.coordinate.generator(X6,4)
ilrX6.orthogonal5 <- Orthogonal.coordinate.generator(X6,5)
ilrX6.orthogonal6 <- Orthogonal.coordinate.generator(X6,6)

ilrX4.orthogonal1 <- Orthogonal.coordinate.generator(X4,1)
ilrX4.orthogonal2 <- Orthogonal.coordinate.generator(X4,2)
ilrX4.orthogonal3 <- Orthogonal.coordinate.generator(X4,3)
ilrX4.orthogonal4 <- Orthogonal.coordinate.generator(X4,4)

## We can use this to make 6 and 4 regressions respectively, taking one poi from each.
## Seems to be a valid approach, but is cumbersome.

## Collecting the two clusterings in one file

subject_tracking_clusters <- left_join(subject_tracking_six_clusters,subject_tracking_four_clusters[,c("userid","cluster","cluster1prob","cluster2prob","cluster3prob","cluster4prob","description","state0prob","state1prob","state2prob","state3prob","ilr1","ilr2","ilr3")],by="userid")
subject_tracking_clusters <- rename(subject_tracking_clusters,ilr1 = ilr1.x, ilr2=ilr2.x, ilr3=ilr3.x,cluster1prob=cluster1prob.x,cluster2prob=cluster2prob.x,cluster3prob=cluster3prob.x,cluster4prob=cluster4prob.x,
                                    state0prob=state0prob.x,state1prob=state1prob.x,state2prob=state2prob.x,state3prob=state3prob.x,cluster=cluster.x,description=description.x)

subject_tracking_clusters$ilrX6.orthogonal1.1 <- ilrX6.orthogonal1[,1]
subject_tracking_clusters$ilrX6.orthogonal1.2 <- ilrX6.orthogonal1[,2]
subject_tracking_clusters$ilrX6.orthogonal1.3 <- ilrX6.orthogonal1[,3]
subject_tracking_clusters$ilrX6.orthogonal1.4 <- ilrX6.orthogonal1[,4]
subject_tracking_clusters$ilrX6.orthogonal1.5 <- ilrX6.orthogonal1[,5]
                                                                   
subject_tracking_clusters$ilrX6.orthogonal2.1 <- ilrX6.orthogonal2[,1]
subject_tracking_clusters$ilrX6.orthogonal2.2 <- ilrX6.orthogonal2[,2]
subject_tracking_clusters$ilrX6.orthogonal2.3 <- ilrX6.orthogonal2[,3]
subject_tracking_clusters$ilrX6.orthogonal2.4 <- ilrX6.orthogonal2[,4]
subject_tracking_clusters$ilrX6.orthogonal2.5 <- ilrX6.orthogonal2[,5]

subject_tracking_clusters$ilrX6.orthogonal3.1 <- ilrX6.orthogonal3[,1]
subject_tracking_clusters$ilrX6.orthogonal3.2 <- ilrX6.orthogonal3[,2]
subject_tracking_clusters$ilrX6.orthogonal3.3 <- ilrX6.orthogonal3[,3]
subject_tracking_clusters$ilrX6.orthogonal3.4 <- ilrX6.orthogonal3[,4]
subject_tracking_clusters$ilrX6.orthogonal3.5 <- ilrX6.orthogonal3[,5]

subject_tracking_clusters$ilrX6.orthogonal4.1 <- ilrX6.orthogonal4[,1]
subject_tracking_clusters$ilrX6.orthogonal4.2 <- ilrX6.orthogonal4[,2]
subject_tracking_clusters$ilrX6.orthogonal4.3 <- ilrX6.orthogonal4[,3]
subject_tracking_clusters$ilrX6.orthogonal4.4 <- ilrX6.orthogonal4[,4]
subject_tracking_clusters$ilrX6.orthogonal4.5 <- ilrX6.orthogonal4[,5]

subject_tracking_clusters$ilrX6.orthogonal5.1 <- ilrX6.orthogonal5[,1]
subject_tracking_clusters$ilrX6.orthogonal5.2 <- ilrX6.orthogonal5[,2]
subject_tracking_clusters$ilrX6.orthogonal5.3 <- ilrX6.orthogonal5[,3]
subject_tracking_clusters$ilrX6.orthogonal5.4 <- ilrX6.orthogonal5[,4]
subject_tracking_clusters$ilrX6.orthogonal5.5 <- ilrX6.orthogonal5[,5]

subject_tracking_clusters$ilrX6.orthogonal6.1 <- ilrX6.orthogonal6[,1]
subject_tracking_clusters$ilrX6.orthogonal6.2 <- ilrX6.orthogonal6[,2]
subject_tracking_clusters$ilrX6.orthogonal6.3 <- ilrX6.orthogonal6[,3]
subject_tracking_clusters$ilrX6.orthogonal6.4 <- ilrX6.orthogonal6[,4]
subject_tracking_clusters$ilrX6.orthogonal6.5 <- ilrX6.orthogonal6[,5]

##

subject_tracking_clusters$ilrX4.orthogonal1.1 <- ilrX4.orthogonal1[,1]
subject_tracking_clusters$ilrX4.orthogonal1.2 <- ilrX4.orthogonal1[,2]
subject_tracking_clusters$ilrX4.orthogonal1.3 <- ilrX4.orthogonal1[,3]

subject_tracking_clusters$ilrX4.orthogonal2.1 <- ilrX4.orthogonal2[,1]
subject_tracking_clusters$ilrX4.orthogonal2.2 <- ilrX4.orthogonal2[,2]
subject_tracking_clusters$ilrX4.orthogonal2.3 <- ilrX4.orthogonal2[,3]

subject_tracking_clusters$ilrX4.orthogonal3.1 <- ilrX4.orthogonal3[,1]
subject_tracking_clusters$ilrX4.orthogonal3.2 <- ilrX4.orthogonal3[,2]
subject_tracking_clusters$ilrX4.orthogonal3.3 <- ilrX4.orthogonal3[,3]

subject_tracking_clusters$ilrX4.orthogonal4.1 <- ilrX4.orthogonal4[,1]
subject_tracking_clusters$ilrX4.orthogonal4.2 <- ilrX4.orthogonal4[,2]
subject_tracking_clusters$ilrX4.orthogonal4.3 <- ilrX4.orthogonal4[,3]

## load baseline data
base_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv")
base_data$mobileUseNight <- factor(base_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
base_data$mobileUseBeforeSleep <- factor(base_data$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))

## load followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")
CSS$mobileUseNight <- factor(CSS$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
CSS$mobileUseBeforeSleep <- factor(CSS$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))

## load population sample
pop_data <-read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv")
pop_data$mobileUseNight <- factor(pop_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
pop_data$mobileUseBeforeSleep <- factor(pop_data$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))

## load clinical data (survey and clinical data)
clin_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample/imp_clinical.csv")
clin_clinical <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/SmartSleep Clinical/Data/Rådata/SmartSleepClinicalData.csv")

clin_data$mobileUseNight <- factor(clin_data$mobileUseNight, levels = c("Never","A few times a month or less","A few times a week","Every night or almost every night"))
clin_data$mobileUseBeforeSleep <- factor(clin_data$mobileUseBeforeSleep, levels = c("Never","Every month or less","Once a week","2-4 times per week","5-7 times per week"))


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
table(base_data$selfScoreCat[base_data$imputation!=0],useNA="always")/25

# bar chart selfScoreCat 
ggplot(base_data, aes(x = factor(selfScoreCat))) +
  geom_bar()

## tjek risk profiles
publish(univariateTable( ~ mobileUseBeforeSleep,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))
publish(univariateTable( ~ mobileCheck,data=base_data, column.percent=TRUE))
table((base_data$pmpuScale<14)*1+(base_data$pmpuScale>=14 & base_data$pmpuScale<17)*2+(base_data$pmpuScale>=17 & base_data$pmpuScale<19)*3+(base_data$pmpuScale>=19)*4)/20

publish(univariateTable( ~ selfScoreCat,data=base_data, column.percent=TRUE))
      
base_data$bmi30 <- (base_data$bmi>=30)
base_data$bmi25 <- (base_data$bmi>=25)

#save(base_data,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/base_data.RData")

base_data_mids <- as.mids(base_data,.imp="imputation")

## BMI kategoriseringer ved baseline
table(base_data$bmi<25, base_data$selfScoreCat)/26
table(base_data$bmi>=25&base_data$bmi<30, base_data$selfScoreCat)/26
table(base_data$bmi>=30, base_data$selfScoreCat)/26

table(base_data$bmi<25, base_data$selfScoreCat)/21
table(base_data$bmi>=25&base_data$bmi<30, base_data$selfScoreCat)/21
table(base_data$bmi>=30, base_data$selfScoreCat)/21


# --------------------------------------------------------------------------- ##
#Followup sample - using quartile levels from baseline sample

## risk profiles for CSS
CSS$selfScore <- (CSS$mobileUseBeforeSleep=="5-7 times per week")*4+(CSS$mobileUseBeforeSleep=="2-4 times per week")*3+(CSS$mobileUseBeforeSleep=="Once a week")*3+(CSS$mobileUseBeforeSleep=="Every month or less")*2+(CSS$mobileUseBeforeSleep=="Never")*1+
  (CSS$mobileUseNight=="Every night or almost every night")*4+(CSS$mobileUseNight=="A few times a week")*3+(CSS$mobileUseNight=="A few times a month or less")*2+(CSS$mobileUseNight=="Never")*1+
  (CSS$mobileCheck==">20 times an hour")*4+(CSS$mobileCheck=="11-20 times an hour")*4+(CSS$mobileCheck=="5-10 times an hour")*3+(CSS$mobileCheck=="1-4 times an hour")*2+(CSS$mobileCheck=="Every 2nd hour")*2+(CSS$mobileCheck=="Several times a day")*1+(CSS$mobileCheck=="Once a day or less")*1+
  (CSS$pmpuScale<=14)*1+(CSS$pmpuScale>14 & CSS$pmpuScale<17)*2+(CSS$pmpuScale>=17 & CSS$pmpuScale<19)*3+(CSS$pmpuScale>=19)*4
summary(CSS$selfScore[CSS$imputation!=0])
CSS$selfScoreCat <- NA
CSS$selfScoreCat[!is.na(CSS$selfScore)]<-"1"
CSS$selfScoreCat[CSS$selfScore>=8]="2"
CSS$selfScoreCat[CSS$selfScore>=10]="3"
CSS$selfScoreCat[CSS$selfScore>=12]="4"
table(CSS$selfScoreCat[CSS$imputation!=0], useNA="always")

# bar chart selfScoreCat 
ggplot(CSS, aes(x = factor(selfScoreCat))) +
  geom_bar()


## merge survey and tracking data 
CSS_track <- inner_join(CSS,subject_tracking_clusters,by="userid")
CSS_track$sample_weights <- as.numeric(CSS_track$sample_weights)

#save(CSS_track,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/CSS_track.RData")

CSS_track_mids<-as.mids(CSS_track,.imp="imputation",.id="userid")

# --------------------------------------------------------------------------- ##
#Merging base and followup


#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup

bmi_followup <- rename(inner_join(CSS,base_data,by=c("CS_ID","imputation")),bmi.base=bmi.y,bmi.fu=bmi.x)

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
bmi_followup$followup_time <- (as.Date(str_c(substr(bmi_followup$EndDate,7,10),"-",substr(bmi_followup$EndDate,4,5),"-",substr(bmi_followup$EndDate,1,2)))-as.Date(str_c(substr(bmi_followup$responseDate.y,7,10),"-",substr(bmi_followup$responseDate.y,4,5),"-",substr(bmi_followup$responseDate.y,1,2))))/365.25

bmi_followup_mids <- as.mids(bmi_followup,.imp="imputation")

#New idea: Try to make long format where followup and baseline are at different time points, and then make an interaction effect with time with bmi (indicators) as response.

long_data <- data.frame("bmi"=c(bmi_followup$bmi.base,bmi_followup$bmi.fu),"userid"=rep(bmi_followup$userid,2),"sample_weights"=rep(bmi_followup$sample_weights.x,2),"gender"=rep(bmi_followup$gender.y,2),"age"=rep(bmi_followup$age.y,2),
                        education=rep(bmi_followup$education.y,2),occupation=rep(bmi_followup$occupation.y,2),selfScoreCat = rep(bmi_followup$selfScoreCat.y,2),"followup"=c(rep(0,length(bmi_followup$bmi.base)),rep(1,length(bmi_followup$bmi.fu))),
                        "imputation"=rep(bmi_followup$imputation,2),"time" = c(rep(0,length(bmi_followup$bmi.base)),bmi_followup$followup_time))

#save(long_data,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/long_bmi.RData")

long_data_mids <- as.mids(long_data,.imp="imputation")



# --------------------------------------------------------------------------- ##
#Population sample

## risk profiles for population sample
pop_data$selfScore <- (pop_data$mobileUseBeforeSleep=="5-7 times per week")*4+(pop_data$mobileUseBeforeSleep=="2-4 times per week")*3+(pop_data$mobileUseBeforeSleep=="Once a week")*3+(pop_data$mobileUseBeforeSleep=="Every month or less")*2+(pop_data$mobileUseBeforeSleep=="Never")*1+
  (pop_data$mobileUseNight=="Every night or almost every night")*4+(pop_data$mobileUseNight=="A few times a week")*3+(pop_data$mobileUseNight=="A few times a month or less")*2+(pop_data$mobileUseNight=="Never")*1+
  (pop_data$mobileCheck==">20 times an hour")*4+(pop_data$mobileCheck=="11-20 times an hour")*4+(pop_data$mobileCheck=="5-10 times an hour")*3+(pop_data$mobileCheck=="1-4 times an hour")*2+(pop_data$mobileCheck=="Every 2nd hour")*2+(pop_data$mobileCheck=="Several times a day")*1+(pop_data$mobileCheck=="Once a day or less")*1+
  (pop_data$pmpuScale<=14)*1+(pop_data$pmpuScale>14 & pop_data$pmpuScale<17)*2+(pop_data$pmpuScale>=17 & pop_data$pmpuScale<19)*3+(pop_data$pmpuScale>=19)*4
summary(pop_data$selfScore[pop_data$imputation!=0])

## categorise selfScoreCat?
pop_data$selfScoreCat <- NA
pop_data$selfScoreCat[!is.na(pop_data$selfScore)]<-"1"
pop_data$selfScoreCat[pop_data$selfScore>=8]="2"
pop_data$selfScoreCat[pop_data$selfScore>=10]="3"
pop_data$selfScoreCat[pop_data$selfScore>=12]="4"
table(pop_data$selfScoreCat[pop_data$imputation!=0])

## merge tracking and survey data for population sample
pop_track <- inner_join(pop_data,subject_tracking_clusters,by="userid")
pop_track$sample_weights<-as.numeric(pop_track$sample_weights)


## Tracking clusters som én numerisk variabel
pop_track$track_severity <- (pop_track$cluster %in% c("Cluster 1"))*1+(pop_track$cluster %in% c("Cluster 2","Cluster 3"))*2+(pop_track$cluster %in% c("Cluster 5","Cluster 6"))*3+(pop_track$cluster %in% c("Cluster 4"))*4


#save(pop_track,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/pop_track.RData")

pop_track_mids<-as.mids(pop_track,.imp="imputation",.id="userid")

## Clinical Data

## merge survey and clinical data
clinical_sample <- rename(inner_join(clin_data,rename(clin_clinical,PNR=cpr),by="PNR"),bmi.self=bmi.x , bmi.clinical=bmi.y)
## merge with tracking data
clinical_sample <- inner_join(clinical_sample,subject_tracking_clusters,by="userid")



#### -------------------------------- ####

#### Write out data files for gamlss bootstrap

#### -------------------------------- ####


# --------------------------------------------------------------------------- #
boot_path <- "S:/SUND-IFSV-SmartSleep/Christoffer/gamlssBootstrap/"
N_imp <- 25
M <- 10

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

## test kontinuert bmi
plot(fitted(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))),residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))))
hist(residuals(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))),xlim=c(-20,20),breaks=200)
hist(simulate(lm(bmi~(mobileUseNight+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1)))$sim_1,breaks=40) #The bell-shape is not that well suited

#Conclusion: The lm is not by itself appropriate for describing the distribution of BMI.


#Alternative: Pretty good fit. A general family of models.


#Confidence intervals and estimates for continuous BMI:

#MobileUseNight:

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ mobileUseNight+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseNight","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG) #May use BCS instead of BCCG which corresponds to using a t distribution instead of normal. This can fit heavier tails, though in this case a very large df is fitted, meaning that there is not much difference.
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_base <- miceadds::pool_mi(qhat = coefs, u = vcovs)
#pool_inf_base$qbar
#pool_inf_base$ubar
#pool_inf_base$ba
#pool_inf_base$pval


#Seems that we can get stable contrats of the mean (taking in varying medians), in spite of skewness.

#One slightly hacky way to achieve this may be to take one of the fitted models created by fit() and replace the stored coefficients with the final pooled estimates. I haven't done detailed testing but it seems to be working on this simple example:
m$mu.coefficients <- pool_inf_base$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_base$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]


#Confidence intervals:
lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,5],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,5],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,5],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

## estimates og 95%CI for mobileUseNight og continous BMI
confints_base_Night <- cbind(c(lowerCat2,lowerCat3,lowerCat4),c(estCat2,estCat3,estCat4),c(upperCat2,upperCat3,upperCat4))


# MobileUseBeforeSleep

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ mobileUseBeforeSleep+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseBeforeSleep","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG) #May use BCS instead of BCCG which corresponds to using a t distribution instead of normal. This can fit heavier tails, though in this case a very large df is fitted, meaning that there is not much difference.
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_base <- miceadds::pool_mi(qhat = coefs, u = vcovs)
#pool_inf_base$qbar
#pool_inf_base$ubar
#pool_inf_base$ba
#pool_inf_base$pval


#Seems that we can get stable contrats of the mean (taking in varying medians), in spite of skewness.

#One slightly hacky way to achieve this may be to take one of the fitted models created by fit() and replace the stored coefficients with the final pooled estimates. I haven't done detailed testing but it seems to be working on this simple example:
m$mu.coefficients <- pool_inf_base$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_base$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]


#Confidence intervals:
lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,5],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,5],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,5],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCat5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[5,5]+1,sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value - integrate(function(y) y*dBCCG(x=y,mu=1,sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCat5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[5,1],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCat5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[5,6],sigma=exp(pool_inf_base$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

## 95%CI for mobileUseBefore sleep og continuos BMI
confints_base_Before <- cbind(c(lowerCat2,lowerCat3,lowerCat4,lowerCat5),c(estCat2,estCat3,estCat4,estCat5),c(upperCat2,upperCat3,upperCat4,upperCat5))



#Profile intervals
#mod <- quote(gamlss())
#prof.term(model=m,criterion="GD",min=-5,max=5,step=1,plot=T)$CI


##

# sigma = exp(pool_inf_base$qbar[20])

# nu = pool_inf_base$qbar[21]

# mu = summary(pool_inf_base)[2,1]

#Density integration: integrate(function(y) (1/(sqrt(2*pi)*sigma))*(y^(nu-1)/mu^nu)*exp(-(((y/mu)^(nu)-1)/(nu*sigma))^2/2),0,Inf)$value

#Untruncated integration: integrate(function(y) y*(1/(sqrt(2*pi)*sigma))*(y^(nu-1)/mu^nu)*exp(-(((y/mu)^(nu)-1)/(nu*sigma))^2/2),0,Inf)$value


##


#trend

#Because the skewness is constant in covariates, the shift between median and mean is constant in covariates.
#We can thus simply make a test for trend on the median to get a p value, as median equality <=> mean equality.

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","mobileUseNight","mobileUseBeforeSleep","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_baseTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_baseTrend$qbar
pool_inf_baseTrend$ubar
pool_inf_baseTrend$pval

summary(pool_inf_baseTrend)

#Going to means:

m$mu.coefficients <- pool_inf_baseTrend$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]

#interval
lowerCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[2,5],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[2,1],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCatTrendNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[2,6],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

lowerCatTrendBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[3,5],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCatTrendBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[3,1],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCatTrendBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[3,6],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 


confints_baseTrend <-rbind(c(lowerCatTrendNight,estCatTrendNight,upperCatTrendNight),
                           c(lowerCatTrendBS,estCatTrendBS,upperCatTrendBS))

summary(pool_inf_baseTrend)$p[2]
summary(pool_inf_baseTrend)$p[3]


# --------------------------------------------------------------------------- ##
## cross-sectional associations between risk profiles and bmi (25, 30)
# --------------------------------------------------------------------------- ##

## Using the mice package with mids objects
mod30 <- with(base_data_mids,glm(bmi30~(mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod25 <- with(base_data_mids,glm(bmi25~(mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))

pool(mod30)
pool(mod25)

## test for trend
TEST <- with(base_data_mids,glm((bmi>=30)~((as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep))+age+gender+education+occupation), weights=sample_weights,family=binomial))
testT <- summary(pool(TEST), conf.int = T)

TEST2 <- with(base_data_mids,glm((bmi>=25)~((as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep))+age+gender+education+occupation), weights=sample_weights,family=binomial))
test2 <- summary(pool(TEST2), conf.int = T)

## OR for BMI>30
model30 <- summary(pool(mod30),conf.int = T)
cbind(exp(model30$estimate),
exp(model30$`2.5 %`),
exp(model30$`97.5 %`))

## OR for BMI >25
model25 <- summary(pool(mod25), conf.int=T)
cbind(exp(model25$estimate),
exp(model25$`2.5 %`),
exp(model25$`97.5 %`))


# --------------------------------------------------------------------------- ##
## longitudinal analysis of risk scores of smartphone behavior and changes in BMI
# --------------------------------------------------------------------------- ##
#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup


## ----- ##
#Final change analyses 
## ----- ##

## from below 25 to above 25
model25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
model_summary25<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=25 ~ (mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))), conf.int = T)
exp(cbind(model_summary25$estimate[1:4],model_summary25$`2.5 %`[1:4],model_summary25$`97.5 %`[1:4]))


## test for trend

test25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (as.numeric(mobileUseNight.y)+as.numeric(mobileUseBeforeSleep.y)+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
testT25 <- summary(pool(test25), conf.int = T)
#anova(test25,model25)


## from below 30 to above 30
model30 <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
model_summary30<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=30 ~ (mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))), conf.int = T)
exp(cbind(model_summary30$estimate[1:4],model_summary30$`2.5 %`[1:4],model_summary30$`97.5 %`[1:4]))

test30 <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (as.numeric(mobileUseNight.y)+as.numeric(mobileUseBeforeSleep.y)+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
testT30 <- summary(pool(test30), conf.int=T)
#anova(test30,model30)



## Modelling numeric difference in bmi between baseline and followup
hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))

plot(fitted(lm(difference~(mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),
     residuals(lm(difference~(mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))))
hist(residuals(lm(difference~(mobileUseNight.y+mobileUseBeforeSleep.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),breaks=50)


#m <- lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*followup_time-selfScoreCat.y-age.y-gender.y-education.y-occupation.y,weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","selfScoreCat.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights")]))
m <- lm(bmi.fu~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","mobileUseNight.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_diff_Night <- summary(pool(with(bmi_followup_mids,lm(bmi.fu~((mobileUseNight.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))), conf.int = T)
#Does it make sense to include tracking information at follow up, or to include self score at followup as well?
#(reduce noise in measurement?)

cbind(model_summary_diff_Night$estimate[17:20],model_summary_diff_Night$`2.5 %`[17:20],model_summary_diff_Night$`97.5 %`[17:20])

## test for trend (mobileUseNight)
test_num <- with(bmi_followup_mids,lm(bmi.fu~((as.numeric(mobileUseNight.y)+as.numeric(mobileUseBeforeSleep.y)):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))
test_Tnum <- summary(pool(test_num), conf.int=T)


# BeforeSleep
m <- lm(bmi.fu~((mobileUseBeforeSleep.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","mobileUseBeforeSleep.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_diff_Before <- summary(pool(with(bmi_followup_mids,lm(bmi.fu~((mobileUseBeforeSleep.y):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))), conf.int = T)
#Does it make sense to include tracking information at follow up, or to include self score at followup as well?
#(reduce noise in measurement?)

cbind(model_summary_diff_Before$estimate[17:21],model_summary_diff_Before$`2.5 %`[17:21],model_summary_diff_Before$`97.5 %`[17:21])

## test for trend (mobileUseBeforeSleep)
test_num <- with(bmi_followup_mids,lm(bmi.fu~((as.numeric(mobileUseNight.y)+as.numeric(mobileUseBeforeSleep.y)):as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))
test_Tnum <- summary(pool(test_num), conf.int=T)


#Generally:
#Wald intervals with Robust=T are better for misspecified models.
#Profile likelihood intervals are better for models that are close to correct (mostly so for smaller sample sizes).



# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#####Tracking data for the followup CSS sample: 

#Mice-based inference for models: 

## BMI indicators

## Tracking

## >=25

## Interpretations are as in Hron et. al (2016), coefficients for effect of doubling proportion with respect to average contribution of other parts (for the remaining mass)
m25_6clustprob <- pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial)))

beta1 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta2 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta3 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta4 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta5 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta6 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]

six.cluster.betas <- rbind(beta1,beta2,beta3,beta4,beta5,beta6)
six.cluster.betas$term <- c("cluster 1","cluster 2", "cluster 3", "cluster 4", "cluster 5", "cluster 6")

#MSE25_6clustprob <- mean((expit(as.numeric(m25_6clustprob$pooled$estimate %*% t(data.frame(1,CSS_track$ilrX6.orthogonal1.1[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.2[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.3[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.4[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.5[CSS_track$imputation!=0],
#                                                                                           CSS_track$age[CSS_track$imputation!=0],CSS_track$gender[CSS_track$imputation!=0]=="Male",CSS_track$gender[CSS_track$imputation!=0]=="Other",CSS_track$education[CSS_track$imputation!=0] == "Primary school",CSS_track$education[CSS_track$imputation!=0] == "medium cycle higher education",
#                                                                                           CSS_track$education[CSS_track$imputation!=0] == "Other",CSS_track$education[CSS_track$imputation!=0] == "Primary school",CSS_track$education[CSS_track$imputation!=0] == "short cycle higher education",CSS_track$education[CSS_track$imputation!=0] == "Technical vocational education",
#                                                                                           CSS_track$education[CSS_track$imputation!=0] == "Upper secondary education", CSS_track$occupation[CSS_track$imputation!=0] == "Other",CSS_track$occupation[CSS_track$imputation!=0] == "outside labor market",CSS_track$occupation[CSS_track$imputation!=0] == "student",CSS_track$occupation[CSS_track$imputation!=0] == "unemployed")))
#)-(CSS_track$bmi[CSS_track$imputation!=0]>=25))^2)

#note correspondence between parametrizations
m <- glm((bmi>=25) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial,data=CSS_track[CSS_track$imputation==1,])
m$coefficients <- m25_6clustprob$pooled$estimate

MSE25_6clustprob <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=25))^2)

##

m25_4clustprob <- pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial)))

beta1 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta2 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta3 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta4 <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]

four.cluster.betas <- rbind(beta1,beta2,beta3,beta4)
four.cluster.betas$term <- c("cluster 1","cluster 2", "cluster 3", "cluster 4")

# Prediction

m <- glm((bmi>=25) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial,data=CSS_track[CSS_track$imputation==1,])
m$coefficients <- m25_4clustprob$pooled$estimate

MSE25_4clustprob <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=25))^2)

## Hence, six clusters seem to distinguish better in this case.

## Using maximal cluster assignments:
## 6 clusters
m25_6clustmax <- pool(with(CSS_track_mids,glm((bmi>=25) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial)))
summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:6,]

#MSE25_6clustmax <- mean((expit(as.numeric(m25_6clustmax$pooled$estimate %*% t(data.frame(1,CSS_track$cluster[CSS_track$imputation!=0]=="Cluster 2",CSS_track$cluster[CSS_track$imputation!=0]=="Cluster 3",CSS_track$cluster[CSS_track$imputation!=0]=="Cluster 4",CSS_track$cluster[CSS_track$imputation!=0]=="Cluster 5",CSS_track$cluster[CSS_track$imputation!=0]=="Cluster 6",
#                                                                                         CSS_track$age[CSS_track$imputation!=0],CSS_track$gender[CSS_track$imputation!=0]=="Male",CSS_track$gender[CSS_track$imputation!=0]=="Other",CSS_track$education[CSS_track$imputation!=0] == "Primary school",CSS_track$education[CSS_track$imputation!=0] == "medium cycle higher education",
#                                                                                         CSS_track$education[CSS_track$imputation!=0] == "Other",CSS_track$education[CSS_track$imputation!=0] == "Primary school",CSS_track$education[CSS_track$imputation!=0] == "short cycle higher education",CSS_track$education[CSS_track$imputation!=0] == "Technical vocational education",
#                                                                                         CSS_track$education[CSS_track$imputation!=0] == "Upper secondary education", CSS_track$occupation[CSS_track$imputation!=0] == "Other",CSS_track$occupation[CSS_track$imputation!=0] == "outside labor market",CSS_track$occupation[CSS_track$imputation!=0] == "student",CSS_track$occupation[CSS_track$imputation!=0] == "unemployed")))
#)-(CSS_track$bmi[CSS_track$imputation!=0]>=25))^2)

m <- glm((bmi>=25) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial,data=CSS_track[CSS_track$imputation==1,])
m$coefficients <- m25_6clustmax$pooled$estimate

MSE25_6clustmax <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=25))^2)


## 4 clusters
m25_4clustmax <- pool(with(CSS_track_mids,glm((bmi>=25) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial)))
summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:6,]

# Prediction

m <- glm((bmi>=25) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial,data=CSS_track[CSS_track$imputation==1,])
m$coefficients <- m25_4clustmax$pooled$estimate

MSE25_4clustmax <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=25))^2)



## >=30
m30_6clustprob <- pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial)))

beta_30_1 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_2 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_3 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_4 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_5 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_6 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]

six.cluster.beta_30 <- rbind(beta_30_1,beta_30_2,beta_30_3,beta_30_4,beta_30_5,beta_30_6)
six.cluster.beta_30$term <- c("cluster 1","cluster 2", "cluster 3", "cluster 4", "cluster 5", "cluster 6")

#mean((assigned probability - observed indicator)^2)
#MSE30_6clustprob <- mean((expit(as.numeric(m30_6clustprob$pooled$estimate %*% t(data.frame(1,CSS_track$ilrX6.orthogonal1.1[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.2[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.3[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.4[CSS_track$imputation!=0],CSS_track$ilrX6.orthogonal1.5[CSS_track$imputation!=0],
#                                                                                         CSS_track$age[CSS_track$imputation!=0],CSS_track$gender[CSS_track$imputation!=0]=="Male",CSS_track$gender[CSS_track$imputation!=0]=="Other",CSS_track$education[CSS_track$imputation!=0] == "Primary school",CSS_track$education[CSS_track$imputation!=0] == "medium cycle higher education",
#                                                                                         CSS_track$education[CSS_track$imputation!=0] == "Other",CSS_track$education[CSS_track$imputation!=0] == "Primary school",CSS_track$education[CSS_track$imputation!=0] == "short cycle higher education",CSS_track$education[CSS_track$imputation!=0] == "Technical vocational education",
#                                                                                         CSS_track$education[CSS_track$imputation!=0] == "Upper secondary education", CSS_track$occupation[CSS_track$imputation!=0] == "Other",CSS_track$occupation[CSS_track$imputation!=0] == "outside labor market",CSS_track$occupation[CSS_track$imputation!=0] == "student",CSS_track$occupation[CSS_track$imputation!=0] == "unemployed")))
#)-(CSS_track$bmi[CSS_track$imputation!=0]>=30))^2)

m <- glm((bmi>=30) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial, data = CSS_track[CSS_track$imputation==1,])
m$coefficients <- beta_30_1$pooled$estimate

MSE30_6clustprob <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=30))^2)


m30_4clustprob <- pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial)))

beta_30_1 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_2 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_3 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]
beta_30_4 <- summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[2,]

four.cluster.beta_30 <- rbind(beta_30_1,beta_30_2,beta_30_3,beta_30_4)
four.cluster.beta_30$term <- c("cluster 1","cluster 2", "cluster 3", "cluster 4")

# Prediction 

m <- glm((bmi>=30) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial, data = CSS_track[CSS_track$imputation==1,])
m$coefficients <- beta_30_1$pooled$estimate

MSE30_4clustprob <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=30))^2)


## Maximal assignment - compare with results from above
## 6 clusters
m30_6clustmax <- pool(with(CSS_track_mids,glm((bmi>=30) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial)))
summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:4,]

# Prediction

m <- glm((bmi>=30) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial, data = CSS_track[CSS_track$imputation==1,])
m$coefficients <- m30_6clustmax$pooled$estimate

MSE30_6clustmax <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=30))^2)


## 4 clusters
m30_4clustmax <- pool(with(CSS_track_mids,glm((bmi>=30) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial)))
summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:6,]

# Prediction

m <- glm((bmi>=30) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial, data = CSS_track[CSS_track$imputation==1,])
m$coefficients <- m30_4clustmax$pooled$estimate

MSE30_4clustmax <- mean((expit(predict(m,newdata = CSS_track[CSS_track$imputation!=0,]))- (CSS_track$bmi[CSS_track$imputation!=0]>=30))^2)


## Self score:

## >=25
summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:4,]
summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:4,]

## >=30
summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:4,]
summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)[1:4,]



## Numeric BMI - Not of interest in this sample at present time (as we already have the cross-sectional look at base data)

#Without adjustment for tracking:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("bmi","mobileUseNight","mobileUseBeforeSleep","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoT <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_CSSTrackNoT$qbar
pool_inf_CSSTrackNoT$ubar
pool_inf_CSSTrackNoT$pval

summary(pool_inf_CSSTrackNoT)

lowerNight2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerNight3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[3,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[3,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[3,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerNight4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[4,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[4,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[4,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[5,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[5,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[5,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[6,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[6,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[6,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[7,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[7,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[7,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[8,5]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[8,1]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoT)[8,6]+10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_CSSTrackNoT <- cbind(c(lowerNight2,lowerNight3,lowerNight4,lowerBS2,lowerBS3,lowerBS4,lowerBS5),
                              c(estNight2,estNight3,estNight4,estBS2,estBS3,estBS4,estBS5),
                              c(upperNight2,upperNight3,upperNight4,upperBS2,upperBS3,upperBS4,upperBS5))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


#Trends:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("mobileUseNight","mobileUseBeforeSleep","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoTTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_CSSTrackNoTTrend$qbar
pool_inf_CSSTrackNoTTrend$ubar
pool_inf_CSSTrackNoTTrend$pval

summary(pool_inf_CSSTrackNoTTrend)

lowerCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoTTrend)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatNight <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoTTrend)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoTTrend)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCatBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoTTrend)[3,5]+10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatBS <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoTTrend)[3,1]+10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoTTrend)[3,6]+10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_CSSTrackNoTTrend <- rbind(c(lowerCatNight,estCatNight,upperCatNight) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value,
                                   c(lowerCatBS,estCatBS,upperCatBS) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value)


#Without adjustment for selfScore (i.e. with adjustment for tracking) with six clusters: 

#Reference 1

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal1.1","ilrX6.orthogonal1.2","ilrX6.orthogonal1.3","ilrX6.orthogonal1.4","ilrX6.orthogonal1.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust1 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 2

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal2.1","ilrX6.orthogonal2.2","ilrX6.orthogonal2.3","ilrX6.orthogonal2.4","ilrX6.orthogonal2.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 3

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal3.1","ilrX6.orthogonal3.2","ilrX6.orthogonal3.3","ilrX6.orthogonal3.4","ilrX6.orthogonal3.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 4

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal4.1","ilrX6.orthogonal4.2","ilrX6.orthogonal4.3","ilrX6.orthogonal4.4","ilrX6.orthogonal4.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


# Reference 5

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal5.1","ilrX6.orthogonal5.2","ilrX6.orthogonal5.3","ilrX6.orthogonal5.4","ilrX6.orthogonal5.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 6

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal6.1","ilrX6.orthogonal6.2","ilrX6.orthogonal6.3","ilrX6.orthogonal6.4","ilrX6.orthogonal6.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

##

confints_CSSTrackNoS_Six <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6),
                                  c(estClust1,estClust2,estClust3,estClust4,estClust5,estClust6),
                                  c(upperClust1,upperClust2,upperClust3,upperClust4,upperClust5,upperClust6))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

## prediction error - matching parametrizations

m <- gamlss(bmi~(ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX6.orthogonal6.1","ilrX6.orthogonal6.2","ilrX6.orthogonal6.3","ilrX6.orthogonal6.4","ilrX6.orthogonal6.5","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_CSSTrackNoS$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_CSSTrackNoS$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_CSSTrackNoS$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmiCSS <- predict(m, newdata = CSS_track[CSS_track$imputation!=0,c("ilrX6.orthogonal6.1","ilrX6.orthogonal6.2","ilrX6.orthogonal6.3","ilrX6.orthogonal6.4","ilrX6.orthogonal6.5","bmi","age","gender","education","occupation","sample_weights","imputation")])

MSEbmiCSSsixprob <- mean((predbmiCSS - CSS_track$bmi[CSS_track$imputation!=0])^2)


## Without adjustment for selfScore (i.e. with adjustment for tracking) with four clusters:

#Reference 1 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX4.orthogonal1.1","ilrX4.orthogonal1.2","ilrX4.orthogonal1.3","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust1 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 2 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX4.orthogonal2.1","ilrX4.orthogonal2.2","ilrX4.orthogonal2.3","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 3 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX4.orthogonal3.1","ilrX4.orthogonal3.2","ilrX4.orthogonal3.3","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 4 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX4.orthogonal4.1","ilrX4.orthogonal4.2","ilrX4.orthogonal4.3","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

##

confints_CSSTrackNoS_Four <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4),
                                   c(estClust1,estClust2,estClust3,estClust4),
                                   c(upperClust1,upperClust2,upperClust3,upperClust4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

## prediction error - matching parametrizations

m <- gamlss(bmi~(ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("ilrX4.orthogonal4.1","ilrX4.orthogonal4.2","ilrX4.orthogonal4.3","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_CSSTrackNoS$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_CSSTrackNoS$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_CSSTrackNoS$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmiCSS <- predict(m, newdata = CSS_track[CSS_track$imputation!=0,c("ilrX4.orthogonal4.1","ilrX4.orthogonal4.2","ilrX4.orthogonal4.3","bmi","age","gender","education","occupation","sample_weights","imputation")])

MSEbmiCSSfourprob  <- mean((predbmiCSS - CSS_track$bmi[CSS_track$imputation!=0])^2)

## Level curves of prediction surface on ternary plot
pop_data <- #constructed prediction surface for level curves, take pop_track and change the probs. Then use predict.
  ggtern(data=pop_data,aes(x=prob1,y=prob2, z=prob3)) +
  geom_point() #color ~ bmi


## Maximal posterior probability assignment 6 clusters

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}
pool_inf_CSSTrackNoS_mp <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[3,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[3,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[3,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[4,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[4,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[4,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[5,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[5,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[5,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[6,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[6,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[6,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_CSSTrackNoS_mpSix <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6),
                                    c(estClust1,estClust2,estClust3,estClust4,estClust5,estClust6),
                                    c(upperClust1,upperClust2,upperClust3,upperClust4,upperClust5,upperClust6))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

## prediction

m <- gamlss(bmi~(cluster+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_CSSTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_CSSTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_CSSTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmiCSS <- predict(m, newdata = CSS_track[CSS_track$imputation!=0,c("cluster","bmi","age","gender","education","occupation","sample_weights","imputation")])

MSEbmiCSSsixmax  <- mean((predbmiCSS - CSS_track$bmi[CSS_track$imputation!=0])^2)


## Maximal posterior probability assignment 4 clusters

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster.y+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster.y","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_CSSTrackNoS_mp <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[2,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[2,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[2,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[3,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[3,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[3,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[4,5]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[4,1]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_CSSTrackNoS_mp)[4,6]+10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_CSSTrackNoS_mpFour <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4),
                                     c(estClust1,estClust2,estClust3,estClust4),
                                     c(upperClust1,upperClust2,upperClust3,upperClust4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_CSSTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

## prediction

m <- gamlss(bmi~(cluster.y+age+gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster.y","bmi","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_CSSTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_CSSTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_CSSTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmiCSS <- predict(m, newdata = CSS_track[CSS_track$imputation!=0,c("cluster.y","bmi","age","gender","education","occupation","sample_weights","imputation")])

MSEbmiCSSfourmax  <- mean((predbmiCSS - CSS_track$bmi[CSS_track$imputation!=0])^2)



# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#Tracking data: Population sample (random sample) - analyses of numerical bmi and of indicators

hist(pop_track$bmi,breaks=50,xlim=c(0,50))

ggplot(pop_track, aes(x = factor(selfScoreCat))) +
  geom_bar()
ggplot(pop_track, aes(x = factor(cluster))) +
  geom_bar()

#analyses

## regression analysis of clusters of night-time smartphone use and overweight/obesity #justeres for selfScoreCat?? ## hvad er de forskellige clusters?? ## fortolkning?? ## inkluderer imp_nr=0?

## Continuous Outcome


#Without adjustment for tracking:

## Night:

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){ #Slow
  m <- gamlss(bmi~(mobileUseNight+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("bmi","mobileUseNight","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoT <- miceadds::pool_mi(qhat = coefs, u = vcovs)


lowerNight2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerNight3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerNight4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estNight4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperNight4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoT <- cbind(c(lowerNight2,lowerNight3,lowerNight4),
                                 c(estNight2,estNight3,estNight4,estBS2),
                                 c(upperNight2,upperNight3,upperNight4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


#Trends:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(as.numeric(mobileUseNight)+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("mobileUseNight","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatNight <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatNight <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackNoTTrend <- rbind(c(lowerCatNight,estCatNight,upperCatNight) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value)

## BeforeSleep

#Without adjustment for tracking:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){ #Slow
  m <- gamlss(bmi~(mobileUseBeforeSleep+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("bmi","mobileUseBeforeSleep","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoT <- miceadds::pool_mi(qhat = coefs, u = vcovs)


lowerBS2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerBS5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[5,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estBS5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[5,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperBS5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[5,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackNoT <- cbind(c(lowerBS2,lowerBS3,lowerBS4,lowerBS5),
                              c(estBS2,estBS3,estBS4,estBS5),
                              c(upperBS2,upperBS3,upperBS4,upperBS5))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


#Trends:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(as.numeric(mobileUseBeforeSleep)+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("mobileUseBeforeSleep","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerCatBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCatBS <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCatBS <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoTTrend <- rbind(c(lowerCatBS,estCatBS,upperCatBS) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value)


#Without adjustment for selfScore (i.e. with adjustment for tracking) with six clusters: 

#Reference 1

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal1.1","ilrX6.orthogonal1.2","ilrX6.orthogonal1.3","ilrX6.orthogonal1.4","ilrX6.orthogonal1.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust1 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 2

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal2.1","ilrX6.orthogonal2.2","ilrX6.orthogonal2.3","ilrX6.orthogonal2.4","ilrX6.orthogonal2.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 3

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal3.1","ilrX6.orthogonal3.2","ilrX6.orthogonal3.3","ilrX6.orthogonal3.4","ilrX6.orthogonal3.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 4

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal4.1","ilrX6.orthogonal4.2","ilrX6.orthogonal4.3","ilrX6.orthogonal4.4","ilrX6.orthogonal4.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


# Reference 5

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal5.1","ilrX6.orthogonal5.2","ilrX6.orthogonal5.3","ilrX6.orthogonal5.4","ilrX6.orthogonal5.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 6

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal6.1","ilrX6.orthogonal6.2","ilrX6.orthogonal6.3","ilrX6.orthogonal6.4","ilrX6.orthogonal6.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

##

confints_PopTrackNoS_Six <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6),
                                 c(estClust1,estClust2,estClust3,estClust4,estClust5,estClust6),
                                 c(upperClust1,upperClust2,upperClust3,upperClust4,upperClust5,upperClust6))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

## prediction error - matching parametrizations

m <- gamlss(bmi~(ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX6.orthogonal6.1","ilrX6.orthogonal6.2","ilrX6.orthogonal6.3","ilrX6.orthogonal6.4","ilrX6.orthogonal6.5","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_PopTrackNoS$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmipop_six <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("ilrX6.orthogonal6.1","ilrX6.orthogonal6.2","ilrX6.orthogonal6.3","ilrX6.orthogonal6.4","ilrX6.orthogonal6.5","bmi","age","Gender","education","occupation","sample_weights","imputation")])

MSEbmipopsixprob <- mean((predbmipop_six - pop_track$bmi[pop_track$imputation!=0])^2)


## Without adjustment for selfScore (i.e. with adjustment for tracking) with four clusters:

#Reference 1 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX4.orthogonal1.1","ilrX4.orthogonal1.2","ilrX4.orthogonal1.3","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust1 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust1 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 2 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX4.orthogonal2.1","ilrX4.orthogonal2.2","ilrX4.orthogonal2.3","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 3 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX4.orthogonal3.1","ilrX4.orthogonal3.2","ilrX4.orthogonal3.3","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

# Reference 4 (4)

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX4.orthogonal4.1","ilrX4.orthogonal4.2","ilrX4.orthogonal4.3","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

##

confints_PopTrackNoS_Four <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4),
                              c(estClust1,estClust2,estClust3,estClust4),
                              c(upperClust1,upperClust2,upperClust3,upperClust4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

## prediction error - matching parametrizations

m <- gamlss(bmi~(ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("ilrX4.orthogonal4.1","ilrX4.orthogonal4.2","ilrX4.orthogonal4.3","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==1)),family = BCCG)

m$mu.coefficients <- pool_inf_PopTrackNoS$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmipop_four <- predict(m,type="response",newdata = pop_track[pop_track$imputation!=0,c("ilrX4.orthogonal4.1","ilrX4.orthogonal4.2","ilrX4.orthogonal4.3","bmi","age","Gender","education","occupation","sample_weights","imputation")])

MSEbmipopfourprob <- mean((predbmipop_four - pop_track$bmi[pop_track$imputation!=0])^2)


## Level curves of prediction surface on ternary plot - four clusters

pop_plot <- pop_track[pop_track$imputation==1,c("bmi","age","Gender","education","occupation","sample_weights","imputation")]#constructed prediction surface for level curves, take pop_track and change the probs (fix other covariates). Then use predict.
pop_plot$Gender <- "Male"
pop_plot$age <- 30
pop_plot$education <- "long cycle higher education"
pop_plot$occupation <- "employed"

pop_plot <- pop_plot[1:2400,-c(6,7)]
pop_plot$prob1 <- 0
pop_plot$prob2 <- c(rep(0.01,600),rep(0.25,600),rep(0.50,600),rep(0.75,600))
pop_plot$prob3 <- 0
pop_plot$prob4 <- 0

pop_plot$prob1[1:200] <- runif(n=200,0,(1-pop_plot$prob2)[1:200])
pop_plot$prob3[1:200] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob2)[1:200])
pop_plot$prob4[1:200] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob3-pop_plot$prob2)[1:200])

pop_plot$prob3[201:400] <- runif(n=200,0,(1-pop_plot$prob2)[201:400])
pop_plot$prob4[201:400] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob2)[201:400])
pop_plot$prob1[201:400] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob4-pop_plot$prob2)[201:400])

pop_plot$prob4[401:600] <- runif(n=200,0,(1-pop_plot$prob2)[401:600])
pop_plot$prob1[401:600] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob2)[401:600])
pop_plot$prob3[401:600] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob1-pop_plot$prob2)[401:600])

#
pop_plot$prob1[601:800] <- runif(n=200,0,(1-pop_plot$prob2)[601:800])
pop_plot$prob3[601:800] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob2)[601:800])
pop_plot$prob4[601:800] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob3-pop_plot$prob2)[601:800])

pop_plot$prob3[801:1000] <- runif(n=200,0,(1-pop_plot$prob2)[801:1000])
pop_plot$prob4[801:1000] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob2)[801:1000])
pop_plot$prob1[801:1000] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob4-pop_plot$prob2)[801:1000])

pop_plot$prob4[1001:1200] <- runif(n=200,0,(1-pop_plot$prob2)[1001:1200])
pop_plot$prob1[1001:1200] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob2)[1001:1200])
pop_plot$prob3[1001:1200] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob1-pop_plot$prob2)[1001:1200])

#
pop_plot$prob1[1201:1400] <- runif(n=200,0,(1-pop_plot$prob2)[1201:1400])
pop_plot$prob3[1201:1400] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob2)[1201:1400])
pop_plot$prob4[1201:1400] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob3-pop_plot$prob2)[1201:1400])

pop_plot$prob3[1401:1600] <- runif(n=200,0,(1-pop_plot$prob2)[1401:1600])
pop_plot$prob4[1401:1600] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob2)[1401:1600])
pop_plot$prob1[1401:1600] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob4-pop_plot$prob2)[1401:1600])

pop_plot$prob4[1601:1800] <- runif(n=200,0,(1-pop_plot$prob2)[1601:1800])
pop_plot$prob1[1601:1800] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob2)[1601:1800])
pop_plot$prob3[1601:1800] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob1-pop_plot$prob2)[1601:1800])

#
pop_plot$prob1[1801:2000] <- runif(n=200,0,(1-pop_plot$prob2)[1801:2000])
pop_plot$prob3[1801:2000] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob2)[1801:2000])
pop_plot$prob4[1801:2000] <- runif(n=200,0,(1-pop_plot$prob1-pop_plot$prob3-pop_plot$prob2)[1801:2000])

pop_plot$prob3[2001:2200] <- runif(n=200,0,(1-pop_plot$prob2)[2001:2200])
pop_plot$prob4[2001:2200] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob2)[2001:2200])
pop_plot$prob1[2001:2200] <- runif(n=200,0,(1-pop_plot$prob3-pop_plot$prob4-pop_plot$prob2)[2001:2200])

pop_plot$prob4[2201:2400] <- runif(n=200,0,(1-pop_plot$prob2)[2201:2400])
pop_plot$prob1[2201:2400] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob2)[2201:2400])
pop_plot$prob3[2201:2400] <- runif(n=200,0,(1-pop_plot$prob4-pop_plot$prob1-pop_plot$prob2)[2201:2400])


ilr_coord <- Orthogonal.coordinate.generator(pop_plot[,c("prob1","prob2","prob3","prob4")],4)
pop_plot$ilrX4.orthogonal4.1 <- ilr_coord[,1]
pop_plot$ilrX4.orthogonal4.2 <- ilr_coord[,2]
pop_plot$ilrX4.orthogonal4.3 <- ilr_coord[,3]  

preds <- predict(m,newdata=pop_plot[,-c(6,7,8,9)])

ggtern(data=pop_plot[1:600,],aes(x=prob1,y=prob3, z=prob4)) +
  geom_point(size=3.0,aes(color=preds[1:600])) +scale_color_gradient2(midpoint=24, low="blue", mid="orange", high="darkred") #color ~ predbmi

ggtern(data=pop_plot[601:1200,],aes(x=prob1,y=prob3, z=prob4)) +
  geom_point(size=3.0,aes(color=preds[601:1200]))+scale_color_gradient2(midpoint=24, low="blue", mid="orange", high="darkred") #color ~ bmi

ggtern(data=pop_plot[1201:1800,],aes(x=prob1,y=prob3, z=prob4)) +
  geom_point(size=3.0,aes(color=preds[1201:1800]))+scale_color_gradient2(midpoint=24, low="blue", mid="orange", high="darkred") #color ~ bmi

ggtern(data=pop_plot[1801:2400,],aes(x=prob1,y=prob3, z=prob4)) +
  geom_point(size=3.0,aes(color=preds[1801:2400]))+scale_color_gradient2(midpoint=24, low="blue", mid="orange", high="darkred") #color ~ bmi



## Maximal posterior probability assignment 6 clusters

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}
pool_inf_PopTrackNoS_mp <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[5,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[5,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[5,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[6,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[6,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[6,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoS_mpSix <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6),
                              c(estClust1,estClust2,estClust3,estClust4,estClust5,estClust6),
                              c(upperClust1,upperClust2,upperClust3,upperClust4,upperClust5,upperClust6))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
## prediction

m <- gamlss(bmi~(cluster+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmipop_sixmax <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("cluster","bmi","age","Gender","education","occupation","sample_weights","imputation")])

MSEbmipopsixmax <- mean((predbmipop_sixmax - pop_track$bmi[pop_track$imputation!=0])^2)


## Maximal posterior probability assignment 4 clusters

coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster.y+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster.y","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS_mp <- miceadds::pool_mi(qhat = coefs, u = vcovs)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[3,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,5]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,1]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS_mp)[4,6]+10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoS_mpFour <- cbind(c(lowerClust1,lowerClust2,lowerClust3,lowerClust4),
                              c(estClust1,estClust2,estClust3,estClust4),
                              c(upperClust1,upperClust2,upperClust3,upperClust4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS_mp$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
## prediction

m <- gamlss(bmi~(cluster.y+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster.y","bmi","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)

m$mu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))] 
m$nu.coefficients <- pool_inf_PopTrackNoS_mp$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))] 

predbmipop_fourmax <- predict(m, newdata = pop_track[pop_track$imputation!=0,c("cluster.y","bmi","age","Gender","education","occupation","sample_weights","imputation")])

MSEbmipopfourmax <- mean((predbmipop_fourmax - pop_track$bmi[pop_track$imputation!=0])^2)


# --------------------------------------------------------------------------- ##
### Binary Outcomes for population sample
# --------------------------------------------------------------------------- ##

## BMI > 25

## clusters and BMI>25 (no adjustment for selfScoreCat) - six clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_1 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_2 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_3 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_4 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_5 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_6 <- modelRandom25No[2,]

beta_25_random_sum <- rbind(beta_25_random_1,beta_25_random_2,beta_25_random_3,beta_25_random_4,beta_25_random_5,beta_25_random_6)

# prediction - matching parametrizations
m <- glm((bmi>=25) ~ (ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random25No)$pooled$estimate
predpopbin25_six <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin25_predsix <- mean((expit(predpopbin25_six)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)


## Four clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_1 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_2 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_3 <- modelRandom25No[2,]

Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int = T)
beta_25_random_4 <- modelRandom25No[2,]

beta_25_random_sum <- rbind(beta_25_random_1,beta_25_random_2,beta_25_random_3,beta_25_random_4)

# prediction - matching parametrizations
m <- glm((bmi>=25) ~ (ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random25No)$pooled$estimate
predpopbin25_four <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin25_predfour <- mean((expit(predpopbin25_four)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)


## Maximal posterior probability assignment


# Six clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No_mpSix <- summary(pool(Random25No), conf.int = T)
m <- glm((bmi>=25) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random25No)$pooled$estimate
predpopbin25_maxsix <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin25_predmaxsix <- mean((expit(predpopbin25_maxsix)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)

# Four clusters
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No_mpFour <- summary(pool(Random25No), conf.int = T)
m <- glm((bmi>=25) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random25No)$pooled$estimate
predpopbin25_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin25_predmaxfour <- mean((expit(predpopbin25_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=25))^2)



## SelfScoreCat and BMI>25 (no adjustment for tracking)
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25NoT <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25NoT <- summary(pool(Random25NoT), conf.int=T)
cbind(exp(modelRandom25NoT$estimate),
exp(modelRandom25NoT$`2.5 %`),
exp(modelRandom25NoT$`97.5 %`))

#test for trend (selfscoreCat uden cluster)
Random25NoTest <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25NoTest), conf.int=T)

# --------------------------------------------------------------------------- ##
## BMI > 30

## BMI >30 

## clusters and BMI>30 (no adjustment for selfScoreCat) - six clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_1 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_2 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_3 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_4 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_5 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_6 <- modelRandom30No[2,]

beta_30_random_sum <- rbind(beta_30_random_1,beta_30_random_2,beta_30_random_3,beta_30_random_4,beta_30_random_5,beta_30_random_6)

# prediction - matching parametrizations
m <- glm((bmi>=30) ~ (ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random30No)$pooled$estimate
predpopbin30_six <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin30_predsix <- mean((expit(predpopbin30_six)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)


## Four clusters

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_1 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_2 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_3 <- modelRandom30No[2,]

Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
beta_30_random_4 <- modelRandom30No[2,]

beta_30_random_sum <- rbind(beta_30_random_1,beta_30_random_2,beta_30_random_3,beta_30_random_4)

# prediction - matching parametrization
m <- glm((bmi>=30) ~ (ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random30No)$pooled$estimate
predpopbin30_four <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin30_predfour <- mean((expit(predpopbin30_four)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)


## Maximal posterior probability assignment

# Six clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No_mpSix <- summary(pool(Random30No), conf.int = T)
m <- glm((bmi>=30) ~ (cluster+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random30No)$pooled$estimate
predpopbin30_maxsix <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin30_predmaxsix <- mean((expit(predpopbin30_maxsix)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)

# Four clusters
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No_mpFour <- summary(pool(Random30No), conf.int = T)
m <- glm((bmi>=30) ~ (cluster.y+age+gender+education+occupation), weights=sample_weights,family=binomial, data=pop_track[pop_track$imputation==1,])
m$coefficients <- pool(Random30No)$pooled$estimate
predpopbin30_maxfour <- predict(m,newdata = pop_track[pop_track$imputation!=0,])
MSEpopbin30_predmaxfour <- mean((expit(predpopbin30_maxfour)-(pop_track$bmi[pop_track$imputation!=0]>=30))^2)


## no adjustment for tracking
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30NoT <- with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30NoT <- summary(pool(Random30NoT), conf.int=T)
cbind(exp(modelRandom30NoT$estimate),
exp(modelRandom30NoT$`2.5 %`),
exp(modelRandom30NoT$`97.5 %`))

#test for trend (selfScoreCat as numeric and no adjustment for clusters)
Random30NoT <- with(pop_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseNight)+as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random30NoT), conf.int = T)



###############################################################################
###############################################################################
###############################################################################

#Analysis of the clinical sample data - interest in biomarkers

# --------------------------------------------------------------------------- ##
## BMI
table(clinical_sample$bmi.clinical)
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
summary(clinical_sample$bmi.clinical, useNA="always")

## BMI clinical
table(clinical_sample$bmi.clinical)
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable( ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

# --------------------------------------------------------------------------- ##
## SelfScoreCat

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
table(clinical_sample$selfScoreCat[clinical_sample$imputation!=0],useNA="always")/25
clinical_sample$selfScoreCat <- as.factor(clinical_sample$selfScoreCat)

# --------------------------------------------------------------------------- ##
## descriptive of clinical sample
## age
publish(univariateTable(selfScoreCat ~ age.x,data=clinical_sample, column.percent=TRUE))/25
clinical_sample$age.x <- as.numeric(clinical_sample$age.x)

## BMI
publish(univariateTable(selfScoreCat ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

## categorize BMI
clinical_sample$bmiCat[clinical_sample$bmi.clinical<25] <- "<25"
clinical_sample$bmiCat[clinical_sample$bmi.clinical>=25 & clinical_sample$bmi.clinical<30] <- "25-30"
clinical_sample$bmiCat[clinical_sample$bmi.clinical>=30] <- ">=30"
table(clinical_sample$bmiCat, useNA="always")/26
prop.table(table(clinical_sample$bmiCat, useNA="always"))
table(clinical_sample$selfScoreCat[clinical_sample$imputation!=0], clinical_sample$bmiCat[clinical_sample$imputation!=0], useNA="always")/25
prop.table(table(clinical_sample$selfScoreCat[clinical_sample$imputation!=0], clinical_sample$bmiCat[clinical_sample$imputation!=0], useNA="always"))
publish(univariateTable(selfScoreCat ~ bmiCat,data=clinical_sample, column.percent=TRUE))


#The subjects are scoring in the high end. Is this an issue or a characteristic of the data?

#Introducing interesting derived variables

clinical_sample$bmi <- as.numeric(clinical_sample$bmi.clinical)
clinical_sample$bmi25 <- as.numeric(clinical_sample$bmi.clinical>=25)
clinical_sample$bmi30 <- as.numeric(clinical_sample$bmi.clinical>=30)

publish(univariateTable(selfScoreCat ~ bmi25,data=clinical_sample, column.percent=TRUE))

## age at clinical examination (OBS. NOGET ER GALT I DENNE VARIABLE)
table(clinical_sample$age.y, useNA="always")
clinical_sample$age<- as.numeric(str_c(substr(clinical_sample$age.y,1,1),substr(clinical_sample$age.y,2+(mod(nchar(clinical_sample$age.y),4)==1),2+(mod(nchar(clinical_sample$age.y),4)==1)),".",
                 substr(clinical_sample$age.y,3+(mod(nchar(clinical_sample$age.y),4)!=3),3+(mod(nchar(clinical_sample$age.y),4)!=3))))


# --------------------------------------------------------------------------- ##

## systolic blood pressure
clinical_sample$sbp<-rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ sbp,data=clinical_sample, column.percent=TRUE))

# diastolic blood pressure
clinical_sample$dbp<-rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T)
publish(univariateTable(selfScoreCat ~ dbp,data=clinical_sample, column.percent=TRUE))

## hip waist ratio

clinical_sample$ratiowaisthip <- as.numeric(clinical_sample$ratiowaisthip)
publish(univariateTable(selfScoreCat ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))


## bmi clinical
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable(selfScoreCat ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))


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
hist(as.numeric(clinical_sample$bmi.clinical),breaks=20)
hist(rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T),breaks=20)
hist(rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T),breaks=20)

## HDL
clinical_sample$hdl <- as.numeric(clinical_sample$hdl)
publish(univariateTable(selfScoreCat ~ hdl,data=clinical_sample, column.percent=TRUE))

## LDL
clinical_sample$ldl <- as.numeric(clinical_sample$ldl)
publish(univariateTable( selfScoreCat~ ldl,data=clinical_sample, column.percent=TRUE))

## VLDL
clinical_sample$vldl <- as.numeric(clinical_sample$vldl)
publish(univariateTable( selfScoreCat~ vldl,data=clinical_sample, column.percent=TRUE))

## total cholesterol
clinical_sample$t_cholesterol <- as.numeric(clinical_sample$t_cholesterol)
publish(univariateTable(selfScoreCat ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))

## triglycerides
clinical_sample$triglycerids <- as.numeric(clinical_sample$triglycerids)
publish(univariateTable(selfScoreCat ~ triglycerids,data=clinical_sample, column.percent=TRUE))

## hba1c
clinical_sample$hba1c <- as.numeric(clinical_sample$hba1c)
publish(univariateTable(selfScoreCat ~ hba1c,data=clinical_sample, column.percent=TRUE))


# --------------------------------------------------------------------------- ##

#Models - multiple testing issue if we are going to 'pick and choose' which responses we would like to look at.
table(clinical_sample$age.x)

#Transforming to mids for modelling and inference

clinical_sample$cluster.y <- factor(clinical_sample$cluster.y,levels = c("Cluster 3", "Cluster 2", "Cluster 4", "Cluster 1"))

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")

#hdl
#hdl_sum<-cbind(summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
lines(seq(from=min(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),length.out=100),dnorm(x=seq(from=min(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),length.out=100),mean=mean(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))),sd=sd(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))))
plot(residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))


#ldl 
#ldl_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res=residuals(lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))
plot(fitted(lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

cbind(confint(glm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#vldl
#vldl_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),
      cbind(coef(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))-1.96*summary(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2],
            coef(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))+1.96*summary(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2]))


#t_cholesterol
#t_cholesterol_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#                summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#triglycerids
#tri_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(triglycerids) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, glm(as.numeric(triglycerids) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(glm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(glm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(triglycerids) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),
      cbind(coef(glm(as.numeric(triglycerids) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))-1.96*summary(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2],
            coef(glm(as.numeric(triglycerids) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))+1.96*summary(glm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2]))

#hba1c
#hba1c_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res <- residuals(lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

# bmi
hist(residuals(lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))


#ratiowaisthip
#wh_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation))))$std.error[2:4])
hist(residuals(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#sbp
#sbp_sum<-cbind(summary(pool(with(data=clinical_mids, lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(sbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(sbp) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(sbp) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#dbp
#dbp_sum<-cbind(summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(dbp) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(dbp) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#Generally the residuals look reasonably centered, with a few positive outliers. The residual distributions on the first imputation actually look reasonably normal, save for the few (extreme) outliers.

#Seems that these models are appropriate, and that normal approximations of confidence interval will be reasonable too. We can however also just use the profile likelihood CI's.
#The vldl measure is the only one that has a right skewed distribution. (Do something particular about it?)

#Wald confidence intervals

## Tracking: Six clusters 

dbp_int1 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int2 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int3 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int4 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int5 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int6 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])

dbp_int <- rbind(dbp_int1,dbp_int2,dbp_int3,dbp_int4,dbp_int5,dbp_int6)

glu_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])

glu_int <- rbind(glu_int1,glu_int2,glu_int3,glu_int4,glu_int5,glu_int6)

hba1c_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])

hba1c_int <- rbind(hba1c_int1,hba1c_int2,hba1c_int3,hba1c_int4,hba1c_int5,hba1c_int6)

hdl_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])

hdl_int <- rbind(hdl_int1,hdl_int2,hdl_int3,hdl_int4,hdl_int5,hdl_int6)

ldl_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])

ldl_int <- rbind(ldl_int1,ldl_int2,ldl_int3,ldl_int4,ldl_int5,ldl_int6)

t_chol_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])

t_chol_int <- rbind(t_chol_int1,t_chol_int2,t_chol_int3,t_chol_int4,t_chol_int5,t_chol_int6)

sbp_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])

sbp_int <- rbind(sbp_int1,sbp_int2,sbp_int3,sbp_int4,sbp_int5,sbp_int6)

tri_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])

tri_int <- rbind(tri_int1,tri_int2,tri_int3,tri_int4,tri_int5,tri_int6)

vldl_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])

vldl_int <- rbind(vldl_int1,vldl_int2,vldl_int3,vldl_int4,vldl_int5,vldl_int6)

wh_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

wh_int <- rbind(wh_int1,wh_int2,wh_int3,wh_int4,wh_int5,wh_int6)

bmi_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal1.1+ilrX6.orthogonal1.2+ilrX6.orthogonal1.3+ilrX6.orthogonal1.4+ilrX6.orthogonal1.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal2.1+ilrX6.orthogonal2.2+ilrX6.orthogonal2.3+ilrX6.orthogonal2.4+ilrX6.orthogonal2.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal3.1+ilrX6.orthogonal3.2+ilrX6.orthogonal3.3+ilrX6.orthogonal3.4+ilrX6.orthogonal3.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal4.1+ilrX6.orthogonal4.2+ilrX6.orthogonal4.3+ilrX6.orthogonal4.4+ilrX6.orthogonal4.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int5 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal5.1+ilrX6.orthogonal5.2+ilrX6.orthogonal5.3+ilrX6.orthogonal5.4+ilrX6.orthogonal5.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int6 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX6.orthogonal6.1+ilrX6.orthogonal6.2+ilrX6.orthogonal6.3+ilrX6.orthogonal6.4+ilrX6.orthogonal6.5+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- rbind(bmi_int1,bmi_int2,bmi_int3,bmi_int4,bmi_int5,bmi_int6)

df_ints_Six <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_Six) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")

## Tracking: Four clusters

dbp_int1 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int2 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int3 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
dbp_int4 <- summary(pool(with(data=clinical_mids, lm(dbp ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])

dbp_int <- rbind(dbp_int1,dbp_int2,dbp_int3,dbp_int4)

glu_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
glu_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])

glu_int <- rbind(glu_int1,glu_int2,glu_int3,glu_int4)

hba1c_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hba1c_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])

hba1c_int <- rbind(hba1c_int1,hba1c_int2,hba1c_int3,hba1c_int4)

hdl_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
hdl_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])

hdl_int <- rbind(hdl_int1,hdl_int2,hdl_int3,hdl_int4)

ldl_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
ldl_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])

ldl_int <- rbind(ldl_int1,ldl_int2,ldl_int3,ldl_int4)

t_chol_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])
t_chol_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_chol_sum[,1]-1.96*t_chol_sum[,2],t_chol_sum[,1]+1.96*t_chol_sum[,2])

t_chol_int <- rbind(t_chol_int1,t_chol_int2,t_chol_int3,t_chol_int4)

sbp_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
sbp_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])

sbp_int <- rbind(sbp_int1,sbp_int2,sbp_int3,sbp_int4)

tri_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
tri_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])

tri_int <- rbind(tri_int1,tri_int2,tri_int3,tri_int4)

vldl_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
vldl_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])

vldl_int <- rbind(vldl_int1,vldl_int2,vldl_int3,vldl_int4)

wh_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
wh_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

wh_int <- rbind(wh_int1,wh_int2,wh_int3,wh_int4)

bmi_int1 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX4.orthogonal1.1+ilrX4.orthogonal1.2+ilrX4.orthogonal1.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int2 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX4.orthogonal2.1+ilrX4.orthogonal2.2+ilrX4.orthogonal2.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int3 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX4.orthogonal3.1+ilrX4.orthogonal3.2+ilrX4.orthogonal3.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])
bmi_int4 <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ ilrX4.orthogonal4.1+ilrX4.orthogonal4.2+ilrX4.orthogonal4.3+age+education+occupation,na.action=na.omit))),conf.int=T)[2,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- rbind(bmi_int1,bmi_int2,bmi_int3,bmi_int4)

df_ints_Four <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_Four) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")

## Maximal posterior probability assignment: Six clusters

dbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(dbp) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

vldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

df_ints_mpSix <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_mpSix) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")

## Maximal posterior probability assignment: Four clusters

dbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(dbp) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

vldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(vldl) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster.y+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(bmi_sum[,1]-1.96*bmi_sum[,2],bmi_sum[,1]+1.96*bmi_sum[,2])

## estimater for four clusters i clinical sample
df_ints_mpFour <- list(hdl_int,ldl_int,vldl_int,t_chol_int,tri_int,hba1c_int,dbp_int,sbp_int,wh_int,glu_int,bmi_int)

names(df_ints_mpFour) <- c("hdl","ldl","vldl","t_chol","tri","hba1c","dbp","sbp","wh","glu","bmi")


## Selfscore

predict.person <- data.frame("age"=18,education="",occupation="Employed")

dbp_intsNight <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseNight+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intsNight <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseNight+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intsNight <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseNight+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_intsNight <- data.frame(rbind(hdl_intsNight[2,],ldl_intsNight[2,],vldl_intsNight[2,],t_chol_intsNight[2,],tri_intsNight[2,],hba1c_intsNight[2,],dbp_intsNight[2,],sbp_intsNight[2,],wh_intsNight[2,],glu_intsNight[2,],bmi_intsNight[2,]),
                      rbind(hdl_intsNight[3,],ldl_intsNight[3,],vldl_intsNight[3,],t_chol_intsNight[3,],tri_intsNight[3,],hba1c_intsNight[3,],dbp_intsNight[3,],sbp_intsNight[3,],wh_intsNight[3,],glu_intsNight[3,],bmi_intsNight[3,]),
                      rbind(hdl_intsNight[4,],ldl_intsNight[4,],vldl_intsNight[4,],t_chol_intsNight[4,],tri_intsNight[4,],hba1c_intsNight[4,],dbp_intsNight[4,],sbp_intsNight[4,],wh_intsNight[4,],glu_intsNight[4,],bmi_intsNight[4,]))

colnames(df_intsNight) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsNight) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsNight


## Before

dbp_intsBefore <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intsBefore <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intsBefore <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseBeforeSleep+age+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_intsBefore <- data.frame(rbind(hdl_intsBefore[2,],ldl_intsBefore[2,],vldl_intsBefore[2,],t_chol_intsBefore[2,],tri_intsBefore[2,],hba1c_intsBefore[2,],dbp_intsBefore[2,],sbp_intsBefore[2,],wh_intsBefore[2,],glu_intsBefore[2,],bmi_intsBefore[2,]),
                                rbind(hdl_intsBefore[3,],ldl_intsBefore[3,],vldl_intsBefore[3,],t_chol_intsBefore[3,],tri_intsBefore[3,],hba1c_intsBefore[3,],dbp_intsBefore[3,],sbp_intsBefore[3,],wh_intsBefore[3,],glu_intsBefore[3,],bmi_intsBefore[3,]),
                                rbind(hdl_intsBefore[4,],ldl_intsBefore[4,],vldl_intsBefore[4,],t_chol_intsBefore[4,],tri_intsBefore[4,],hba1c_intsBefore[4,],dbp_intsBefore[4,],sbp_intsBefore[4,],wh_intsBefore[4,],glu_intsBefore[4,],bmi_intsBefore[4,]))

colnames(df_intsBefore) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsBefore) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsBefore
