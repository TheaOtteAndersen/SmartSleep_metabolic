
## Analyses for Theas 2nd PhD paper on night-time smartphone use and obesity and metabolic biomarkers

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
## Coefficients are: Intercept, No activity vs rest, little activity vs more activity, little offset vs little onset, moderate activity vs much activity, moderate offset vs moderate onset
SBP4 <- matrix(c(1,1,-1,1,1,-1,0,1,1,0,0,-1),byrow=T,nrow=3,ncol=4)
## Coefficients are: Intercept, No activity vs rest, little activity vs more activity, offset with activity vs onset with activity


X6 <- subject_tracking_six_clusters[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob")]
X4 <- subject_tracking_four_clusters[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob")]

Phi.generator <- function(X) {
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

## Collecting the two clusterings in one file

subject_tracking_clusters <- left_join(subject_tracking_six_clusters,subject_tracking_four_clusters[,c("userid","cluster","cluster1prob","cluster2prob","cluster3prob","cluster4prob","description","state0prob","state1prob","state2prob","state3prob","ilr1","ilr2","ilr3")],by="userid")
subject_tracking_clusters <- rename(subject_tracking_clusters,ilr1 = ilr1.x, ilr2=ilr2.x, ilr3=ilr3.x,cluster1prob=cluster1prob.x,cluster2prob=cluster2prob.x,cluster3prob=cluster3prob.x,cluster4prob=cluster4prob.x,
                                    state0prob=state0prob.x,state1prob=state1prob.x,state2prob=state2prob.x,state3prob=state3prob.x,cluster=cluster.x,description=description.x)

## load baseline data
setwd("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation")
base_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Experiment/imp_Experiment.csv")

## load followup sample
CSS <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Citizen Science Sample/imp_citizenScience.csv")

## load population sample
pop_data <-read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Population Sample/imp_population.csv")

## load clinical data (survey and clinical data)
clin_data <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/Data imputation/Data/Renset imputation/Clinical Sample/imp_clinical.csv")
clin_clinical <- read.csv2("S:/SUND-IFSV-SmartSleep/Data cleaning/SmartSleep Clinical/Data/Rådata/SmartSleepClinicalData.csv")

# --------------------------------------------------------------------------- ##

#Baseline data with self-reports

## if no mobile phone = NA
base_data$mobileUseBeforeSleep[base_data$mobilephone=="No mobile phone"] <- NA
base_data$mobileUseNight[base_data$mobilephone=="No mobile phone"] <- NA

## Smartphone use before sleep onset
publish(univariateTable( ~ mobileUseBeforeSleep,data=base_data, column.percent=TRUE))
base_data$mobileUseBeforeSleep <- factor(base_data$mobileUseBeforeSleep, levels = c("Never", "Every month or less", "Once a week", "2-4 times per week", "5-7 times per week"))

## Night-time smartphone use
publish(univariateTable( ~ mobileUseNight,data=base_data, column.percent=TRUE))
base_data$mobileUseNight <- factor(base_data$mobileUseNight, levels = c("Never", "A few times a month or less", "A few times a week", "Every night or almost every night"))

## daytime smartphone use
publish(univariateTable( ~ mobileCheck,data=base_data, column.percent=TRUE))

## bmi  
publish(univariateTable( ~ bmi,data=base_data, column.percent=TRUE))
base_data$bmi30 <- (base_data$bmi>=30)
base_data$bmi25 <- (base_data$bmi>=25)

#save(base_data,file="H:/SmartSleep backup IT Issues/gamlssBootstrap/base_data.RData")

base_data_mids <- as.mids(base_data,.imp="imputation")

## BMI kategoriseringer ved baseline
table(base_data$bmi<25, base_data$mobileUseNight)/26
table(base_data$bmi>=25&base_data$bmi<30, base_data$mobileUseNight)/26
table(base_data$bmi>=30, base_data$mobileUseNight)/26

table(base_data$bmi<25, base_data$mobileUseBeforeSleep)/26
table(base_data$bmi>=25&base_data$bmi<30, base_data$mobileUseBeforeSleep)/26
table(base_data$bmi>=30, base_data$mobileUseBeforeSleep)/26


# --------------------------------------------------------------------------- ##
#Followup sample (CSS sample)

## Smartphone use before sleep onset
publish(univariateTable( ~ mobileUseBeforeSleep,data=CSS, column.percent=TRUE))
CSS$mobileUseBeforeSleep <- factor(CSS$mobileUseBeforeSleep, levels = c("Never", "Every month or less", "Once a week", "2-4 times per week", "5-7 times per week"))

## night-time smartphone use
publish(univariateTable( ~ mobileUseNight,data=CSS, column.percent=TRUE))
CSS$mobileUseNight <- factor(CSS$mobileUseNight, levels = c("Never", "A few times a month or less", "A few times a week", "Every night or almost every night"))

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

## smartphone use before sleep onset
table(pop_data$mobileUseBeforeSleep, useNA="always")
pop_data$mobileUseBeforeSleep <- factor(pop_data$mobileUseBeforeSleep, levels = c("Never", "Every month or less", "Once a week", "2-4 times per week", "5-7 times per week"))

## night-time smartphone use 
publish(univariateTable( ~ mobileUseNight,data=pop_data, column.percent=TRUE))
pop_data$mobileUseNight <- factor(pop_data$mobileUseNight, levels = c("Never", "A few times a month or less", "A few times a week", "Every night or almost every night"))

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
plot(fitted(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))),residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))))
hist(residuals(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1))),xlim=c(-20,20),breaks=200)
hist(simulate(lm(bmi~(selfScoreCat+age+gender+education+occupation), weights=sample_weights, data=subset(base_data,imputation==1)))$sim_1,breaks=40) #The bell-shape is not that well suited

#Conclusion: The lm is not by itself appropriate for describing the distribution of BMI.


#Alternative: Pretty good fit. A general family of models.


#Confidence intervals and estimates:

coefs <- list()
ses <- list()
vcovs <- list()
models <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi ~ selfScoreCat+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG) #May use BCS instead of BCCG which corresponds to using a t distribution instead of normal. This can fit heavier tails, though in this case a very large df is fitted, meaning that there is not much difference.
  m_sum <- summary(m)
  models[[i]] <- m
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_base <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_base$qbar
pool_inf_base$ubar
pool_inf_base$ba
pool_inf_base$pval


#Seems that we can get stable contrats of the mean (taking in varying medians), in spite of skewness.

#One slightly hacky way to achieve this may be to take one of the fitted models created by fit() and replace the stored coefficients with the final pooled estimates. I haven't done detailed testing but it seems to be working on this simple example:
m$mu.coefficients <- pool_inf_base$qbar[1:length(m$mu.coefficients)]
m$sigma.coefficients <- pool_inf_base$qbar[(length(m$mu.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients))]
m$nu.coefficients <- pool_inf_base$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1):(length(m$mu.coefficients)+length(m$sigma.coefficients)+length(m$nu.coefficients))]

summary(pool_inf_base,conf.int=T)

#Median contrasts:
mus<- c(predict(m,what="mu",type="response")[(base_data$selfScoreCat=="1" & base_data$age==35 & base_data$gender=="Female" & base_data$education=="long cycle higher education" & base_data$occupation=="employed")[base_data$imputation==i & rowSums(is.na(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)))==0]][1], 
  predict(m,what="mu",type="response")[(base_data$selfScoreCat=="2" & base_data$age==35 & base_data$gender=="Female" & base_data$education=="long cycle higher education" & base_data$occupation=="employed")[base_data$imputation==i & rowSums(is.na(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)))==0]][1], 
  predict(m,what="mu",type="response")[(base_data$selfScoreCat=="3" & base_data$age==35 & base_data$gender=="Female" & base_data$education=="long cycle higher education" & base_data$occupation=="employed")[base_data$imputation==i & rowSums(is.na(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)))==0]][1], 
  predict(m,what="mu",type="response")[(base_data$selfScoreCat=="4" & base_data$age==35 & base_data$gender=="Female" & base_data$education=="long cycle higher education" & base_data$occupation=="employed")[base_data$imputation==i & rowSums(is.na(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)))==0]][1]) 

mus-mus[1]

#So how do we get distribution of means from the BCCG? Can we even get a single number to characterize the means and differences in means?
#Additive median contrasts appear to yield additive mean contrasts. (Why is it so theoretically? Shouldn't the lower value of zero influence the relationsship?)
#The distribution adjusts for the truncation...
#The truncation makes sense to have for BMI.

#Finding particular means (contrasts) by integration
m1 <- integrate(function(y) y*dBCCG(x=y,mu=mus[1],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value  
m2 <- integrate(function(y) y*dBCCG(x=y,mu=mus[2],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value
m3 <- integrate(function(y) y*dBCCG(x=y,mu=mus[3],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value
m4 <- integrate(function(y) y*dBCCG(x=y,mu=mus[4],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value

ms <- c(m1,m2,m3,m4)
c(m2,m3,m4)-m1

plot(mus,ms)

#Now what about confidence regions and p values? We simulate from the fitted BCCG distributions?

#Profile intervals:
lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,5],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,1],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[2,6],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,5],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,1],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[3,6],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,5],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,1],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_base)[4,6],sigma=exp(pool_inf_base$qbar[20]),nu=pool_inf_base$qbar[21]),0,Inf)$value 

confints_base <- cbind(c(lowerCat2,lowerCat3,lowerCat4),c(estCat2,estCat3,estCat4),c(upperCat2,upperCat3,upperCat4))

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
  m <- gamlss(bmi ~ as.numeric(selfScoreCat)+age+gender+education+occupation, sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(base_data[,c("bmi","selfScoreCat","age","gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
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

integrate(function(y) y*dBCCG(x=y,mu=m$mu.coefficients[1],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)  
integrate(function(y) y*dBCCG(x=y,mu=m$mu.coefficients[1]+m$mu.coefficients[2],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)

integrate(function(y) y*dBCCG(x=y,mu=m$mu.coefficients[1],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value  -
  integrate(function(y) y*dBCCG(x=y,mu=m$mu.coefficients[1]+m$mu.coefficients[2],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value

integrate(function(y) y*dBCCG(x=y,mu=m$mu.coefficients[2],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)


#interval
lowerCatTrend <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[2,5],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
estCatTrend <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[2,1],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 
upperCatTrend <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_baseTrend)[2,6],sigma=exp(pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+1)]),nu=pool_inf_baseTrend$qbar[(length(m$mu.coefficients)+length(m$sigma.coefficients)+1)]),0,Inf)$value 

confints_baseTrend <-c(lowerCatTrend,estCatTrend,upperCatTrend)

summary(pool_inf_baseTrend)$p[2]


# --------------------------------------------------------------------------- ##
## cross-sectional associations between risk profiles and bmi (25, 30 & continous)
# --------------------------------------------------------------------------- ##

## Using the mice package with mids objects (night-time smartphone use and bmi)

## virker åbenbart ikke? (02062022)
mod25 <- with(base_data_mids,glm((bmi>=25)~(mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod30 <- with(base_data_mids,glm((bmi>=30)~(mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))

## test for trend
TEST <- with(base_data_mids,glm((bmi>=30)~((as.numeric(mobileUseNight))+age+gender+education+occupation), weights=sample_weights,family=binomial))
testT <- summary(pool(TEST), conf.int = T)

TEST2 <- with(base_data_mids,glm((bmi>=25)~((as.numeric(mobileUseNight))+age+gender+education+occupation), weights=sample_weights,family=binomial))
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

## smartphone use before sleep onset and bmi
mod25Before <- with(base_data_mids,glm(bmi25~(mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
mod30Before <- with(base_data_mids,glm(bmi30~(mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))

## test for trend
## BMI > 25
TEST2 <- with(base_data_mids,glm((bmi>=25)~((as.numeric(mobileUseBeforeSleep))+age+gender+education+occupation), weights=sample_weights,family=binomial))
test2 <- summary(pool(TEST2), conf.int = T)

## BMI > 30
TEST <- with(base_data_mids,glm((bmi>=30)~((as.numeric(mobileUseBeforeSleep))+age+gender+education+occupation), weights=sample_weights,family=binomial))
testT <- summary(pool(TEST), conf.int = T)

## OR for BMI >25
model25Before <- summary(pool(mod25Before), conf.int=T)
cbind(exp(model25$estimate),
      exp(model25$`2.5 %`),
      exp(model25$`97.5 %`))

## OR for BMI>30
model30Before <- summary(pool(mod30Before),conf.int = T)
cbind(exp(model30$estimate),
      exp(model30$`2.5 %`),
      exp(model30$`97.5 %`))


# --------------------------------------------------------------------------- ##
## longitudinal analysis of risk scores of smartphone behavior and changes in BMI
# --------------------------------------------------------------------------- ##
#BMI followup difference - match with emailAddress or CS_ID
#y: base, x: followup


## ----- ##
#Final change analyses 
## ----- ##

## from below 25 to above 25
model25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
model_summary25<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=25 ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))), conf.int = T)
exp(cbind(model_summary25$estimate[1:4],model_summary25$`2.5 %`[1:4],model_summary25$`97.5 %`[1:4]))


## test for trend

test25 <- with(bmi_followup_mids,glm(bmi.fu>=25 ~ (as.numeric(selfScoreCat.y)+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=25), weights=sample_weights.y,family=binomial))
testT25 <- summary(pool(test25), conf.int = T)
#anova(test25,model25)


## from below 30 to above 30
model30 <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
model_summary30<-summary(pool(with(bmi_followup_mids,glm(bmi.fu>=30 ~ (selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))), conf.int = T)
exp(cbind(model_summary30$estimate[1:4],model_summary30$`2.5 %`[1:4],model_summary30$`97.5 %`[1:4]))

test30 <- with(bmi_followup_mids,glm(bmi.fu>=30 ~ (as.numeric(selfScoreCat.y)+age.y+gender.y+education.y+occupation.y+bmi.base)*(bmi.base>=30), weights=sample_weights.y,family=binomial))
testT30 <- summary(pool(test30), conf.int=T)
#anova(test30,model30)



## Modelling numeric difference in bmi between baseline and followup
hist(bmi_followup$difference,xlim=c(-10,10),breaks=600,ylim=c(0,2500))

plot(fitted(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),
     residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))))
hist(residuals(lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y+bmi.base), weights=sample_weights.y, data=subset(bmi_followup,imputation==5))),breaks=50)


#m <- lm(difference~(selfScoreCat.y+age.y+gender.y+education.y+occupation.y)*followup_time-selfScoreCat.y-age.y-gender.y-education.y-occupation.y,weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","selfScoreCat.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights")]))
m <- lm(bmi.fu~(selfScoreCat.y:as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights,data=na.omit(bmi_followup[bmi_followup$imputation==1,c("difference","selfScoreCat.y","age.y","gender.y","education.y","occupation.y","followup_time","sample_weights","bmi.base","bmi.fu")]))

model_summary_diff <- summary(pool(with(bmi_followup_mids,lm(bmi.fu~(selfScoreCat.y:as.numeric(followup_time)+as.numeric(followup_time)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))), conf.int = T)
#Does it make sense to include tracking information at follow up, or to include self score at followup as well?
#(reduce noise in measurement?)

cbind(model_summary_diff$estimate[1:4],model_summary_diff$`2.5 %`[1:4],model_summary_diff$`97.5 %`[1:4])


test_num <- with(bmi_followup_mids,lm(difference/as.numeric(followup_time)~(as.numeric(selfScoreCat.y)+age.y+gender.y+education.y+occupation.y+bmi.base),weights=sample_weights))
test_Tnum <- summary(pool(test_num), conf.int=T)




#Generally:
#Wald intervals with Robust=T are better for misspecified models.
#Profile likelihood intervals are better for models that are close to correct (mostly so for smaller sample sizes).



# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#####Tracking data for the followup CSS sample

#Mice-based inference for models:

#BMI indicators

summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilr1+ilr2+ilr3+ilr4+ilr5+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)

clr.beta <- summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilr1+ilr2+ilr3+ilr4+ilr5+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)$estimate[2:6]%*%Phi6
#clr.beta = log(beta/gm(beta)), and so contrasts between clr.beta entries are contrasts between log(beta) entries.

beta <- exp(summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilr1+ilr2+ilr3+ilr4+ilr5+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)$estimate[2:6]%*%Phi6)#/sum(exp(summary(pool(with(CSS_track_mids,glm((bmi>=25) ~ (ilr1+ilr2+ilr3+ilr4+ilr5+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)$estimate[2:6]%*%Phi6))
#Now we can get contrasts between beta's in some simplex - does the particular version of the simplex matter?

## How should we present the contrast parameters?

summary(pool(with(CSS_track_mids,glm((bmi>=30) ~ (ilr1+ilr2+ilr3+ilr4+ilr5+selfScoreCat+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)

#m <- gamlss(bmi~(cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+selfScoreCat+age+gender+education+occupation), sigma.formula = ~ 1,
#            nu.formula = ~ 1,weights=sample_weights, data=na.omit(subset(CSS_track[,c("cluster1prob","cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","selfScoreCat","age","gender","education","occupation","bmi","sample_weights","imputation")],imputation==10)),family=BCCG,method=RS(100))
#robust=TRUE? Doesn't change much when the fit is this good.

#Numeric BMI

#Confidence intervals:


# --------------------------------------------------------------------------- ##
# --------------------------------------------------------------------------- ##

#Tracking data: Population sample (random sample) - same analysis

hist(pop_track$bmi,breaks=50,xlim=c(0,50))

ggplot(pop_track, aes(x = factor(selfScoreCat))) +
  geom_bar()
ggplot(pop_track, aes(x = factor(cluster))) +
  geom_bar()

#analyses

## regression analysis of clusters of night-time smartphone use and overweight/obesity #justeres for selfScoreCat?? ## hvad er de forskellige clusters?? ## fortolkning?? ## inkluderer imp_nr=0?

## Continuous Outcome

#With gamlss:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+selfScoreCat+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG) #Alternative is BCT, but fits with many many degrees of freedom, suggesting no real difference between BCCG and BCT.
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}
plot(m)

pool_inf_PopTrack <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrack$qbar
pool_inf_PopTrack$ubar
pool_inf_PopTrack$pval

summary(pool_inf_PopTrack)

#Profile intervals: 
lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[2,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[2,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[2,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[3,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[3,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[3,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[4,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[4,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[4,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[5,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[5,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[5,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[6,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[6,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[6,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[7,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[7,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[7,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[8,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[8,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[8,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[9,5]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[9,1]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrack)[9,6]+10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrack <- cbind(c(lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6,lowerCat2,lowerCat3,lowerCat4),
                           c(estClust2,estClust3,estClust4,estClust5,estClust6,estCat2,estCat3,estCat4),
                           c(upperClust2,upperClust3,upperClust4,upperClust5,upperClust6,upperCat2,upperCat3,upperCat4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrack$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrack$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
rownames(confints_PopTrack) <- names(pool_inf_PopTrack$qbar)[2:9]

#Trends:


#SelfScore
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+as.numeric(selfScoreCat)+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackSTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrackSTrend$qbar
pool_inf_PopTrackSTrend$ubar
pool_inf_PopTrackSTrend$pval

summary(pool_inf_PopTrackSTrend) #Get the p value for trend from here, as sigma and nu are constant.

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[3,5]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[3,1]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[3,6]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[4,5]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[4,1]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[4,6]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[5,5]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[5,1]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[5,6]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[6,5]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[6,1]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[6,6]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[7,5]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[7,1]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackSTrend)[7,6]+10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackSTrend <- cbind(c(lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6,lowerCat),
                           c(estClust2,estClust3,estClust4,estClust5,estClust6,estCat),
                           c(upperClust2,upperClust3,upperClust4,upperClust5,upperClust6,upperCat))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

#Track
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(track_severity+selfScoreCat+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation","track_severity")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackTTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrackTTrend$qbar
pool_inf_PopTrackTTrend$ubar
pool_inf_PopTrackTTrend$pval #Get the p value for trend from here, as sigma and nu are constant.

summary(pool_inf_PopTrackTTrend)#Get the p value for trend from here, as sigma and nu are constant.

lowerTrack <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estTrack <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperTrack <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[3,5]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[3,1]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[3,6]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[4,5]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[4,1]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[4,6]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[5,5]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[5,1]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackTTrend)[5,6]+10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackTTrend <- cbind(c(lowerTrack,lowerCat2,lowerCat3,lowerCat4),
                                 c(estTrack,estCat2,estCat3,estCat4),
                                 c(upperTrack,upperCat2,upperCat3,upperCat4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 






#Without adjustment for tracking:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(selfScoreCat+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoT <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrackNoT$qbar
pool_inf_PopTrackNoT$ubar
pool_inf_PopTrackNoT$pval

summary(pool_inf_PopTrackNoT)

lowerCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[2,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[3,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,5]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,1]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoT)[4,6]+10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackNoT <- cbind(c(lowerCat2,lowerCat3,lowerCat4),
                                 c(estCat2,estCat3,estCat4),
                                 c(upperCat2,upperCat3,upperCat4))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoT$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


#Trend:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(as.numeric(selfScoreCat)+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoTTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrackNoTTrend$qbar
pool_inf_PopTrackNoTTrend$ubar
pool_inf_PopTrackNoTTrend$pval

summary(pool_inf_PopTrackNoTTrend)

lowerCat <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estCat <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperCat <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoTTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoTTrend <- c(lowerCat,estCat,upperCat) -  integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoTTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 



#Without adjustment for selfScore:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoS <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrackNoS$qbar
pool_inf_PopTrackNoS$ubar
pool_inf_PopTrackNoS$pval

summary(pool_inf_PopTrackNoS)

lowerClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust2 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust2 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[2,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[3,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust3 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[3,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust3 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[3,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[4,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust4 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[4,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust4 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[4,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[5,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust5 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[5,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust5 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[5,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

lowerClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[6,5]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estClust6 <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[6,1]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperClust6 <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoS)[6,6]+10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


confints_PopTrackNoS <- cbind(c(lowerClust2,lowerClust3,lowerClust4,lowerClust5,lowerClust6),
                                 c(estClust2,estClust3,estClust4,estClust5,estClust6),
                                 c(upperClust2,upperClust3,upperClust4,upperClust5,upperClust6))-integrate(function(y) y*dBCCG(x=y,mu=10,sigma=exp(pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoS$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 


#Trend:
coefs <- list()
ses <- list()
vcovs <- list()

for (i in 1:N_imp){
  m <- gamlss(bmi~(track_severity+age+Gender+education+occupation), sigma.formula = ~1, nu.formula =~ 1, weights=sample_weights, data=na.omit(subset(pop_track[,c("cluster2prob","cluster3prob","cluster4prob","cluster5prob","cluster6prob","bmi","selfScoreCat","age","Gender","education","occupation","sample_weights","imputation","track_severity")],imputation==i)),family = BCCG)
  m_sum <- summary(m)
  coefs[[i]] <- m_sum[,1]
  ses[[i]] <- m_sum[,2]
  vcovs[[i]] <- vcov(m)
}

pool_inf_PopTrackNoSTrend <- miceadds::pool_mi(qhat = coefs, u = vcovs)
pool_inf_PopTrackNoSTrend$qbar
pool_inf_PopTrackNoSTrend$ubar
pool_inf_PopTrackNoSTrend$pval

summary(pool_inf_PopTrackNoSTrend)

lowerTrack <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoSTrend)[2,5]+10,sigma=exp(pool_inf_PopTrackNoSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
estTrack <-  integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoSTrend)[2,1]+10,sigma=exp(pool_inf_PopTrackNoSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 
upperTrack <- integrate(function(y) y*dBCCG(x=y,mu=summary(pool_inf_PopTrackNoSTrend)[2,6]+10,sigma=exp(pool_inf_PopTrackNoSTrend$qbar[length(m$mu.coefficients)+1]),nu=pool_inf_PopTrackNoSTrend$qbar[length(m$mu.coefficients)+length(m$sigma.coefficients)+1]),0,Inf)$value 

confints_PopTrackNoSTrend <- c(lowerTrack,estTrack,upperTrack)


# --------------------------------------------------------------------------- ##
### Binary Outcomes for population sample
# --------------------------------------------------------------------------- ##

## BMI > 25

## Clusters and BMI<25 
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25No <- with(pop_track_mids,glm((bmi>=25) ~ (cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25No <- summary(pool(Random25No), conf.int=T)
cbind(exp(modelRandom25No$estimate),
exp(modelRandom25No$`2.5 %`),
exp(modelRandom25No$`97.5 %`))

#test for trend  (clusters as numeric)
Random25NoSest <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(track_severity)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25NoSest), conf.int=T)

## self-reported smartphone use

## night-time smartphone use and BMI >25
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25Night <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25Night <- summary(pool(Random25Night), conf.int=T)
cbind(exp(modelRandom25Night$estimate),
      exp(modelRandom25Night$`2.5 %`),
      exp(modelRandom25Night$`97.5 %`))

## test for trend night-time use and BMI >25
Random25NightTest <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseNight)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25NightTest), conf.int=T)

## smartphone use before sleep onset and BMI >25
summary(pool(with(pop_track_mids,glm((bmi>=25) ~ (mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random25Before <- with(pop_track_mids,glm((bmi>=25) ~ (mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom25Before <- summary(pool(Random25Before), conf.int=T)
cbind(exp(modelRandom25Before$estimate),
      exp(modelRandom25Before$`2.5 %`),
      exp(modelRandom25Before$`97.5 %`))

## test for trend smartphone use before sleep and BMI >25
Random25BeforeTest <- with(pop_track_mids,glm((bmi>=25) ~ (as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random25BeforeTest), conf.int=T)

# --------------------------------------------------------------------------- ##
## BMI > 30

## clusters and BMI>30 
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30No <- with(pop_track_mids,glm((bmi>=30) ~ (cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30No <- summary(pool(Random30No), conf.int = T)
cbind(exp(modelRandom30No$estimate),
exp(modelRandom30No$`2.5 %`),
exp(modelRandom30No$`97.5 %`))

## night-time smartphone use and BMI >30
summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30Night <- with(pop_track_mids,glm((bmi>=30) ~ (mobileUseNight+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30Night <- summary(pool(Random30Night), conf.int=T)
cbind(exp(modelRandom30Night$estimate),
      exp(modelRandom30Night$`2.5 %`),
      exp(modelRandom30Night$`97.5 %`))

#test for trend (selfScoreCat as numeric and no adjustment for clusters)
Random30Night <- with(pop_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseNight)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random30Night), conf.int = T)

## smartphone use before sleep onset and BMI >30

summary(pool(with(pop_track_mids,glm((bmi>=30) ~ (mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))),conf.int=T)
Random30Before <- with(pop_track_mids,glm((bmi>=30) ~ (mobileUseBeforeSleep+age+gender+education+occupation), weights=sample_weights,family=binomial))
modelRandom30Before <- summary(pool(Random30Before), conf.int=T)
cbind(exp(modelRandom30Before$estimate),
exp(modelRandom30Before$`2.5 %`),
exp(modelRandom30Before$`97.5 %`))

#test for trend (selfScoreCat as numeric and no adjustment for clusters)
Random30Before <- with(pop_track_mids,glm((bmi>=30) ~ (as.numeric(mobileUseBeforeSleep)+age+gender+education+occupation), weights=sample_weights,family=binomial))
summary(pool(Random30Before), conf.int = T)


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

## night-time smartphone use
table(clinical_sample$mobileUseNight)
clinical_sample$mobileUseNight <- factor(clinical_sample$mobileUseNight, levels = c("Never", "A few times a month or less", "A few times a week", "Every night or almost every night"))

## smartphone use before sleep
table(clinical_sample$mobileUseBeforeSleep)
clinical_sample$mobileUseBeforeSleep <- factor(clinical_sample$mobileUseBeforeSleep, levels = c("Every month or less", "Once a week", "2-4 times per week", "5-7 times per week"))

# --------------------------------------------------------------------------- ##
## descriptive of clinical sample
## age
publish(univariateTable(mobileUseNight ~ age.x,data=clinical_sample, column.percent=TRUE))/25
clinical_sample$age.x <- as.numeric(clinical_sample$age.x)

## BMI
publish(univariateTable(mobileUseNight ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))

#The subjects are scoring in the high end. Is this an issue or a characteristic of the data?

#Introducing interesting derived variables

clinical_sample$bmi <- as.numeric(clinical_sample$bmi.clinical)
clinical_sample$bmi25 <- as.numeric(clinical_sample$bmi.clinical>=25)
clinical_sample$bmi30 <- as.numeric(clinical_sample$bmi.clinical>=30)

publish(univariateTable(mobileUseNight ~ bmi25,data=clinical_sample, column.percent=TRUE))

## age at clinical examination 
table(clinical_sample$age.y, useNA="always")
clinical_sample$age<- as.numeric(str_c(substr(clinical_sample$age.y,1,1),substr(clinical_sample$age.y,2+(mod(nchar(clinical_sample$age.y),4)==1),2+(mod(nchar(clinical_sample$age.y),4)==1)),".",
                 substr(clinical_sample$age.y,3+(mod(nchar(clinical_sample$age.y),4)!=3),3+(mod(nchar(clinical_sample$age.y),4)!=3))))


# --------------------------------------------------------------------------- ##

## systolic blood pressure
clinical_sample$sbp<-rowMeans(cbind(clinical_sample$sbp1,clinical_sample$sbp2,clinical_sample$sbp3),na.rm=T)
publish(univariateTable(mobileUseNight ~ sbp,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ sbp,data=clinical_sample, column.percent=TRUE))

# diastolic blood pressure
clinical_sample$dbp<-rowMeans(cbind(clinical_sample$dbp1,clinical_sample$dbp2,clinical_sample$dbp3),na.rm=T)
publish(univariateTable(mobileUseNight ~ dbp,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ dbp,data=clinical_sample, column.percent=TRUE))

## hip waist ratio

clinical_sample$ratiowaisthip <- as.numeric(clinical_sample$ratiowaisthip)
publish(univariateTable(mobileUseNight ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ ratiowaisthip,data=clinical_sample, column.percent=TRUE))

## bmi clinical
clinical_sample$bmi.clinical <- as.numeric(clinical_sample$bmi.clinical)
publish(univariateTable(mobileUseNight ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ bmi.clinical,data=clinical_sample, column.percent=TRUE))


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
publish(univariateTable(mobileUseNight ~ hdl,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ hdl,data=clinical_sample, column.percent=TRUE))


## LDL
clinical_sample$ldl <- as.numeric(clinical_sample$ldl)
publish(univariateTable(mobileUseNight~ ldl,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ ldl,data=clinical_sample, column.percent=TRUE))

## VLDL
clinical_sample$vldl <- as.numeric(clinical_sample$vldl)
publish(univariateTable(mobileUseNight~ vldl,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ vldl,data=clinical_sample, column.percent=TRUE))

## total cholesterol
clinical_sample$t_cholesterol <- as.numeric(clinical_sample$t_cholesterol)
publish(univariateTable(mobileUseNight ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ t_cholesterol,data=clinical_sample, column.percent=TRUE))

## triglycerides
clinical_sample$triglycerids <- as.numeric(clinical_sample$triglycerids)
publish(univariateTable(mobileUseNight ~ triglycerids,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ triglycerids,data=clinical_sample, column.percent=TRUE))

## hba1c
clinical_sample$hba1c <- as.numeric(clinical_sample$hba1c)
publish(univariateTable(mobileUseNight ~ hba1c,data=clinical_sample, column.percent=TRUE))
publish(univariateTable(mobileUseBeforeSleep ~ hba1c,data=clinical_sample, column.percent=TRUE))

# --------------------------------------------------------------------------- ##

#Models - multiple testing issue if we are going to 'pick and choose' which responses we would like to look at.
table(clinical_sample$age.x)

#Transforming to mids for modelling and inference

clinical_mids <- as.mids(clinical_sample,.imp="imputation",.id="userid")

#hdl
#hdl_sum<-cbind(summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(clinical_mids,lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
lines(seq(from=min(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),length.out=100),dnorm(x=seq(from=min(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),max(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),na.rm=T),length.out=100),mean=mean(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))),sd=sd(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))))
plot(residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))


#ldl 
#ldl_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res=residuals(lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))
plot(fitted(lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

cbind(confint(glm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#vldl
#vldl_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),
      cbind(coef(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))-1.96*summary(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2],
            coef(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))+1.96*summary(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2]))


#t_cholesterol
#t_cholesterol_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#                summary(pool(with(data=clinical_mids, glm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(lm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(lm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#triglycerids
#tri_sum<-cbind(summary(pool(with(data=clinical_mids, glm(as.numeric(triglycerids) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, glm(as.numeric(triglycerids) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(glm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),breaks=20,prob=T)
res <- residuals(glm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(glm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))
plot(fitted(glm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),residuals(glm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)))

cbind(confint(glm(as.numeric(triglycerids) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit)),
      cbind(coef(glm(as.numeric(triglycerids) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))-1.96*summary(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2],
            coef(glm(as.numeric(triglycerids) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))+1.96*summary(glm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1), na.action=na.omit))$coefficients[,2]))

#hba1c
#hba1c_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),breaks=20,prob=T)
res <- residuals(lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)),residuals(lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1),na.action=na.omit)))

# bmi
hist(residuals(lm(as.numeric(bmi.clinical) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(bmi.clinical) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(bmi.clinical) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(bmi.clinical) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))


#ratiowaisthip
#wh_sum<-cbind(summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation))))$std.error[2:4])
hist(residuals(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#sbp
#sbp_sum<-cbind(summary(pool(with(data=clinical_mids, lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(sbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(sbp) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(sbp) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#dbp
#dbp_sum<-cbind(summary(pool(with(data=clinical_mids, lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$estimate[2:4],
#      summary(pool(with(data=clinical_mids, lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))))$std.error[2:4])
hist(residuals(lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),breaks=20,prob=T)
res <- residuals(lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1)))
res_seq=seq(from=min(res),to=max(res),length.out=100)
lines(res_seq,dnorm(res_seq,mean=mean(res),sd=sd(res)))
plot(residuals(lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))
plot(fitted(lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))),residuals(lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation, data=subset(clinical_sample,imputation==1))))

cbind(confint(glm(as.numeric(dbp) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1))),
      confint(lm(as.numeric(dbp) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age.x+education+occupation,na.action=na.omit,data=subset(clinical_sample,imputation==1)),type="Wald"))

#Generally the residuals look reasonably centered, with a few positive outliers. The residual distributions on the first imputation actually look reasonably normal, save for the few (extreme) outliers.

#Seems that these models are appropriate, and that normal approximations of confidence interval will be reasonable too. We can however also just use the profile likelihood CI's.

#Wald confidence intervals

dbp_int <- summary(pool(with(data=clinical_mids, lm(dbp ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_int <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_int <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ cluster2prob+cluster3prob+cluster4prob+cluster5prob+cluster6prob+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])


df_ints <- data.frame(rbind(hdl_int[2,],ldl_int[2,],vldl_int[2,],t_chol_int[2,],tri_int[2,],hba1c_int[2,],dbp_int[2,],sbp_int[2,],wh_int[2,],glu_int[2,],bmi_int[2,]),
                            rbind(hdl_int[3,],ldl_int[3,],vldl_int[3,],t_chol_int[3,],tri_int[3,],hba1c_int[3,],dbp_int[3,],sbp_int[3,],wh_int[3,],glu_int[3,],bmi_int[3,]),
                            rbind(hdl_int[4,],ldl_int[4,],vldl_int[4,],t_chol_int[4,],tri_int[4,],hba1c_int[4,],dbp_int[4,],sbp_int[4,],wh_int[4,],glu_int[4,],bmi_int[4,]),
                            rbind(hdl_int[5,],ldl_int[5,],vldl_int[5,],t_chol_int[5,],tri_int[5,],hba1c_int[5,],dbp_int[5,],sbp_int[5,],wh_int[5,],glu_int[5,],bmi_int[5,]),
                            rbind(hdl_int[6,],ldl_int[6,],vldl_int[6,],t_chol_int[6,],tri_int[6,],hba1c_int[6,],dbp_int[6,],sbp_int[6,],wh_int[6,],glu_int[6,],bmi_int[6,]))

colnames(df_ints) <- c("type.2.estimate","type.2.lower","type.2.upper","type.2.pvalue","type.3.estimate","type.3.lower","type.3.upper","type.3.pvalue","type.4.estimate","type.4.lower","type.4.upper","type.4.pvalue",
                             "type.5.estimate","type.5.lower","type.5.upper","type.5.pvalue","type.6.estimate","type.6.lower","type.6.upper","type.6.pvalue")
rownames(df_ints) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_ints

## Night-time smartphone use and biomarkers

dbp_intS <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intS <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseNight+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_intsS <- data.frame(rbind(hdl_intS[2,],ldl_intS[2,],vldl_intS[2,],t_chol_intS[2,],tri_intS[2,],hba1c_intS[2,],dbp_intS[2,],sbp_intS[2,],wh_intS[2,],glu_intS[2,],bmi_intS[2,]),
                      rbind(hdl_intS[3,],ldl_intS[3,],vldl_intS[3,],t_chol_intS[3,],tri_intS[3,],hba1c_intS[3,],dbp_intS[3,],sbp_intS[3,],wh_intS[3,],glu_intS[3,],bmi_intS[3,]),
                      rbind(hdl_intS[4,],ldl_intS[4,],vldl_intS[4,],t_chol_intS[4,],tri_intS[4,],hba1c_intS[4,],dbp_intS[4,],sbp_intS[4,],wh_intS[4,],glu_intS[4,],bmi_intS[4,]))

colnames(df_intsS) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsS) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsS


## Smartphone use before sleep onset and biomarkers
dbp_intS <- summary(pool(with(data=clinical_mids, lm(dbp ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(dbp_sum[,1]-1.96*dbp_sum[,2],dbp_sum[,1]+1.96*dbp_sum[,2])
glu_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(glucose) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(glu_sum[,1]-1.96*glu_sum[,2],glu_sum[,1]+1.96*glu_sum[,2])
hba1c_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hba1c) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hba1c_sum[,1]-1.96*hba1c_sum[,2],hba1c_sum[,1]+1.96*hba1c_sum[,2])
hdl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(hdl) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(hdl_sum[,1]-1.96*hdl_sum[,2],hdl_sum[,1]+1.96*hdl_sum[,2])
ldl_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ldl) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(ldl_sum[,1]-1.96*ldl_sum[,2],ldl_sum[,1]+1.96*ldl_sum[,2])
t_chol_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(t_cholesterol) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(t_cholesterol_sum[,1]-1.96*t_cholesterol_sum[,2],t_cholesterol_sum[,1]+1.96*t_cholesterol_sum[,2])
sbp_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(sbp) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(sbp_sum[,1]-1.96*sbp_sum[,2],sbp_sum[,1]+1.96*sbp_sum[,2])
tri_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(triglycerids) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(tri_sum[,1]-1.96*tri_sum[,2],tri_sum[,1]+1.96*tri_sum[,2])
vldl_intS <- summary(pool(with(data=clinical_mids,lm(as.numeric(vldl) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit,family=Gamma))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(vldl_sum[,1]-1.96*vldl_sum[,2],vldl_sum[,1]+1.96*vldl_sum[,2])
wh_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(ratiowaisthip) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])
bmi_intS <- summary(pool(with(data=clinical_mids, lm(as.numeric(bmi.clinical) ~ mobileUseBeforeSleep+age+gender+education+occupation,na.action=na.omit))),conf.int=T)[,c("estimate","2.5 %", "97.5 %","p.value")]#cbind(wh_sum[,1]-1.96*wh_sum[,2],wh_sum[,1]+1.96*wh_sum[,2])

df_intsS <- data.frame(rbind(hdl_intS[2,],ldl_intS[2,],vldl_intS[2,],t_chol_intS[2,],tri_intS[2,],hba1c_intS[2,],dbp_intS[2,],sbp_intS[2,],wh_intS[2,],glu_intS[2,],bmi_intS[2,]),
                       rbind(hdl_intS[3,],ldl_intS[3,],vldl_intS[3,],t_chol_intS[3,],tri_intS[3,],hba1c_intS[3,],dbp_intS[3,],sbp_intS[3,],wh_intS[3,],glu_intS[3,],bmi_intS[3,]),
                       rbind(hdl_intS[4,],ldl_intS[4,],vldl_intS[4,],t_chol_intS[4,],tri_intS[4,],hba1c_intS[4,],dbp_intS[4,],sbp_intS[4,],wh_intS[4,],glu_intS[4,],bmi_intS[4,]))

colnames(df_intsS) <- c("cat.2.estimate","cat.2.lower","cat.2.upper","cat.2.pvalue","cat.3.estimate","cat.3.lower","cat.3.upper","cat.3.pvalue","cat.4.estimate","cat.4.lower","cat.4.upper","cat.4.pvalue")
rownames(df_intsS) <- c("hdl","ldl","vldl","total cholesterol","triglycerids","hba1c","dbp","sbp","waist-hip-ratio","glucose","bmi")
df_intsS

