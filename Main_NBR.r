###### Main program to illustrate how to perform nbreg using a simulated censored observational data ###
# tested under R version 4.0.4

# if the two packages are not installed, install them first
# install.packages("geepack")
# install.packages("dplyr")

require(geepack)
require(dplyr)

#set the path to the folder with Est_NBR.r and dataset file (use "/" instead of "\" in path)
setwd("Your Path")
source("Est_NBR.r")	# load functions used to estimate

############### Read in data ###########

data=read.csv("Censored_CEdata.csv")

head(data)	#show first few rows of data

# make the 3 categorical variables in covariates as factor before running analysis 
# not neccesary since they are 0 or 1 binary, but useful for 3+ level variables
data=data %>% mutate(across(c(Age65,LBBB,Female), as.factor))

# check cost history data
head(data[,9:23])
head(data[,9:24])

# QALY history data
head(data[,25:39])

# check covariates
head(data[,6:8])

# divide all costs by 1000 so that cost is in $1000
data[,9:24]=data[,9:24]/1000
head(data[,9:24])	#check cost again

table(data$Trt)
#   0   1 
# 303 297 
# 297 patients in considered treatment group, 303 in comparison 

#check distribution of follow-up time (before truncted by L)
summary(data$survival)
hist(data$survival)	#check histogram


################################# Fit net-benefit reg (NBR) ########################

# Categorical covariates with numeric values need to set as factor before running analysis
# Covariates in interactions between covariates and treatment must be a subset of main effects of covariates
# Patients with missing values will be excluded by nbreg function automatically

#lambda=c(0,2,5,6)	# can select a few WTP values $0, $2,000, $5,000, $6,000
lambda=seq(0,6,0.5)	# choose WTP as a sequence from $0 to $6,000 for 1 additional QALY

########## An example using SW method ##########
##  fit covariate-adjusted NBR using Simple Weighted (SW) method, effectiveness is quality-adjusted life years (QALY)
# use cost and effectiveness history to better truncate them by L
fit1<-nbreg(Followup=data$survival,		#follow-up time	
	delta=data$dead,				#death indicator (1=dead,0=censor)
	group=data$Trt,				#treatment indicator (1=considered treatment, 0=comparison)
	Cost=data[,9:23],				#a matrix of cost history (recommended to better truncate cost by L)
	Eff=data[,25:39],				#a matrix of effectiveness history (recommended to better truncate by L),default is survival if not provided
	Patition.times=1:15,			#end timepoints of intervals, meaning that intervals are [0,1],(1,2],...,(14,15] years for grouped costs
	Method='SW',				#use Simple Weighted Method 
	Z=data[,6:8],				#adjust 3 covariates in regressions (main effects), will do unadjusted regression if not provided
	Eff.only=TRUE,				#also fit a effectiveness only regression (WTP lambda=NA, not belong to net-benefit regression)
	Cost.only=TRUE,				#also fit a Cost only regression (equivalent to setting WTP lambda=0 and then switching sign of coef est)
	lambda=lambda,				#WTP for 1 additional QALY
	L=10						#time limit within L years
)
# main results printed out, detailed results saved in fit1

## can simply fit cost-only and effect-only regressions to quickly calculate adjusted ICER, with lambda=NULL option
fit1_0<-nbreg(Followup=data$survival, delta=data$dead, group=data$Trt, Cost=data[,9:23],
	Eff=data[,25:39], Patition.times=1:15, Method='SW', Z=data[,6:8],	Eff.only=TRUE,Cost.only=TRUE,	L=10,
	lambda=NULL 		#equivalently, can remove this option						
)

## special case of no covariates (eg, for randomized studies), by removing the option Z= 
fit1_1<-nbreg(Followup=data$survival, delta=data$dead, group=data$Trt, Cost=data[,9:23],
	Eff=data[,25:39], Patition.times=1:15, Method='SW', Eff.only=TRUE,Cost.only=TRUE,	L=10,
	lambda=lambda						
)

# obtain covariate-adjusted ICER based on SW method
# need to change sign for cost-only regression (lambda=0) because dependent variable = -Cost
2.54654804/0.9482411 
#[1] 2.685549


## use total costs/eff in 15 years, not very accurate to truncate them by L=10 (nbreg simply prorate them)
# hence providing cost/eff histories is more accurate
fit2<-nbreg(Followup=data$survival, delta=data$dead,group=data$Trt,		
	Cost=data$tot.cost,			# observed total cost (recommend to use cost history to better truncate cost in L)
	Eff=data$tot.QALY,			# observed total effectiveness, default is survival if not provided
	Method='SW', Z=data[,6:8], Eff.only=TRUE, lambda=lambda, L=10)


## fit covariate-adjusted NBR using SW method, effectiveness is life years
# does not need to provide Eff, Followup will be used to calculate effectiveness
fit3<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],
	Eff=NULL,	#equivalently, this option can be removed
	Patition.times=1:15,Method='SW',Z=data[,6:8], Eff.only=TRUE,lambda=lambda,L=10)

## too large L=15 leads to error, need to choose L where a relatively large number of patients are still under observation
fit3_bigL<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Patition.times=1:15,
	Method='SW',Z=data[,6:8], Eff.only=TRUE,lambda=lambda,L=15)

## fit covariate-adjusted NBR using Partition (PT) method with QALY as Eff, which leads to smaller SE than SW method
fit4<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],
	Patition.times=1:15,Method='PT',Z=data[,6:8], Eff.only=TRUE,lambda=lambda,L=10)

# compare with naive Complete-case only (CC) method 
fit4_cc<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='CC',Z=data[,6:8],Eff.only=TRUE,lambda=lambda,L=10)

# compare with naive All data ignoring censoring (AL) method 
fit4_al<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='AL',Z=data[,6:8],Eff.only=TRUE,lambda=lambda,L=10)


########## code example for dataset with unequal time intervals ############# 
# assume time intervals in dataset do not have same length (although not true for this dataset)
# assume cost.1 (QALY.1) is cost/QALY in first 2 years, cost.2 (QALY.2) and cost.3 (QALY.3) are cost/QALY in the following 6 months (other time intervals keep the same)
# thus, intervals are [0,2],(2,2.5],(2.5,3],(3,4],...,(14,15] for cost and QALY histories
fit4_unequal<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],
	Patition.times=c(2,2.5,3:15),		# set the end times of first 3 intervals as 2, 2.5, and 3, and following end times are increasing by 1 until 15
	Method='PT',Z=data[,6:8], Eff.only=TRUE,lambda=lambda,L=10)

########## Possible issue in variable name when only one covariate is provided #######
## one covariate only, may drop covariate name if input covariate name is not correctly kept in R 
fit5_1<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='PT',
	Z=data[,6],					#1 covariate in regressions may drop the covariate name by R automatically and show the default name Z 
	Eff.only=TRUE, lambda=lambda,L=10)

fit5_2<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='PT',
	Z=data.frame(Age65=data[,6]),					#can fix it by re-assign the name
	Eff.only=TRUE, lambda=lambda,L=10)

fit5_3<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='PT',
	Z=data[,6,drop=FALSE],						#or fix it by using drop=FALSE option
	Eff.only=TRUE, lambda=lambda,L=10)


######### fit covariate-adjusted NBR with interactions between treatment and covariates  #######
fit6<-nbreg(Followup=data$survival, delta=data$dead, group=data$Trt, Cost=data[,9:23], Eff=data[,25:39], Patition.times=1:15, 
	Method='PT',Z=data[,6:8],
	interaction=c("LBBB"),			#Treatment-LBBB interaction, must belong main effects Z, can add more, eg, interaction=c("Age65","LBBB"), or interaction=names(data[,6:8])
	Eff.only=TRUE, lambda=lambda, L=10)

# include both Treatment x LBBB and Treatment x Female interactions 
fit6_2<-nbreg(Followup=data$survival, delta=data$dead, group=data$Trt, Cost=data[,9:23], Eff=data[,25:39], Patition.times=1:15, 
	Method='PT',Z=data[,6:8],
	interaction=c("LBBB","Female"),
	Eff.only=TRUE, lambda=lambda, L=10)

# include interactions between treatment and all 3 covariate
fit6_3<-nbreg(Followup=data$survival, delta=data$dead, group=data$Trt, Cost=data[,9:23], Eff=data[,25:39], Patition.times=1:15, 
	Method='PT',Z=data[,6:8],
	interaction=names(data[,6:8]),
	Eff.only=TRUE, lambda=lambda, L=10)

############# Doubly Robust method ############
fit7<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='PT',Z=data[,6:8],interaction=c("LBBB"),
	PS.Z=data[,6:8],			#all covariates are used to estimate propensity scores using logistic regression, will fit unadjusted logistic if not provided
	Doubly.Robust=TRUE,		#do DR method to estimate causal average incremental net benefit (INB)
	Eff.only=TRUE,lambda=lambda,L=10)

#view saved fitted models for 1st WTP value, can change to fit7[[2]] to look at results for 2nd WTP
#last one or two are WTP=NA if Eff.only=TRUE and/or Cost.only =TRUE
fit7[[1]]$lambda		#WTP of the 1st one
fit7[[1]]$Reg.type		#Regression type of the 1st one (Effect, Cost, or net-benefit reg (NBR))
fit7[[1]]$coef.table
fit7[[1]]$Regmodel
fit7[[1]]$PSmodel

# examine propensity scores
summary(fit7[[1]]$PS)

# histograms for propensity scores
require(ggplot2)	
ggplot(data,aes(x=fit7[[1]]$PS))+geom_histogram()+facet_grid(rows = fit7[[1]]$group)+ 
	xlab("Estimated Propensity Score")+xlim(c(0,1))


# DR method with regressions including all possible treatment-covariate interactions, equivalent to fitting regressions separately in each treatment group  
fit8<-nbreg(Followup=data$survival,delta=data$dead,group=data$Trt,Cost=data[,9:23],Eff=data[,25:39],Patition.times=1:15,
	Method='PT',Z=data[,6:8],
	interaction=names(data[,6:8]),	#include all treatment-covariate interactions 
	PS.Z=data[,6:8], Doubly.Robust=TRUE, Eff.only=TRUE,lambda=lambda,L=10)


############### Compare results using the true uncensored data, 4 methods are equivalent to OLS when no censoring ###########

data.uncensored=read.csv("True_CEdata.csv")
head(data.uncensored)		#show first 5 rows of data
data.uncensored=data.uncensored %>% mutate(across(c(Age65,LBBB,Female), as.factor))
data.uncensored[,9:24]=data.uncensored[,9:24]/1000	# divide all costs by 1000 so that cost is in $1000

# covariate-adjusted regression without interaction
fit9<-nbreg(Followup=data.uncensored$survival,delta=data.uncensored$dead,group=data.uncensored$Trt,
	Cost=data.uncensored[,9:23],Eff=data.uncensored[,25:39],Patition.times=1:15,
	Method='CC',Z=data.uncensored[,6:8],Eff.only=TRUE,lambda=lambda,L=10)

# covariate-adjusted regression with interaction
fit10<-nbreg(Followup=data.uncensored$survival,delta=data.uncensored$dead,group=data.uncensored$Trt,
	Cost=data.uncensored[,9:23],Eff=data.uncensored[,25:39],Patition.times=1:15,
	Method='CC',Z=data.uncensored[,6:8],interaction=c("LBBB"),Eff.only=TRUE,lambda=lambda,L=10)

# DR method
fit11<-nbreg(Followup=data.uncensored$survival,delta=data.uncensored$dead,group=data.uncensored$Trt,
	Cost=data.uncensored[,9:23],Eff=data.uncensored[,25:39],Patition.times=1:15,
	Method='CC',Z=data.uncensored[,6:8],interaction=c("LBBB"),PS.Z=data.uncensored[,6:8], Doubly.Robust=TRUE,Eff.only=TRUE,lambda=lambda,L=10)

############ CEAC plot ################


# For non-DR method, if no interaction with treatment, there is only 1 CEAC curve 
# meaning that cost-effectiveness is the same across all subgroups
plot(fit4,ylab="Probability new treatment is cost-effective",xlab="WTP (in $1000) for one additional QALY",lwd=2, pch=19,cex=1.2)	#PT adjusted

# if there is interaction between one or more covariates and treatment,
# cost-effectiveness is heterougenous across subgroups defined by the covariate(s) in interaction, each subgroup has 1 curve 
# If there are 2 covariates in interactions, can provide more subgroup values, eg, "subgroup=list(LBBB=0, Age65=1)"
# add=TRUE to add curve to existing figure instead of creating a new one
plot(fit6,subgroup=list(LBBB=0),add=TRUE,col="gray50",lwd=2,lty=2,pch=15,cex=1.2)		#subgroup of non-LBBB, adjusted for age and gender
plot(fit6,subgroup=list(LBBB=1),add=TRUE,col="gray50",lwd=2,lty=3,pch=17,cex=1.2)		#subgroup of LBBB, adjusted for age and gender

# for DR method 
plot(fit7,add=TRUE,lty=4,col="gray70",lwd=2,pch=0,cex=1.2)

# add legend to plot
legend('right',c("Adjusted","Subgroup:non-LBBB","Subgroup:LBBB","Doubly Robust"),
	lty=c(1,2,3,4),lwd=c(2,2,2,2),col=c(1,"gray50","gray50","gray70"),pch=c(19,15,17,0),cex=c(1.2,1.2,1.2,1.2),seg.len =2.8)

