###### Main program to create a simulated censored observational data ###

# if the two packages are not installed, install them first
# install.packages("mvtnorm")
# install.packages("survival")

require(mvtnorm)
require(survival)

#set the path to the folder with the following R file (use "/" instead of "\" in path)
setwd("You Path")
source("Data_Gen.r")	# load functions used to generate data


####### simulate a dataset (in year) with grouped costs and LY/QALY, and 3 binary covariates (age>=65, LBBB, female) ###########

set.seed(123)	#set random seed

n=2000	#sample size

data<-sim.CostEff.Obs(n=n,		#sample size
	Coef.Surv.MainEff=c(2.5,-0.2,0.2,0.1), #coef for survival model, 1st for intercept, others for main effects of covariates	
	Coef.Surv.InterEff=c(0.05,0,0.6,0),  #coef for survival model, 1st for main effect of treatment, others for treatment-covariate interactions
	Par.Censor=15,			#Parameter for max censoring time=15 years
	Coef.Trt=c(-0.5,-0.5,1.5,0),		#coef for treatment assignment model using logistic regression
	p=3,					#3 covariates
	k=1,					#k=1-year interval to group outcomes (ie, costs within (0,1], costs within (1,2], etc.)
	max.time=15,					# choose 15 years as max time 
	cov.name=c("Age65","LBBB","Female")		#if not provided, covariates will be X1, X2, etc. 
)

# save the observed censored data
datacsv=data.frame(id=1:n,survival=data$Follow.up,dead=data$Death,Trt=data$Trt,data$Covariate,
	cost=data$Cost.grp.obs,tot.cost=apply(data$Cost.grp.obs,1,sum),
	QALY=data$QASurv.grp.obs,tot.QALY=apply(data$QASurv.grp.obs,1,sum)
	)
write.csv(datacsv,"Censored_CEdata.csv",row.names=FALSE)

# save true uncensored data for comparison
datacsv.true=data.frame(id=1:n,survival=data$TL,dead=rep(1,n),Trt=data$Trt,data$Covariate,
	cost=data$Cost_true,tot.cost=apply(data$Cost_true,1,sum),
	QALY=data$QALY_true,tot.QALY=apply(data$QALY_true,1,sum)
	)
write.csv(datacsv.true,"True_CEdata.csv",row.names=FALSE)


####  Use potential outcomes (counterfactuals) to calculate true ICER and true INBs to evaluate methods ##########
lambda=seq(0,6,3)

#True ICER within 10-year
mean(apply(data$Cost_true1[,1:10]-data$Cost_true0[,1:10],1,sum))/1000/mean(apply(data$QALY_true1[,1:10]-data$QALY_true0[,1:10],1,sum))

#true causal average INB within 10-year, the last row is for WTP=NA (Effect only) and thus is mean extra effect instead of INB
cbind(lambda=c(lambda,NA),
	INB=c(lambda*(mean(apply(data$QALY_true1[,1:10],1,sum))-mean(apply(data$QALY_true0[,1:10],1,sum)))-
		(mean(apply(data$Cost_true1[,1:10],1,sum))-mean(apply(data$Cost_true0[,1:10],1,sum)))/1000,
		mean(apply(data$QALY_true1[,1:10],1,sum))-mean(apply(data$QALY_true0[,1:10],1,sum)))
)

