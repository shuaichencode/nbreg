########## This file provides fuctions to fit NBR, doubly robust method combining NBR and PS, and create CEAC plot ##################
# Version 1.0.3 (date: 10/16/2022)
# Author: Shuai Chen
# developed under R version 4.0.4

# Note:
# The nbreg allows factor covariates Z (but need to set as factor before using nbreg)
# Interactions between covariates and treatment must be a subset of main effects (included in Z)
# Will exclude patients with missing values 
# Although patients with missing covariates will be excluded, censoring rate is calculated using as many as possible patients.

require(geepack)	#used by CC and AL method
require(dplyr)	

###### K-M estimator (used to calculate weights of inverse probability of censoring) 		    ############	
###### assume censoring is not covariate-dependent, otherwise need to use Cox model to estimate K instead ######

KM<-function(X,delta)
# X=follow-up time
# delta=indicator of event such as death (=1 if event occurs, =0 if event does not occur)
{
  n=length(X)
  K=numeric(n)
  k<-1
  ord<-order(X)
  u=1
  while(u<=n)
  {
      s<-u+1      
      while((s<=n)&(X[ord[u]]==X[ord[s]])) s=s+1
      k<-k*(1-sum(delta[ord[u:(s-1)]]==1)/(n-u+1))
      K[ord[u:(s-1)]]<-k
      u=s
  }
  return(K)
}


####### Naive method: use complete case only for net-benefit regression, not recommended but included for comparison	#####

nbreg.CC<-function(X,delta,Cost.total,Eff.total,group,Z=NULL,interaction=NULL,lambda=NULL,Cost.only=FALSE)
# X=follow-up time
# delta=indicator of event (=1 if complete, =0 if censored)
# Cost.total=observed total costs
# Eff.total=observed total effectiveness, if not provided, assume effectiveness is survival 
# group=treatment
# lambda=WTP
# Z=covariate matrix (if not provided, will do unadjusted analysis using simple regression)
# interaction=variable names to be included in interaction terms (must be a subset of Z, otherwise will be discarded)
# Cost.only=TRUE to fit a cost-only regression
{
  lambda1=lambda
  if(is.null(lambda))lambda1=1	
  if(Cost.only)lambda1=0

  #total cost and surv among complete cases only
  if(Cost.only)DV=Cost.total else DV=lambda1*Eff.total-Cost.total
  id=1:length(X)

  if(is.null(Z))	Z=matrix(0,length(X),0) 
 
  data=data.frame(DV,X,delta,group,Z,id) %>% filter(if_all(,  ~ !is.na(.))) %>%	#keep non-missing only
		filter(`delta`==1) %>% 	#keep delta==1 only
		mutate(across(where(is.factor), droplevels)) 	#drop some levels unused in final data
  data=data	%>% select(which(sapply(data,nlevels)!=1))	#drop variables with only 1 level
  main.name=names(data)[-c(1:4,dim(data)[2])]
  Z.name=paste(main.name,collapse = "+")

  if(!is.null(interaction))
  {
	int.name= intersect(main.name,interaction)	#interaction terms must be in Z
	AZ.name=paste(paste("group",int.name,sep=":"),collapse = "+")
	if(length(int.name)>0)form=paste("DV~group+",Z.name,"+",AZ.name,sep="") else interaction=NULL
  }else int.name=character(0)

  if(is.null(interaction)){form=paste("DV~group+",Z.name,sep="")}
  if(dim(data)[2]==5){form="DV~group"}

  #use GEE to obtain robust SE
  fit.gee<-geeglm(as.formula(form),id=id,data=data)

  cat("All n =",length(X),", Used n =",dim(data)[1],"\n")
  if(is.null(lambda)) {
	if(!Cost.only) {cat("WTP (lambda) = NA (for Effect)\n");Reg.type="Effect"} else {cat("WTP (lambda) = NA (for Cost)\n");Reg.type="Cost"}
	lambda=NA
  } else {cat("WTP (lambda) =",lambda,"\n");Reg.type="NBR"}
  print(summary(fit.gee)$coef)
  cat("\n")

  est=summary(fit.gee)$coef[,1]
  se=summary(fit.gee)$coef[,2]
  if((is.na(lambda))|length(int.name)>0) CEAC=NA else CEAC=pnorm(est[2]/se[2])	#if there is interaction, does not provide CEAC here

  return(list(Method='CC',lambda=lambda,est=est,se=se,covariance=summary(fit.gee)$cov.unscaled,
	coef.table=summary(fit.gee)$coef,CEAC=CEAC,int.name=int.name,covar1st=data[1,-c(1:4,dim(data)[2]),drop=FALSE],
	Reg.type=Reg.type))			
}

####### Naive method: all data ignoring censoring status for net-benefit regression, not recommended but included for comparison	#####

nbreg.AL<-function(X,delta,Cost.total,Eff.total,group,Z=NULL,interaction=NULL,lambda=NULL,Cost.only=FALSE)
# X=follow-up time
# delta=indicator of event (=1 if complete, =0 if censored)
# Cost.total=observed total costs
# Eff.total=observed total effectiveness, if not provided, assume effectiveness is survival 
# group=treatment
# lambda=WTP
# Z=covariate matrix (if not provided, will do unadjusted analysis using simple regression)
# interaction=variable names to be included in interaction terms (must be a subset of Z, otherwise will be discarded)
# Cost.only=TRUE to fit a cost-only regression
{
  lambda1=lambda
  if(is.null(lambda))lambda1=1
  if(Cost.only)lambda1=0

  #assume observed total cost and surv are not censored
   if(Cost.only)DV=Cost.total else DV=lambda1*Eff.total-Cost.total
  id=1:length(delta)

  if(is.null(Z)) Z=matrix(0,length(delta),0) 
 
  data=data.frame(DV,X,delta,group,Z,id) %>% filter(if_all(,  ~ !is.na(.))) %>%	#keep non-missing only
		mutate(across(where(is.factor), droplevels)) 	#drop some levels unused
  data=data	%>% select(which(sapply(data,nlevels)!=1))	#drop variables with only 1 level
  main.name=names(data)[-c(1:4,dim(data)[2])]
  Z.name=paste(main.name,collapse = "+")

  if(!is.null(interaction))
  {
	int.name= intersect(main.name,interaction)	#interaction terms must be in Z
	AZ.name=paste(paste("group",int.name,sep=":"),collapse = "+")
	if(length(int.name)>0)form=paste("DV~group+",Z.name,"+",AZ.name,sep="") else interaction=NULL
  }else int.name=character(0)

  if(is.null(interaction)){form=paste("DV~group+",Z.name,sep="")}
  if(dim(data)[2]==5){form="DV~group"}

  #use GEE to obtain robust SE
  fit.gee<-geeglm(as.formula(form),id=id,data=data)

  cat("All n =",length(X),", Used n =",dim(data)[1],"\n")
  if(is.null(lambda)) {
	if(!Cost.only) {cat("WTP (lambda) = NA (for Effect)\n");Reg.type="Effect"} else {cat("WTP (lambda) = NA (for Cost)\n");Reg.type="Cost"}
	lambda=NA
  } else {cat("WTP (lambda) =",lambda,"\n");Reg.type="NBR"}
  print(summary(fit.gee)$coef)
  cat("\n")

  est=summary(fit.gee)$coef[,1]
  se=summary(fit.gee)$coef[,2]
  if((is.na(lambda))|length(int.name)>0) CEAC=NA else CEAC=pnorm(est[2]/se[2])	#if there is interaction, does not provide CEAC here

  return(list(Method='AL',lambda=lambda,est=est,se=se,covariance=summary(fit.gee)$cov.unscaled,
	coef.table=summary(fit.gee)$coef,CEAC=CEAC,int.name=int.name,covar1st=data[1,-c(1:4,dim(data)[2]),drop=FALSE],
	Reg.type=Reg.type))			
}

#######calculate Simple Weighted (SW, unpartitioned version) estimators

nbreg.SW<-function(X,delta,Cost.total,Eff.total,group,Z=NULL,interaction=NULL,Sep.K=TRUE,lambda=NULL,Cost.only=FALSE)
# X=follow-up time
# delta=indicator of event (=1 if complete, =0 if censored)
# Cost.total=observed total costs
# Eff.total=observed total effectiveness, if not provided, assume effectiveness is survival 
# group=treatment
# lambda=WTP
# Z=covariate matrix (if not provided, will do unadjusted analysis using simple regression)
# interaction=variable names to be included in interaction terms (must be a subset of Z, otherwise will be discarded)
# Sep.K=TRUE to estimate K within each treatment separately
# Cost.only=TRUE to fit a cost-only regression
{
  lambda1=lambda
  if(is.null(lambda))lambda1=1
  if(Cost.only)lambda1=0

  if(Cost.only)DV=Cost.total else DV=lambda1*Eff.total-Cost.total
  if(is.null(Z)) Z=matrix(0,length(delta),0) 
  data=data.frame(DV,X,delta,group,Z) %>% filter(if_all(,  ~ !is.na(.))) %>%	#keep non-missing only
		mutate(across(where(is.factor), droplevels)) 	#drop some levels unused
  data=data	%>% select(which(sapply(data,nlevels)!=1))	#drop variables with only 1 level
  main.name=names(data)[-c(1:4)]
  Z.name=paste(main.name,collapse = "+")
  covar1st=data[1,-c(1:4),drop=FALSE]

  #create covariates matrix, need to create dummy for factor
  if(!is.null(interaction))
  {
	int.name= intersect(main.name,interaction)	#interaction terms must be in Z
	AZ.name=paste(paste("group",int.name,sep=":"),collapse = "+")
	if(length(int.name)>0)form=paste("~group+",Z.name,"+",AZ.name,sep="") else interaction=NULL
  } else int.name=character(0)
  if(is.null(interaction)){form=paste("~group+",Z.name,sep="")}
  if(dim(data)[2]==4){form="~group"}
  covar=model.matrix(as.formula(form), data = data)

  # estimate K by K-M (assume no covariates, otherwise need to use Cox model)
  if(!Sep.K)K=KM(data$X, 1-data$delta)
  if(Sep.K){
	K=numeric(dim(data)[1])
	K[data$group==0]=KM(data$X[data$group==0], 1-data$delta[data$group==0])
	K[data$group==1]=KM(data$X[data$group==1], 1-data$delta[data$group==1])
  }
  if (sum(K[data$delta==1]==0)>0) stop("Estimate of probability of censoring becomes 0 at some time point. Time limit L may be too large.")
  if (min(K[data$delta==1])<0.1) warning("Estimate of probability of censoring <10% at some time point. To have more stable results, recommend to choose a smaller time limit L.")
  data$w0=data$delta/K				#weight
  data$w0[data$delta==0]=0

  #coef
  #to prevent singular matrix, drop variables with \sum w0*X^2=0 from covar, and use 0 as their coef est
  index.drop=which(apply(covar^2*data$w0,2,sum)==0)
  if(length(index.drop)>0)
  {
  	covar1=covar[,-index.drop]
	est2=numeric(dim(covar)[2])  
  	est2[-index.drop]=solve(t(covar1)%*%(covar1*data$w0))%*%(apply(covar1*data$w0*data$DV,2,sum))
  	est2[index.drop]=0
  } else est2=solve(t(covar)%*%(covar*data$w0))%*%(apply(covar*data$w0*data$DV,2,sum))
  
  #use formula to obtain robust SE
  A_mat=t(covar)%*%covar
  Q_mat=tmp1=matrix(0,dim(data)[1],dim(covar)[2])
  N.largerX=numeric(dim(data)[1])
  for(i in 1:dim(data)[1])
  {
	N.largerX[i]=sum(data$X>=data$X[i])
	Q_mat[i,]=apply(covar*data$w0*(data$X>data$X[i])*c(data$DV-covar%*%est2),2,sum)/N.largerX[i]
  }
  for(i in 1:dim(data)[1])
	tmp1[i,]=apply(Q_mat*(1-data$delta)*(data$X<=data$X[i])/N.largerX,2,sum)	
  tmp2=covar*data$w0*c(data$DV-covar%*%est2)+Q_mat*(1-data$delta)-tmp1
  B_mat=t(tmp2)%*%tmp2
  cov_mat=solve(A_mat)%*%B_mat%*%solve(A_mat)
  se2=sqrt(diag(cov_mat))

  coef.table=data.frame(Estimate=est2,Std.err=se2,Wald=(est2/se2)^2,p=2-2*pnorm(abs(est2/se2)))
  row.names(coef.table)=colnames(covar)
  cat("All n =",length(X),", Used n =",dim(data)[1],"\n")
  if(is.null(lambda)) {
	if(!Cost.only) {cat("WTP (lambda) = NA (for Effect)\n");Reg.type="Effect"} else {cat("WTP (lambda) = NA (for Cost)\n");Reg.type="Cost"}
	lambda=NA
  } else {cat("WTP (lambda) =",lambda,"\n");Reg.type="NBR"}
  print(coef.table)
  cat("\n")

  if((is.na(lambda))|length(int.name)>0) CEAC=NA else CEAC=pnorm(est2[2]/se2[2])	#if there is interaction, does not provide CEAC here

  return(list(Method='SW',lambda=lambda,est=est2,se=se2,covariance=cov_mat,coef.table=coef.table,CEAC=CEAC,
	int.name=int.name,covar1st=covar1st,Reg.type=Reg.type))
}



#######calculate partitioned (PT) estimators############

nbreg.PT<-function(X,delta,Cost.grp,Eff.grp,Patition.times,group,Z=NULL,interaction=NULL,Sep.K=TRUE,lambda=NULL,Cost.only=FALSE)
# X=follow-up time
# delta=indicator of event (=1 if complete, =0 if censored)
# Cost.grp=observed grouped 30-day costs, Cost.grp[i,j]=observed cost of ith people accumulated in jth interval
# Eff.grp=observed grouped 30-day effectiveness, if not provided, assume effectiveness is survival 
# Patition.times=time points used to partition the time
# group=treatment
# lambda=WTP
# Z=covariate matrix (if not provided, will do unadjusted analysis using simple regression)
# interaction=variable names to be included in interaction terms (must be a subset of Z, otherwise will be discarded)
# Sep.K=TRUE to estimate K within each treatment separately
# Cost.only=TRUE to fit a cost-only regression
{
  lambda1=lambda
  if(is.null(lambda))lambda1=1
  if(Cost.only)lambda1=0

  if(is.null(Z)) Z=matrix(0,length(delta),0) 

  data=data.frame(Cost.grp,Eff.grp,X,delta,group,Z) %>% filter(if_all(,  ~ !is.na(.))) %>%	#keep non-missing only
		mutate(across(where(is.factor), droplevels)) 	#drop some levels unused
  data=data	%>% select(which(sapply(data,nlevels)!=1))	#drop variables with only 1 level
  main.name=names(data)[-c(1:(dim(Cost.grp)[2]+dim(Eff.grp)[2]+3))]
  Z.name=paste(main.name,collapse = "+")
  covar1st=data[1,-c(1:(dim(Cost.grp)[2]+dim(Eff.grp)[2]+3)),drop=FALSE]

  #create covariates matrix, need to create dummy for factor
  if(!is.null(interaction))
  {
	int.name= intersect(main.name,interaction)	#interaction terms must be in Z
	AZ.name=paste(paste("group",int.name,sep=":"),collapse = "+")
	if(length(int.name)>0)form=paste("~group+",Z.name,"+",AZ.name,sep="") else interaction=NULL
  } else int.name=character(0)
  if(is.null(interaction)){form=paste("~group+",Z.name,sep="")}
  if(dim(data)[2]==(3+(dim(Cost.grp)[2]+dim(Eff.grp)[2]))){form="~group"}
  covar=model.matrix(as.formula(form), data = data)

  x=w0=d=totalQ=DV=matrix(0,dim(data)[1],length(Patition.times)-1)
  K=numeric(dim(data)[1])
  coef.split=matrix(0,length(Patition.times)-1,dim(covar)[2])
  Ksi=array(0,c(length(Patition.times)-1,dim(data)[1],dim(covar)[2]))
  B.split=array(0,c(length(Patition.times)-1,length(Patition.times)-1,dim(covar)[2],dim(covar)[2]))

  Cost.grp1=data[,1:dim(Cost.grp)[2]]
  Eff.grp1=data[,(dim(Cost.grp)[2]+1):(dim(Cost.grp)[2]+dim(Eff.grp)[2])]

## loop for each month
for(j in 2:length(Patition.times))
{
   trunc.time=Patition.times[j] 
 
   # truncate follow-up time for each interval
   # set death indicator to be 1 if truncated
   for(l in 1:dim(data)[1]) 
   {
      x[l,j-1]<-min(data$X[l],trunc.time)          #new follow-up x[l,j-1] of lth patient for j-1th interval
      if(data$X[l]>=trunc.time)  d[l,j-1]<-1  else d[l,j-1]=data$delta[l]         #d[l,j]=new death indicator of lth patient for j-1th interval
   }

# re-estimate K by K-M using new data(assume no covariates, otherwise need to use Cox model)
# (theoretically not needed but better coverage in simulation) 
  if(!Sep.K)K=KM(x[,j-1],1-d[,j-1])
  if(Sep.K){
	K[data$group==0]=KM(x[data$group==0,j-1], 1-d[data$group==0,j-1])
	K[data$group==1]=KM(x[data$group==1,j-1], 1-d[data$group==1,j-1])
  }
  if (sum(K[d[,j-1]==1]==0)>0) stop("Estimate of probability of censoring becomes 0 at some time point. Time limit L may be too large.")
  if (min(K[d[,j-1]==1])<0.1) warning("Estimate of probability of censoring <10% at some time point. To have more stable results, recommend to choose a smaller time limit L.")

  #weight need to be adjusted for x[l,j-1]<-min(X[l],trunc.time) due to truncated time 
  w0[,j-1]=d[,j-1]/K	
  w0[d[,j-1]==0,j-1]=0

  #coef
  if(Cost.only)DV[,j-1]=Cost.grp1[,j-1] else DV[,j-1]=lambda1*Eff.grp1[,j-1]-Cost.grp1[,j-1]	#new dependent variable within this interval

  #to prevent singular matrix, drop variables with \sum w0*X^2=0 from covar, and use 0 as their coef est
  index.drop=which(apply(covar^2*w0[,j-1],2,sum)==0)
  if(length(index.drop)>0)
  {
  	covar1=covar[,-index.drop]  
  	coef.split[j-1,-index.drop]=solve(t(covar1)%*%(covar1*w0[,j-1]))%*%(apply(covar1*w0[,j-1]*DV[,j-1],2,sum))
  	coef.split[j-1,index.drop]=0
  } else coef.split[j-1,]=solve(t(covar)%*%(covar*w0[,j-1]))%*%(apply(covar*w0[,j-1]*DV[,j-1],2,sum))

}

#sum up all coef across intervals
est=apply(coef.split,2,sum)

#use formula to obtain robust SE
## loop for each split time interval
for(j in 2:length(Patition.times))
{
  Q_mat=tmp1=matrix(0,dim(data)[1],dim(covar)[2])
  N.largerX=numeric(dim(data)[1])
  for(i in 1:dim(data)[1])
  {
	N.largerX[i]=sum(x[,j-1]>=x[i,j-1])
	Q_mat[i,]=apply(covar*w0[,j-1]*(x[,j-1]>x[i,j-1])*c(DV[,j-1]-covar%*%coef.split[j-1,]),2,sum)/N.largerX[i]
  }
  for(i in 1:dim(data)[1])
	tmp1[i,]=apply(Q_mat*(1-d[,j-1])*(x[,j-1]<=x[i,j-1])/N.largerX,2,sum)	
  Ksi[j-1,,]=covar*w0[,j-1]*c(DV[,j-1]-covar%*%coef.split[j-1,])+Q_mat*(1-d[,j-1])-tmp1
}

for(j in 2:length(Patition.times))
for(k in 2:length(Patition.times))
{
  B.split[j-1,k-1,,]=t(Ksi[j-1,,])%*%Ksi[k-1,,]
}

A_mat=t(covar)%*%covar
B_mat=apply(B.split,c(3,4),sum)
cov_mat=solve(A_mat)%*%B_mat%*%solve(A_mat)

  coef.table=data.frame(Estimate=est,Std.err=sqrt(diag(cov_mat)),
	Wald=(est/sqrt(diag(cov_mat)))^2,p=2-2*pnorm(abs(est/sqrt(diag(cov_mat)))))
  row.names(coef.table)=colnames(covar)
  cat("All n =",length(X),", Used n =",dim(data)[1],"\n")
  if(is.null(lambda)) {
	if(!Cost.only) {cat("WTP (lambda) = NA (for Effect)\n");Reg.type="Effect"} else {cat("WTP (lambda) = NA (for Cost)\n");Reg.type="Cost"}
	lambda=NA
  } else {cat("WTP (lambda) =",lambda,"\n");Reg.type="NBR"}
  print(coef.table)
  cat("\n")

  if((is.na(lambda))|length(int.name)>0) CEAC=NA else CEAC=pnorm(est[2]/sqrt(diag(cov_mat))[2])	#if there is interaction, does not provide CEAC here

return(list(Method='PT',lambda=lambda,est=est,se=sqrt(diag(cov_mat)),covariance=cov_mat,coef.table=coef.table,
	CEAC=CEAC,int.name=int.name,covar1st=covar1st,Reg.type=Reg.type))#,A_mat=A_mat,Ksi=Ksi))
}



#######calculate DR Simple Weighted (SW, unpartitioned version) estimators

nbreg.SW.DR<-function(X,delta,Cost.total,Eff.total,group,Z=NULL,PS.Z=NULL,interaction=NULL,Sep.K=TRUE,PS.trim=0.01,lambda=NULL,DR.Reg.SE=TRUE,Cost.only=FALSE)
# X=follow-up time
# delta=indicator of event (=1 if complete, =0 if censored)
# Cost.total=observed total costs
# Eff.total=observed total effectiveness, if not provided, assume effectiveness is survival 
# group=treatment
# lambda=WTP
# Z=covariate matrix for NBR model (if not provided, will do unadjusted analysis using simple regression)
# PS.Z=covariate matrix for PS model (if not provided, will do unadjusted analysis)
# interaction=variable names to be included in interaction with treatment in NBR model (must be a subset of Z, otherwise will be discarded)
# Sep.K=TRUE to estimate K within each treatment separately
# PS.trim= threshold to trim extreme PS values outside (PS.trim, 1-PS.trim)
# DR.Reg.SE=TRUE to request SE (and p-value) reported for the regression part of doubly robust method, otherwise only coef is reported
# Cost.only=TRUE to fit a cost-only regression
{
  lambda1=lambda
  if(is.null(lambda))lambda1=1
  if(Cost.only)lambda1=0

  ### create design matrix for NBR and PS, need to create dummy for factor
  if(Cost.only)DV=Cost.total else DV=lambda1*Eff.total-Cost.total
  if(is.null(Z)) Z=matrix(0,length(delta),0) 
  if(is.null(PS.Z)) PS.Z=matrix(0,length(delta),0) 

  data=data.frame(DV,X,delta,group,Z,PS.Z) %>% filter(if_all(,  ~ !is.na(.))) %>%	#keep non-missing only
		mutate(across(where(is.factor), droplevels)) 	#drop some levels unused
  data=data	%>% select(which(sapply(data,nlevels)!=1))	#drop variables with only 1 level
  if(dim(Z)[2]>0)Z.name=paste(names(data)[5:(dim(Z)[2]+4)],collapse = "+")else Z.name=character(0)
  if(!is.null(interaction))
  {
	int.name= intersect(names(data)[5:(dim(Z)[2]+4)],interaction)	#interaction terms must be in Z
	AZ.name=paste(paste("group",int.name,sep=":"),collapse = "+")
	if(length(int.name)>0)form=paste("DV~group+",Z.name,"+",AZ.name,sep="") else interaction=NULL
  } else int.name=character(0)
  if(is.null(interaction)){form=paste("DV~group+",Z.name,sep="")}
  if(dim(Z)[2]==0){form="DV~group"}
  covar=model.matrix(as.formula(form), data = data)

  # estimate K by K-M (assume no covariates, otherwise need to use Cox model)
  if(!Sep.K)data$K=KM(data$X, 1-data$delta)
  if(Sep.K){
	data$K=0
	data$K[data$group==0]=KM(data$X[data$group==0], 1-data$delta[data$group==0])
	data$K[data$group==1]=KM(data$X[data$group==1], 1-data$delta[data$group==1])
  }
  if (sum(data$K[data$delta==1]==0)>0) stop("Estimate of probability of censoring becomes 0 at some time point. Time limit L may be too large.")
  if (min(data$K[data$delta==1])<0.1) warning("Estimate of probability of censoring <10% at some time point. To have more stable results, recommend to choose a smaller time limit L.")
  data$w0=data$delta/data$K				#weight
  data$w0[data$delta==0]=0

  #create new data used in prediction
  datanew1=datanew0=data
  datanew1$group=1
  datanew0$group=0
  new1=model.matrix(as.formula(form), data = datanew1)
  new0=model.matrix(as.formula(form), data = datanew0) 

  #estimate PS by logistic reg
  #create formula for PS
  if(dim(PS.Z)[2]>0){
	data.PS=data[,(dim(Z)[2]+5):(dim(Z)[2]+dim(PS.Z)[2]+4)]
	names(data.PS)=colnames(PS.Z)
	PS.Z.name=paste(names(data.PS),collapse = "+")
  }else PS.Z.name=character(0)
  form=paste("group~",PS.Z.name,sep="")
  if(dim(PS.Z)[2]==0){form="group~1"}
  PS.covar=model.matrix(as.formula(form), data = data.PS)

  # Estimate propensity score model
  PSmodel=glm(as.formula(form), data = data.PS,family="binomial")
  PS=predict(PSmodel, type='response')
  PScoef=summary(PSmodel)$coef

  #trim PS to prevent extreme value
  PS[PS<PS.trim]=PS.trim
  PS[PS>1-PS.trim]=1-PS.trim

  #NBR coef
  #to prevent singular matrix, drop variables with \sum w0*X^2=0 from covar, and use 0 as their coef est
  index.drop=which(apply(covar^2*data$w0,2,sum)==0)
  if(length(index.drop)>0)
  {
  	covar1=covar[,-index.drop]
	est.coef=numeric(dim(covar)[2])  
  	est.coef[-index.drop]=solve(t(covar1)%*%(covar1*data$w0))%*%(apply(covar1*data$w0*data$DV,2,sum))
  	est.coef[index.drop]=0
  } else est.coef=solve(t(covar)%*%(covar*data$w0))%*%(apply(covar*data$w0*data$DV,2,sum))
  colnames(est.coef)="Estimate"

  m0=new0%*%est.coef
  m1=new1%*%est.coef

  #robust SE for regression part
  if(DR.Reg.SE){
	A_mat=t(covar)%*%covar
	Q_mat=tmp1=matrix(0,dim(data)[1],dim(covar)[2])
	N.largerX=numeric(dim(data)[1])
	for(i in 1:dim(data)[1])
	{
		N.largerX[i]=sum(data$X>=data$X[i])
		Q_mat[i,]=apply(covar*data$w0*(data$X>data$X[i])*c(data$DV-covar%*%est.coef),2,sum)/N.largerX[i]
  	}
	for(i in 1:dim(data)[1])
		tmp1[i,]=apply(Q_mat*(1-data$delta)*(data$X<=data$X[i])/N.largerX,2,sum)	
 	 tmp2=covar*data$w0*c(data$DV-covar%*%est.coef)+Q_mat*(1-data$delta)-tmp1
  	B_mat=t(tmp2)%*%tmp2
  	cov_mat=solve(A_mat)%*%B_mat%*%solve(A_mat)
  	se2=sqrt(diag(cov_mat))

  	reg.coef.table=data.frame(Estimate=est.coef,Std.err=se2,Wald=(c(est.coef)/se2)^2,p=2-2*pnorm(abs(c(est.coef)/se2)))
  	row.names(reg.coef.table)=colnames(covar)
  }else reg.coef.table=est.coef

  #DR est
  mu1=sum(data$w0*(data$group*data$DV/PS-(data$group-PS)/PS*m1))/sum(data$w0)  
  mu0=sum(data$w0*((1-data$group)*data$DV/(1-PS)+(data$group-PS)/(1-PS)*m0))/sum(data$w0)
  est=mu1-mu0	

  #use our formula to obtain SE (seems overestimate SE slightly)
  W=data$group*data$DV/PS-(data$group-PS)/PS*m1-(1-data$group)*data$DV/(1-PS)-(data$group-PS)/(1-PS)*m0-est
  W0=apply(c(data$w0*data$group/PS*(data$DV-m1)*(1-PS))*PS.covar,2,sum)/sum(data$w0*data$group/PS)+
	apply(c(data$w0*(1-data$group)/(1-PS)*(data$DV-m0)*PS)*PS.covar,2,sum)/sum(data$w0*(1-data$group)/(1-PS))
  S=c(data$group-PS)*PS.covar
  ESS=(t(PS.covar)%*%diag(PS*(1-PS))%*%PS.covar)/dim(data)[1]

  G0=G1=G02=G12=N.largerX=numeric(dim(data)[1])
  Y_bar=matrix(NA,dim(data)[1],dim(data)[1])
  for(i in 1:dim(data)[1])
  {
	Y_bar[,i]=(data$X>=data$X[i])
	G1[i]=mean(data$w0*data$group*(data$X>=data$X[i])*(W-sum(data$group*data$w0*(data$X>=data$X[i])*W)/sum(data$group*data$w0*(data$X>=data$X[i]))))
	G0[i]=mean(data$w0*(1-data$group)*(data$X>=data$X[i])*(W-sum((1-data$group)*data$w0*(data$X>=data$X[i])*W)/sum((1-data$group)*data$w0*(data$X>=data$X[i]))))
	G12[i]=mean(data$w0*data$group*(data$X>=data$X[i])*(W-sum(data$group*data$w0*(data$X>=data$X[i])*W)/sum(data$group*data$w0*(data$X>=data$X[i])))^2)
	G02[i]=mean(data$w0*(1-data$group)*(data$X>=data$X[i])*(W-sum((1-data$group)*data$w0*(data$X>=data$X[i])*W)/sum((1-data$group)*data$w0*(data$X>=data$X[i])))^2)
  }
  phi=data$w0*(c(W)-c(W0%*%solve(ESS)%*%t(S)))+data$group*(1-data$delta)*G1/data$K/apply(data$group*Y_bar,2,sum)+
	(1-data$group)*(1-data$delta)*G0/data$K/apply((1-data$group)*Y_bar,2,sum)
  DRvar=sum(phi^2)/dim(data)[1]/dim(data)[1]

  coef.table=data.frame(Estimate=est,Std.err=sqrt(DRvar),Wald=(est/sqrt(DRvar))^2,p=2-2*pnorm(abs(est/sqrt(DRvar))))
  row.names(coef.table)="group"

  cat("All n =",length(X),", Used n =",dim(data)[1],"\n")
  if(is.null(lambda)) {
	if(!Cost.only) {cat("WTP (lambda) = NA (for Effect)\n");Reg.type="Effect"} else {cat("WTP (lambda) = NA (for Cost)\n");Reg.type="Cost"}
	lambda=NA
  } else {cat("WTP (lambda) =",lambda,"\n");Reg.type="NBR"}
  print(coef.table)
  cat("\n")

  if(is.na(lambda))CEAC=NA else CEAC=pnorm(est/sqrt(DRvar))

  return(list(Method='DR.SW',lambda=lambda,est=est,se=sqrt(DRvar),coef.table=coef.table,CEAC=CEAC,
	Regmodel=reg.coef.table,PSmodel=PScoef,PS=PS,group=data$group,Reg.type=Reg.type))
}


#######calculate DR partitioned estimators############

nbreg.PT.DR<-function(X,delta,Cost.grp,Eff.grp,Patition.times,group,Z=NULL,PS.Z=NULL,interaction=NULL,Sep.K=TRUE,PS.trim=0.01,lambda=NULL,DR.Reg.SE=TRUE,Cost.only=FALSE)
# X=follow-up time
# delta=indicator of event (=1 if complete, =0 if censored)
# Cost.grp=observed grouped 30-day costs, Cost.grp[i,j]=observed cost of ith people accumulated in jth interval
# Eff.grp=observed grouped 30-day effectiveness, if not provided, assume effectiveness is survival 
# Patition.times=time points used to partition the time
# group=treatment
# lambda=WTP
# Z=covariate matrix for NBR model (if not provided, will do unadjusted analysis using simple regression)
# PS.Z=covariate matrix for PS model (if not provided and PS is not provided, will do unadjusted analysis)
# interaction=variable names to be included in interaction terms (must be a subset of Z, otherwise will be discarded)
# Sep.K=TRUE to estimate K within each treatment separately
# PS.trim= threshold to trim extreme PS values outside (PS.trim, 1-PS.trim)
# DR.Reg.SE=TRUE to request SE (and p-value) reported for the regression part of doubly robust method, otherwise only coef is reported
# Cost.only=TRUE to fit a cost-only regression
{
  lambda1=lambda
  if(is.null(lambda))lambda1=1
  if(Cost.only)lambda1=0

  #create covariates matrix, need to create dummy for factor
  if(is.null(Z)) Z=matrix(0,length(delta),0) 
  if(is.null(PS.Z)) PS.Z=matrix(0,length(delta),0) 

  data=data.frame(Cost.grp,Eff.grp,X,delta,group,Z,PS.Z) %>% filter(if_all(,  ~ !is.na(.))) %>%	#keep non-missing only
		mutate(across(where(is.factor), droplevels)) 	#drop some levels unused
  data=data	%>% select(which(sapply(data,nlevels)!=1))	#drop variables with only 1 level
  if(dim(Z)[2]>0)Z.name=paste(names(data)[(dim(Cost.grp)[2]+dim(Eff.grp)[2]+4):
	(dim(Z)[2]+dim(Cost.grp)[2]+dim(Eff.grp)[2]+3)],collapse = "+")else Z.name=character(0)

  if(!is.null(interaction))
  {
	int.name= intersect(names(data)[(dim(Cost.grp)[2]+dim(Eff.grp)[2]+4):(dim(Z)[2]+dim(Cost.grp)[2]+dim(Eff.grp)[2]+3)],interaction)	#interaction terms must be in Z
	AZ.name=paste(paste("group",int.name,sep=":"),collapse = "+")
	if(length(int.name)>0)form=paste("~group+",Z.name,"+",AZ.name,sep="") else interaction=NULL
  } else int.name=character(0)
  if(is.null(interaction)){form=paste("~group+",Z.name,sep="")}
  if(dim(Z)[2]==0){form="~group"}
  covar=model.matrix(as.formula(form), data = data)

  #create new data used in prediction
  datanew1=datanew0=data
  datanew1$group=1
  datanew0$group=0
  new1=model.matrix(as.formula(form), data = datanew1)
  new0=model.matrix(as.formula(form), data = datanew0) 

  #estimate PS by logistic reg
  #create formula for PS
  if(dim(PS.Z)[2]>0){
	data.PS=data[,(dim(Z)[2]+dim(Cost.grp)[2]+dim(Eff.grp)[2]+4):(dim(Z)[2]+dim(PS.Z)[2]+dim(Cost.grp)[2]+dim(Eff.grp)[2]+3)]
	names(data.PS)=colnames(PS.Z)
	PS.Z.name=paste(names(data.PS),collapse = "+")
  }else PS.Z.name=character(0)
  form=paste("group~",PS.Z.name,sep="")
  if(dim(PS.Z)[2]==0){form="group~1"}
  PS.covar=model.matrix(as.formula(form), data = data.PS)

  # Estimate propensity score model
  PSmodel=glm(as.formula(form), data = data.PS,family="binomial")
  PS=predict(PSmodel, type='response')
  PScoef=summary(PSmodel)$coef

  #trim PS to prevent extreme value
  PS[PS<PS.trim]=PS.trim
  PS[PS>1-PS.trim]=1-PS.trim

coef.split=coef.split2=matrix(0,length(Patition.times)-1,dim(covar)[2])
Ksi=array(0,c(length(Patition.times)-1,dim(data)[1],dim(covar)[2]))
B.split=array(0,c(length(Patition.times)-1,length(Patition.times)-1,dim(covar)[2],dim(covar)[2]))
x=w0=w2=d=K.new=m1=m0=totalQ=DV=matrix(0,dim(data)[1],length(Patition.times)-1)
mu1=mu0=est.IPW=numeric(length(Patition.times)-1)
K0=K=numeric(dim(data)[1])

  Cost.grp1=data[,1:dim(Cost.grp)[2]]
  Eff.grp1=data[,(dim(Cost.grp)[2]+1):(dim(Cost.grp)[2]+dim(Eff.grp)[2])]

## loop for each split time interval 
for(j in 2:length(Patition.times))
{
   trunc.time=Patition.times[j]
 
   # truncate follow-up time for each interval
   # set death indicator to be 1 if truncated
   for(l in 1:dim(data)[1]) 
   {
      x[l,j-1]<-min(data$X[l],trunc.time)          #new follow-up x[l,j-1] of lth patient for j-1th interval
      if(data$X[l]>=trunc.time)  d[l,j-1]<-1  else d[l,j-1]=data$delta[l]         #d[l,j]=new death indicator of lth patient for j-1th interval
   }

# re-estimate K by K-M using new data(assume no covariates, otherwise need to use Cox model)
# (theoretically not needed but better coverage in simulation) 
  if(!Sep.K)K=KM(x[,j-1],1-d[,j-1])
  if(Sep.K){
	K[data$group==0]=KM(x[data$group==0,j-1], 1-d[data$group==0,j-1])
	K[data$group==1]=KM(x[data$group==1,j-1], 1-d[data$group==1,j-1])
  }
  if (sum(K[d[,j-1]==1]==0)>0) stop("Estimate of probability of censoring becomes 0 at some time point. Time limit L may be too large.")
  if (min(K[d[,j-1]==1])<0.1) warning("Estimate of probability of censoring <10% at some time point. To have more stable results, recommend to choose a smaller time limit L.")

  #weight need to be adjusted for x[l,j-1]<-min(X[l],trunc.time) due to truncated time 
  w0[,j-1]=d[,j-1]/K	#weight re-estimated
  w0[d[,j-1]==0,j-1]=0
  K.new[,j-1]=K

  if(Cost.only)DV[,j-1]=Cost.grp1[,j-1] else DV[,j-1]=lambda1*Eff.grp1[,j-1]-Cost.grp1[,j-1]	#new dependent variable within this interval

  #coef, to prevent singular matrix, drop variables with \sum w0*X^2=0 from covar, and use 0 as their coef est
  index.drop=which(apply(covar^2*w0[,j-1],2,sum)==0)
  if(length(index.drop)>0)
  {
  	covar1=covar[,-index.drop]  
  	coef.split[j-1,-index.drop]=solve(t(covar1)%*%(covar1*w0[,j-1]))%*%(apply(covar1*w0[,j-1]*DV[,j-1],2,sum))
  	coef.split[j-1,index.drop]=0
  } else coef.split[j-1,]=solve(t(covar)%*%(covar*w0[,j-1]))%*%(apply(covar*w0[,j-1]*DV[,j-1],2,sum))

  m0[,j-1]=new0%*%coef.split[j-1,]
  m1[,j-1]=new1%*%coef.split[j-1,]

  mu1[j-1]=sum(w0[,j-1]*(data$group*DV[,j-1]/PS-(data$group-PS)/PS*m1[,j-1]))/sum(w0[,j-1])  
  mu0[j-1]=sum(w0[,j-1]*((1-data$group)*DV[,j-1]/(1-PS)+(data$group-PS)/(1-PS)*m0[,j-1]))/sum(w0[,j-1])
}

  mu1.sum=sum(mu1)  
  mu0.sum=sum(mu0)
  est=mu1.sum-mu0.sum	#est by DR

  # coef of regression part
  est.coef=data.frame(Estimate=apply(coef.split,2,sum))
  rownames(est.coef)=colnames(covar)

  #robust SE for regression part
  if(DR.Reg.SE){
	## loop for each split time interval
	for(j in 2:length(Patition.times))
	{
		Q_mat=tmp1=matrix(0,dim(data)[1],dim(covar)[2])
		N.largerX=numeric(dim(data)[1])
		for(i in 1:dim(data)[1])
		{
			N.largerX[i]=sum(x[,j-1]>=x[i,j-1])
			Q_mat[i,]=apply(covar*w0[,j-1]*(x[,j-1]>x[i,j-1])*c(DV[,j-1]-covar%*%coef.split[j-1,]),2,sum)/N.largerX[i]
  		}
  		for(i in 1:dim(data)[1])
			tmp1[i,]=apply(Q_mat*(1-d[,j-1])*(x[,j-1]<=x[i,j-1])/N.largerX,2,sum)	
  		Ksi[j-1,,]=covar*w0[,j-1]*c(DV[,j-1]-covar%*%coef.split[j-1,])+Q_mat*(1-d[,j-1])-tmp1
	}

	for(j in 2:length(Patition.times))
	for(k in 2:length(Patition.times))
	{
  		B.split[j-1,k-1,,]=t(Ksi[j-1,,])%*%Ksi[k-1,,]
	}

	A_mat=t(covar)%*%covar
	B_mat=apply(B.split,c(3,4),sum)
	cov_mat=solve(A_mat)%*%B_mat%*%solve(A_mat)

  	reg.coef.table=data.frame(Estimate=est.coef,Std.err=sqrt(diag(cov_mat)),
		Wald=(est.coef/sqrt(diag(cov_mat)))^2,p=2-2*pnorm(abs(apply(coef.split,2,sum)/sqrt(diag(cov_mat)))))
  	row.names(reg.coef.table)=colnames(covar)
  }else reg.coef.table=est.coef

#use formula to obtain SE
## loop for each split time interval
S=c(data$group-PS)*PS.covar
ESS=(t(PS.covar)%*%diag(PS*(1-PS))%*%PS.covar)/dim(data)[1]

G0=G1=W=phi=matrix(0,dim(data)[1],length(Patition.times)-1)
W0=matrix(0,dim(PS.covar)[2],length(Patition.times)-1)
Y_bar=array(NA,c(dim(data)[1],dim(data)[1],length(Patition.times)-1))
for(j in 2:length(Patition.times))
{
  W[,j-1]=data$group*DV[,j-1]/PS-(data$group-PS)/PS*m1[,j-1]-(1-data$group)*DV[,j-1]/(1-PS)-(data$group-PS)/(1-PS)*m0[,j-1]-(mu1[j-1]-mu0[j-1])
  W0[,j-1]=apply(c(w0[,j-1]*data$group/PS*(DV[,j-1]-m1[,j-1])*(1-PS))*PS.covar,2,sum)/sum(w0[,j-1]*data$group/PS)+
	apply(c(w0[,j-1]*(1-data$group)/(1-PS)*(DV[,j-1]-m0[,j-1])*PS)*PS.covar,2,sum)/sum(w0[,j-1]*(1-data$group)/(1-PS))
  for(i in 1:dim(data)[1])
  {
	Y_bar[,i,j-1]=(x[,j-1]>=x[i,j-1])
	G1[i,j-1]=mean(w0[,j-1]*data$group*(x[,j-1]>=x[i,j-1])*(W[,j-1]-sum(data$group*w0[,j-1]*(x[,j-1]>=x[i,j-1])*W[,j-1])/sum(data$group*w0[,j-1]*(x[,j-1]>=x[i,j-1]))))
	G0[i,j-1]=mean(w0[,j-1]*(1-data$group)*(x[,j-1]>=x[i,j-1])*(W[,j-1]-sum((1-data$group)*w0[,j-1]*(x[,j-1]>=x[i,j-1])*W[,j-1])/sum((1-data$group)*w0[,j-1]*(x[,j-1]>=x[i,j-1]))))
  }
  phi[,j-1]=w0[,j-1]*(c(W[,j-1])-c(W0[,j-1]%*%solve(ESS)%*%t(S)))+
	data$group*(1-d[,j-1])*G1[,j-1]/K.new[,j-1]/apply(data$group*Y_bar[,,j-1],2,sum)+(1-data$group)*(1-d[,j-1])*G0[,j-1]/K.new[,j-1]/apply((1-data$group)*Y_bar[,,j-1],2,sum)

}

DRvar=sum(t(phi)%*%phi)/dim(data)[1]/dim(data)[1]

  coef.table=data.frame(Estimate=est,Std.err=sqrt(DRvar),Wald=(est/sqrt(DRvar))^2,p=2-2*pnorm(abs(est/sqrt(DRvar))))
  row.names(coef.table)="group"

  cat("All n =",length(X),", Used n =",dim(data)[1],"\n")
  if(is.null(lambda)) {
	if(!Cost.only) {cat("WTP (lambda) = NA (for Effect)\n");Reg.type="Effect"} else {cat("WTP (lambda) = NA (for Cost)\n");Reg.type="Cost"}
	lambda=NA
  } else {cat("WTP (lambda) =",lambda,"\n");Reg.type="NBR"}
  print(coef.table)
  cat("\n")

if(is.na(lambda))CEAC=NA else CEAC=pnorm(est/sqrt(DRvar))

return(list(Method='DR.PT',lambda=lambda,est=est,se=sqrt(DRvar),coef.table=coef.table,CEAC=CEAC,
	Regmodel=reg.coef.table,PSmodel=PScoef,PS=PS,group=data$group,Reg.type=Reg.type))
}



###### function to perform net-benefit regression (NBR) which will call above estimation functions		    ############

nbreg<-function(Followup,delta,group,Cost=NULL,Eff=NULL,Patition.times=NULL,Z=NULL,PS.Z=NULL,interaction=NULL,
	Method=c('SW','PT','CC','AL'),Sep.K=TRUE,PS.trim=0.05,Doubly.Robust=FALSE,DR.Reg.SE=TRUE,Eff.only=FALSE,Cost.only=FALSE,lambda=NULL,L)
# Followup=follow-up time of length n
# delta=indicator of event (=1 if complete, =0 if censored) of length n
# group=treatment indicator (considered treatment=1, comparison=0) of length n
# Cost=n x m matrix with observed grouped costs, Cost.grp[i,j]=observed cost of ith people accumulated in jth interval; can be a vector or one column matrix for only total costs available
# Eff=n x m matrix with observed grouped effectiveness; if not provided, assume effectiveness is survival; can be a vector or one column matrix for only total effectiveness available 
# Patition.times= m-length vector for time points (in days) used to partition the time (ie, end day of each time interval), needed to be monontonically increasing (eg, 1st interval is Day 1 to Patition.times[1], 2nd is Day Patition.times[1]+1 to Patition.times[2]); 
#	Patition.times is needed if using PT method without Eff provided, Patition.times is also needed to truncate grouped costs and effectiveness if outside L. 
# lambda=WTP, can be a vector
# Z=n x p covariate matrix for NBR model (if not provided, will do unadjusted analysis using simple regression)
# PS.Z=covariate matrix for PS model using logistic regression, needed only for doubly robust method (if not provided, will do unadjusted analysis)
# interaction=variable names to be included in interaction terms (must be a subset of Z, otherwise will be discarded)
# Method= one of 'SW' (simple weighted), 'PT' (partitioned), 'CC' (naive complete case), 'AL' (naive all data); SW and PT also can be used for doubly robust method 
# Sep.K=TRUE to estimate K within each treatment separately
# PS.trim= threshold to trim extreme PS values outside (PS.trim, 1-PS.trim), used only for doubly robust method
# Doubly.Robust=TRUE to estimate INB using doubly robust method 
# DR.Reg.SE=TRUE to request SE (and p-value) reported for the regression part of doubly robust method, otherwise only coef is reported
# L=time limit horizon, will be used to truncate costs and effectiveness if they are outside this time limit (assume cost & Eff are evenly spreaded within each time interval)
# Eff.only=TRUE to fit a regression for effectiveness only (does not belong to NBR)
# Cost.only=TRUE to fit a regression for cost only (equivalent to setting WTP lambda=0 and then switching sign of coef est)
{
  ##### check data for possible errors in input data

  if(!is.numeric(L)) stop("L must be numeric.\n")
  if((!Eff.only)&(!Cost.only)){
	if(is.null(lambda)) stop("lambda is required if not fitting effectiveness or cost-only regression.\n") 
  }
  if(!is.null(lambda))if(!is.numeric(lambda)) stop("lambda must be numeric.\n")
  if(!is.numeric(Followup)) stop("Followup must be a numeric vector.\n")
#  if(!is.numeric(Cost)) stop("Cost must be numeric.\n")
#  if(!is.null(Eff))if(!is.numeric(Eff)) stop("Eff must be numeric.\n")

  Followup=unlist(Followup)
  delta=unlist(as.numeric(as.character(delta)))
  group=unlist(as.numeric(as.character(group)))
  lambda=sort(unlist(lambda))
  
  n=length(Followup)
  if(length(delta)!=n) stop("Length of delta is different to length of Followup.\n")
  if(length(group)!=n) stop("Length of group is different to length of Followup.\n") 
  if(PS.trim<0) {warning("PS.trim is negative. Propensity scores are not trimmed.\n");PS.trim=0}
  if(PS.trim>=0.5) stop("PS.trim is too large. Recommend to be between 0 and 0.1.\n")

  if(!identical(levels(as.factor(group)),c("0","1"))) stop("Values in group must be 0 or 1 and cannot be in one group only.\n")

  if(length(unique(delta)) == 1){
		if(!(levels(as.factor(delta))%in%c("0","1")))stop("Values in delta must be 0 or 1.\n")
  }else if(!identical(levels(as.factor(delta)),c("0","1"))) stop("Values in delta must be 0 or 1.\n")
  
  if(min(Followup,na.rm=TRUE)<=0) stop("Values in Followup cannot be zero or negative.\n")
  if(!is.null(lambda))if(min(lambda)<0) stop("Values in lambda cannot be negative.\n")
  if(!is.null(Cost))if(min(Cost,na.rm=TRUE)<0) warning("There is negative value in Cost.\n")

  if(((Cost.only)|(!is.null(lambda)))&(is.null(Cost))) stop("Cost is required if fitting cost-only or net-benefit regression.\n") 

  if(L<=0) stop("Time limit L must be positive.\n") 
  if(L>max(Followup,na.rm=TRUE)) stop("Time limit L is greater than maximun of follow-up times. Choose a smaller L.\n") 

  if(!is.null(Cost)){
  	if(is.vector(Cost))Cost=matrix(Cost,ncol=1)
  	Cost=as.matrix(Cost)
	if(dim(Cost)[1]!=n) stop("Number of rows of Cost is different to length of Followup.\n")
  }
  m=dim(Cost)[2]

  if(!is.null(Patition.times)){
	Patition.times=unlist(Patition.times)
	if(min(Patition.times)<0) stop("Values in Patition.times cannot be negative.\n")
	if(!is.null(Cost)){
	  if(length(Patition.times)!=m){ 
		if(m>1)stop("Length of Patition.times is different to number of columns of Cost.\n") else {
			Patition.times=NULL; warning("Patition.times is not used since only total Cost is provided.\n")}
	  }
	}else m=length(Patition.times)
	if (!all(Patition.times == cummax(Patition.times))) stop("Patition.times are not monontonically increasing.\n")
	if(length(Patition.times)>1)if (Patition.times[length(Patition.times)]<L) stop("Time limit L is greater than last value of Patition.times. Choose a smaller L.\n")
	Patition.times=c(0,Patition.times)	#add time 0
  }

  # if Eff not provided, assume effectiveness is survival 
  if(is.null(Eff)) {
    if((!is.null(lambda))|Eff.only)cat('Note: Eff is not provided. Assume effectiveness is survival time.\n')
    if(!is.null(Patition.times)){
	Eff=matrix(rep(Patition.times[-1]-Patition.times[-(m+1)],each=n),n,m)
	for(i in 1:n)
	{
		if(is.na(Followup[i]))Eff[i,]=0 else
  		for(j in 1:m)
		{
  			# truncate Effect by follow-up time to obtain observed grouped life time
			if(Followup[i]<Patition.times[j+1])
			{
				if (Followup[i]<Patition.times[j]) Eff[i,j]=0 
				else Eff[i,j]=Followup[i]-Patition.times[j]
			}
		}  
  	}
    }else Eff=matrix(Followup,n,1)
  } else
  {
	if(is.vector(Eff))Eff=matrix(Eff,ncol=1)else
	{
		Eff=as.matrix(Eff)
	}
	if(dim(Eff)[1]!=n) stop("Number of rows of Eff is different to length of Followup.\n")
	if(dim(Eff)[2]!=m) stop("Length of Patition.times is different to number of columns of Eff.\n")
  }

  if(Method=='PT')
  {
	if(m==1) {
		cat('Note: Only 1 time interval is provided for PT method. SW method is used.\n')
		Method='SW'
	} #else if((is.null(Eff))&(is.null(Patition.times)))stop("Either Eff or Patition.times is needed for PT method.\n")
  }

  if((is.null(Patition.times))&(m>1)) stop("When cost history is available, Patition.times is required to truncate them by time limit L.\n")

#  if((is.null(Patition.times))&(!is.null(Eff))) cat("Note: Patition.times is not provided. Assume all costs and effectiveness are evenly spreaded until follow-up time when truncated by time limit L.\n")
#  if((is.null(Patition.times))&(is.null(Eff))) cat("Note: Patition.times is not provided. Assume all costs are evenly spreaded until follow-up time when truncated by time limit L.\n")

  if(!is.null(Z)){
	if(is.vector(Z))Z=data.frame(Z)
	Z=as.data.frame(Z) 
	if(dim(Z)[1]!=n) stop("Number of rows of covariates Z is different to length of Followup.\n")
	Z=Z%>% mutate(across(where(is.character), as.factor))	#change those character as factor
  }
  

  ###### truncate data by L ####

  # truncate follow-up time & death indicator
  deltaL=delta
  FollowupL=Followup
  deltaL[Followup>=L]=1
  FollowupL[Followup>=L]=L

  # truncate Cost and Eff. Assume Cost and Eff evenly spreaded in each time interval
  if(m==1){		#only total cost/eff are provided, assume they spreaded in follow-up time
	if(!is.null(Cost))Cost[,1]=Cost[,1]*FollowupL/Followup
	Eff[,1]=Eff[,1]*FollowupL/Followup
  }else
  {
  	for(i in 1:n)
	{
  		for(j in 1:m)
		{
  			# truncate Cost and Effect by follow-up time
			if(!is.na(FollowupL[i]))
			if(FollowupL[i]<Patition.times[j+1])
			{
				if (FollowupL[i]<=Patition.times[j]) {Eff[i,j]=0;if(!is.null(Cost))Cost[i,j]=0}
				else {	#prorate
				  Eff[i,j]=Eff[i,j]*(FollowupL[i]-Patition.times[j])/(min(Followup[i],Patition.times[j+1])-Patition.times[j])	
				  if(!is.null(Cost))Cost[i,j]=Cost[i,j]*(FollowupL[i]-Patition.times[j])/(min(Followup[i],Patition.times[j+1])-Patition.times[j])	
				}
			}
		}  
  	}
	
	# also revise intervals if needed
	if(Patition.times[m+1]>L)
	{
		mL=min(which(Patition.times>=L))-1
		Patition.times=c(Patition.times[1:mL],L)
		if(!is.null(Cost))Cost=Cost[,1:mL]
		Eff=Eff[,1:mL]
		if(mL==1){
			if(!is.null(Cost))Cost=matrix(Cost,ncol=1)
			Eff=matrix(Eff,ncol=1)
  			if(Method=='PT' & m>1){
				cat('Note: Only 1 time interval within time limit L is provided for PT method. SW method is used.\n')
				Method='SW'
			}
		}
		m=mL
	}
  }

  cat("Time limit horizon L =",L,"\n")
  cat("Censoring rate =",round(1-mean(deltaL,na.rm =TRUE),3)*100,"%\n")

  if(sum(1-deltaL,na.rm =TRUE)==0) {
	cat('Note: Data is not censored within time limit L. All methods are equivalent to OLS.\n')
	if(!Doubly.Robust) Method='CC' else Method='SW' 
  }

  if(Doubly.Robust){
	if(!is.null(PS.Z)){
		if(is.vector(PS.Z))PS.Z=matrix(PS.Z,ncol=1)
		PS.Z=as.matrix(PS.Z) 
		if(dim(PS.Z)[1]!=n) stop("Number of rows of covariates PS.Z is different to length of Followup.\n")
	}
	if((Method!='SW')&(Method!='PT')) stop("SW or PT method is required for doubly robust method.\n")
  }

  ###### estimate

  if (is.null(lambda))n.lam=0 else n.lam=length(lambda)
  if(Eff.only)n.lam=n.lam+1
  if(Cost.only)n.lam=n.lam+1
  results <- as.list(1:n.lam)

  if(Method=="CC"){

    if(sum(1-deltaL,na.rm =TRUE)==0)cat("Method: Ordinary Least Squares\n")
    else {
	cat("Method: Complete Case Only\n")
	cat('Note: Complete case only method is not recommended for censored data due to bias and efficiency loss.\n')
    }
    cat("\n")  
   
    Eff.total=apply(Eff,1,sum)
    if((!is.null(lambda))|Cost.only)Cost.total=apply(Cost,1,sum)
    if(!is.null(lambda)){
	for(j in 1:length(lambda))
	{
		fit=nbreg.CC(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,
			interaction=interaction,lambda=lambda[j])
		results[[j]]<-fit
  	}
    }
    if(Eff.only){
	fit=nbreg.CC(X=FollowupL,delta=deltaL,Cost.total=rep(0,n),Eff.total=Eff.total,group=group,Z=Z,interaction=interaction)
	results[[length(lambda)+1]]<-fit
    }
    if(Cost.only){
	fit=nbreg.CC(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,interaction=interaction,Cost.only=TRUE)
	results[[n.lam]]<-fit
    }
  }else
  if(Method=="AL"){

    if(sum(1-deltaL,na.rm =TRUE)==0)cat("Method: Ordinary Least Squares\n")
    else {
	cat("Method: All Data Ignoring Censoring\n")
    	cat('Note: All data method is not recommended for censored data due to bias.\n')
    }
    cat("\n")  

    Eff.total=apply(Eff,1,sum)
    if((!is.null(lambda))|Cost.only)Cost.total=apply(Cost,1,sum)
    if(!is.null(lambda)){
	for(j in 1:length(lambda))
	{
		fit=nbreg.AL(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,
			interaction=interaction,lambda=lambda[j])
		results[[j]]<-fit
  	}
    }
    if(Eff.only){
	fit=nbreg.AL(X=FollowupL,delta=deltaL,Cost.total=rep(0,n),Eff.total=Eff.total,group=group,Z=Z,interaction=interaction)
	results[[n.lam]]<-fit
    }
    if(Cost.only){
	fit=nbreg.AL(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,interaction=interaction,Cost.only=TRUE)
	results[[n.lam]]<-fit
    }
  }else
  if(Method=="SW"){
    if(!Doubly.Robust){
    if(sum(1-deltaL,na.rm =TRUE)==0)cat("Method: Ordinary Least Squares\n\n")
    else cat("Method: Simple Weighted\n\n")

	Eff.total=apply(Eff,1,sum)
      if((!is.null(lambda))|Cost.only)Cost.total=apply(Cost,1,sum)
	if(!is.null(lambda)){
		for(j in 1:length(lambda))
		{
			fit=nbreg.SW(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,
				interaction=interaction,lambda=lambda[j],Sep.K=Sep.K)
			results[[j]]<-fit
  		}
    	}
    	if(Eff.only){
		fit=nbreg.SW(X=FollowupL,delta=deltaL,Cost.total=rep(0,n),Eff.total=Eff.total,group=group,Z=Z,interaction=interaction,Sep.K=Sep.K)
		results[[n.lam]]<-fit
    	}
      if(Cost.only){
		fit=nbreg.SW(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,interaction=interaction,Sep.K=Sep.K,Cost.only=TRUE)
		results[[n.lam]]<-fit
      }
    }else{
	if(sum(1-deltaL,na.rm =TRUE)==0)cat("Method: Doubly Robust Ordinary Least Squares (estimate is causal average INB)\n")
      else cat("Method: Doubly Robust Simple Weighted (estimate is causal average INB)\n")

	if(is.null(PS.Z))cat("Note: PS.Z is not provided and propensity scores are estimated by unadjusted logistic regression.\n")
	if(is.null(Z))cat("Note: Z is not provided and simple net-benefit regressions are fitted.\n")
	cat("\n")

	Eff.total=apply(Eff,1,sum)
      if((!is.null(lambda))|Cost.only)Cost.total=apply(Cost,1,sum)
	if(!is.null(lambda)){
		for(j in 1:length(lambda))
		{
			fit=nbreg.SW.DR(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,
				PS.Z=PS.Z,interaction=interaction,lambda=lambda[j],Sep.K=Sep.K,PS.trim=PS.trim,DR.Reg.SE=DR.Reg.SE)
			results[[j]]<-fit
  		}
    	}
    	if(Eff.only){
		fit=nbreg.SW.DR(X=FollowupL,delta=deltaL,Cost.total=rep(0,n),Eff.total=Eff.total,group=group,Z=Z,
			PS.Z=PS.Z,interaction=interaction,Sep.K=Sep.K,PS.trim=PS.trim,DR.Reg.SE=DR.Reg.SE)
		results[[n.lam]]<-fit
    	}
      if(Cost.only){
		fit=nbreg.SW.DR(X=FollowupL,delta=deltaL,Cost.total=Cost.total,Eff.total=Eff.total,group=group,Z=Z,
			PS.Z=PS.Z,interaction=interaction,Sep.K=Sep.K,PS.trim=PS.trim,DR.Reg.SE=DR.Reg.SE,Cost.only=TRUE)
		results[[n.lam]]<-fit
      }
    }
  }else
  if(Method=="PT"){
    if(!Doubly.Robust){
	if(sum(1-deltaL,na.rm =TRUE)==0)cat("Method: Ordinary Least Squares\n\n")
      else cat("Method: Partitioned\n\n")

	if(!is.null(lambda)){
		for(j in 1:length(lambda))
		{
			fit=nbreg.PT(X=FollowupL,delta=deltaL,Cost.grp=Cost,Eff.grp=Eff,Patition.times=Patition.times,group=group,Z=Z,
				interaction=interaction,lambda=lambda[j],Sep.K=Sep.K)
			results[[j]]<-fit
  		}
    	}
    	if(Eff.only){
		fit=nbreg.PT(X=FollowupL,delta=deltaL,Cost.grp=matrix(0,n,m),Eff=Eff,Patition.times=Patition.times,group=group,Z=Z,
			interaction=interaction,Sep.K=Sep.K)
		results[[n.lam]]<-fit
    	}
      if(Cost.only){
		fit=nbreg.PT(X=FollowupL,delta=deltaL,Cost.grp=Cost,Eff.grp=Eff,Patition.times=Patition.times,group=group,Z=Z,
			interaction=interaction,Sep.K=Sep.K,Cost.only=TRUE)
		results[[n.lam]]<-fit
      }
    }else{
	if(sum(1-deltaL,na.rm =TRUE)==0)cat("Method: Doubly Robust Ordinary Least Squares (estimate is causal average INB)\n")
      else cat("Method: Doubly Robust Partitioned (estimate is causal average INB)\n")

	if(is.null(PS.Z))cat("Note: PS.Z is not provided and propensity scores are estimated by unadjusted logistic regression.\n")
	if(is.null(Z))cat("Note: Z is not provided and simple net-benefit regressions are fitted.\n")
	cat("\n")

	if(!is.null(lambda)){
		for(j in 1:length(lambda))
		{
			fit=nbreg.PT.DR(X=FollowupL,delta=deltaL,Cost.grp=Cost,Eff.grp=Eff,Patition.times=Patition.times,group=group,Z=Z,
				PS.Z=PS.Z,interaction=interaction,lambda=lambda[j],Sep.K=Sep.K,PS.trim=PS.trim,DR.Reg.SE=DR.Reg.SE)
			results[[j]]<-fit
  		}
    	}
    	if(Eff.only){
		fit=nbreg.PT.DR(X=FollowupL,delta=deltaL,Cost.grp=matrix(0,n,m),Eff.grp=Eff,Patition.times=Patition.times,group=group,Z=Z,
			PS.Z=PS.Z,interaction=interaction,Sep.K=Sep.K,PS.trim=PS.trim,DR.Reg.SE=DR.Reg.SE)
		results[[n.lam]]<-fit
    	}
      if(Cost.only){
		fit=nbreg.PT.DR(X=FollowupL,delta=deltaL,Cost.grp=Cost,Eff.grp=Eff,Patition.times=Patition.times,group=group,Z=Z,
			PS.Z=PS.Z,interaction=interaction,Sep.K=Sep.K,PS.trim=PS.trim,DR.Reg.SE=DR.Reg.SE,Cost.only=TRUE)
		results[[n.lam]]<-fit
      }
    }
  }

  if(!Doubly.Robust) class(results) <- "NetBenefitReg" else class(results) <- "DRNetBenefitReg" 
  return(results)
}


########### function to create CEAC plot based on fitted models ############ 

#### for DR method, CEAC is for marginal cost-effectiveness

plot.DRNetBenefitReg<-function(object, add=FALSE,xlab = "WTP",ylab="Probability Treatment 1 is cost-effective",pch=20,cex=1,lty=1,lwd=1,type='o', ...)
# object = a fitted "DR.NetBenefitReg" model object
# add=TRUE to add the curve instead of creating a new plot
# other parameters and ... = other graphical parameters for the plot
{

  n.lam=length(object)
  if(is.na(object[[n.lam]]$lambda)) n.lam=n.lam-1	
  if(is.na(object[[n.lam]]$lambda)) n.lam=n.lam-1	#one for Cost.only and one for Eff.only
  if(n.lam<1)stop("No lambda provided.\n")

  lambda=CEAC=rep(NA,n.lam)
  for(i in 1:n.lam)
  {
	res=object[[i]]
	lambda[i]=res$lambda
	CEAC[i]=res$CEAC
  }

  if(n.lam==1) warning("One 1 value in lambda and CEAC is a point only.\n")

  if(!add){
	plot(lambda,CEAC,ylim = c(0,1), pch=pch,cex = cex,
             xlab = xlab, ylab=ylab,type=type,lty=lty,lwd=lwd,...)
  }else {
	lines(lambda,CEAC,lty=lty,lwd=lwd,type=type,pch=pch,cex = cex,...)
  }

#  return(data.frame(lambda,CEAC))	#also return CEAC values
}


#### for non-DR regression method, CEAC can be for either same or heterougenous cost-effectiveness 

plot.NetBenefitReg<-function(object, subgroup=NULL, add=FALSE,xlab = "WTP",ylab="Probability Treatment 1 is cost-effective",pch=20,cex=1,lty=1,lwd=1,type='o', ...)
# object = a fitted "NetBenefitReg" model object
# subgroup = a list of covariates to define subgroup, if there is interaction in NBR, values for those covariates in interaction terms are required 
# add=TRUE to add the curve instead of creating a new plot
# other parameters and ... = other graphical parameters for the plot
{

  n.lam=length(object)
  if(is.na(object[[n.lam]]$lambda)) n.lam=n.lam-1
  if(is.na(object[[n.lam]]$lambda)) n.lam=n.lam-1	#one for Cost.only and one for Eff.only
  if(n.lam<1)stop("No lambda provided.\n")
  if(n.lam<2) warning('One 1 value in lambda and CEAC is a point only.\n')

  lambda=CEAC=rep(NA,n.lam)

  if(length(object[[1]]$int.name)==0){  #no interaction
	for(i in 1:n.lam)
	{
		res=object[[i]]
		lambda[i]=res$lambda
		CEAC[i]=res$CEAC
	}
  }else{	#there is interaction
	if(is.null(subgroup)) stop("subgroup is required since there is interaction between treatment and covariates.\n")
	if(!is.list(subgroup)) stop("subgroup must be a list for covariate values.\n")

	# check whether subgroup Covariate is correct
	n.sub=length(subgroup)
	sub1=data.frame(subgroup)
	if(dim(sub1)[1]!=1) stop("Only one value should be provided for each covariate in subgroup.")

	sub0=object[[1]]$covar1st
	int1=object[[1]]$int.name
	int2=intersect(names(sub1),int1) 
	if(!identical(sort(int1),sort(int2))) stop("Values for all covariates in interactions need to be provided.\n")
	for(j in 1:length(int2)){
		if(is.factor(sub0[[int1[j]]]))if(!(sub1[[int1[j]]][1]%in%levels(sub0[[int1[j]]]))) 	stop("Invalid factor level in subgroup\n.") 
		sub0[[int1[j]]][1]=sub1[[int1[j]]][1]
	}
	sub3<-c(model.matrix(as.formula(paste("~",paste(int1,collapse = "+"),sep="")),data=sub0))		#get design matrix for interaction
	p=length(object[[1]]$est)	#2nd is for group, last few ones are for interactions

	for(i in 1:n.lam)
	{
		res=object[[i]]
		lambda[i]=res$lambda
		est=c(res$est)
		cov_mat=res$covariance
		est.subgroup=est[c(2,(p-length(sub3)+2):p)]%*%sub3
		se.subgroup=sqrt(sub3%*%cov_mat[c(2,(p-length(sub3)+2):p),c(2,(p-length(sub3)+2):p)]%*%sub3)		
		CEAC[i]=pnorm(est.subgroup/se.subgroup)
	}
  }

  if(!add){
	plot(lambda,CEAC,ylim = c(0,1), pch=pch,cex = cex,
             xlab = xlab, ylab=ylab,type=type,lty=lty,lwd=lwd,...)
  }else {
	lines(lambda,CEAC,lty=lty,lwd=lwd,type=type,pch=pch,cex = cex,...)
  }

#  return(data.frame(lambda,CEAC))	#also return CEAC values
}

