
### Function to generate observational data	#######

# generate time in years 
# assume cost evenly spreaded within each interval (e.g., year) 
# The quality of life was generated randomly between 0.3 to 1 before a heart failure or death event
# and between 0 to 0.3 after a heart failure event before death for patients experiencing heart failure. 
# generate a binary treatment (eg, 1 for new treatment, 0 for comparison)
# include 3 binary covariates (age>=65 vs <65, LBBB vs non-LBBB, female vs male)
# also generate true uncensored outcomes, and potential outcomes (counterfactuals) under both treatment 0 and 1 to help evaluate the methods

require(mvtnorm)
require(survival)

sim.CostEff.Obs<-function(n,Coef.Surv.MainEff,Coef.Surv.InterEff,Par.Censor,Coef.Trt,p=3,k=30,max.time,cov.name=NULL)
#n: sample size
#Coef.Surv.MainEff and Coef.Surv.InterEff: coefficients for T=Survival time until death (following exp distribution)
#Par.Censor: parameter for max of C=censoring time (following uniform distribution [1,Par.Censor])
#Coef.Trt: coefficients for treatment assignment model using logistic regression
#p: number of covariates
#k interval to group outcomes (ie, costs within Day(or Year) 1~k, costs within Day(or Year) (k+1)~2k, etc.)
#max.time: time limit horizon 
#cov.name: names of covariates 
{
  # generate p covariates 
  Sigma=matrix(0.3,p,p)
  diag(Sigma)=1
  Z=rmvnorm(n, mean = rep(0, p), sigma = Sigma)	
  Z[,1]=round(Z[,1]*10+60)	#generate age
  Z[,1]=(Z[,1]>=65)+0	#dichotomize age
  Z[,-1]=(Z[,-1]>0)+0	#dichotomize other covariates

  # generate treatment assignment Trt from logistic model
  PS.lin=cbind(1,Z)%*%Coef.Trt 
  Trt.Prob=exp(PS.lin)/(1+exp(PS.lin))  		
  Trt=rbinom(n,1,prob=Trt.Prob)		#Trt=1 or 0

  # generate survival time T in days, also get potential uncensored outcomes (counterfactuals) used to evaluate methods
  T1=rexp(n,1/exp(cbind(1,Z)%*%Coef.Surv.MainEff+(cbind(1,Z)%*%Coef.Surv.InterEff))) 
  T0=T1/exp(cbind(1,Z)%*%Coef.Surv.InterEff)
  T=T1
  T[Trt==0]=T0[Trt==0]

  # generate censoring time C in days
  C=runif(n)*(Par.Censor-1)+1  #independent censored time 

  L=max.time
  X=delta=TL=TL1=TL0=numeric(n)
  for(i in 1:n)
  {
	TL[i]=min(T[i],L)		#L-truncated survival
	TL1[i]=min(T1[i],L)
	TL0[i]=min(T0[i],L)
	X[i]=min(TL[i],C[i])	#L-truncated follow-up time
	delta[i]=(TL[i]<=C[i])	#L-truncated death indicator
  }
  delta[which(max(X)==X)]=1 # force the last observation to be dead to prevent possibly error due to inversed probability of 1/0 
	#this error can be prevented by reducing L in model fitting
 
#### generate grouped costs, each within a k-time interval, also generate potential costs #####

  t.int=seq(0,L,k)                        #time points partition time 0-L with k intervals
  if (t.int[length(t.int)]<L) t.int=c(t.int,L)
	
  # parameters for fixed costs in a time interval, following log-normal distribution, depending on treatment and LBBB
  par_fix0=rep(7,n)	
  par_fix0[(Z[,2]==1)]=6.6 	
  par_fix1=rep(6,n) 	# If in treatment group, has lower fixed cost
  par_fix1[(Z[,2]==1)]=4.5 	
  Cost_ann1=matrix(rep(rlnorm(n,par_fix1,0.2),length(t.int)-1),n,length(t.int)-1) 
  Cost_ann0=matrix(rep(rlnorm(n,par_fix0,0.2),length(t.int)-1),n,length(t.int)-1) 
  Cost_ran=matrix(rlnorm((length(t.int)-1)*n,4,0.2),n,length(t.int)-1)

  Cost0=Cost_ann0+Cost_ran      
  Cost1=Cost_ann1+Cost_ran 
  
  # parameters for diagnosis costs
  par_diag0=rep(8.5,n)
  par_diag1=9.5	# If in treatment group, has more diagnosis cost
  Cost_diag0=rlnorm(n,par_diag0,0.2)
  Cost_diag1=rlnorm(n,par_diag1,0.2)
  Cost0[,1]=Cost_diag0+Cost0[,1]	#add diagnosis cost into the first time interval
  Cost1[,1]=Cost_diag1+Cost1[,1]
  Cost=Cost0	
  Cost[Trt==1,]=Cost1[Trt==1]

  Cost_term=rlnorm(n,9,0.6)	# terminal cost at death time T

  # add terminal cost, truncate cost within L, and censor cost at the censoring time
  Cost_obs=Cost_true=Cost
  Cost_true0=Cost0
  Cost_true1=Cost1
  for(i in 1:n)
  {
  	for(j in 2:length(t.int))
	{
  		# truncate cost by follow-up time X to obtain observed cost until X
		if(X[i]<=t.int[j])
		{
			if (X[i]<=t.int[j-1]) Cost_obs[i,j-1]=0 
			else 
			{
				Cost_obs[i,j-1]=Cost[i,j-1]*(X[i]-t.int[j-1])/(t.int[j]-t.int[j-1])	
				# prorate cost if X is in middle of an interval 

				#add random terminal cost to the interval where death time T[i] falls
				if(T[i]==X[i])  Cost_obs[i,j-1]=Cost_obs[i,j-1]+Cost_term[i]
			}
		}

  		# truncate cost by truncated survival time TL to obtain true cost until TL
		if(TL[i]<=t.int[j])
		{
			if (TL[i]<=t.int[j-1]) Cost_true[i,j-1]=0 
			else 
			{
				Cost_true[i,j-1]=Cost[i,j-1]*(TL[i]-t.int[j-1])/(t.int[j]-t.int[j-1])	
				# prorate cost if TL is in middle of an interval 

				#add random terminal cost to the interval where death time T[i] falls
				if(T[i]==TL[i])  Cost_true[i,j-1]=Cost_true[i,j-1]+Cost_term[i]
			}
		}

  		# truncate Cost_true1 by truncated survival time TL1 to obtain true potential cost until TL1
		if(TL1[i]<=t.int[j])
		{
			if (TL1[i]<=t.int[j-1]) Cost_true1[i,j-1]=0 
			else 
			{
				Cost_true1[i,j-1]=Cost1[i,j-1]*(TL1[i]-t.int[j-1])/(t.int[j]-t.int[j-1])	
				# prorate cost if TL1 is in middle of an interval 

				#add random terminal cost to the interval where death time T1[i] falls
				if(T1[i]==TL1[i])  Cost_true1[i,j-1]=Cost_true1[i,j-1]+Cost_term[i]
			}
		}

  		# truncate Cost_true0 by truncated survival time TL0 to obtain true potential cost until TL0
		if(TL0[i]<=t.int[j])
		{
			if (TL0[i]<=t.int[j-1]) Cost_true0[i,j-1]=0 
			else 
			{
				Cost_true0[i,j-1]=Cost0[i,j-1]*(TL0[i]-t.int[j-1])/(t.int[j]-t.int[j-1])	
				# prorate cost if TL0 is in middle of an interval 

				#add random terminal cost to the interval where death time T1[i] falls
				if(T0[i]==TL0[i])  Cost_true0[i,j-1]=Cost_true0[i,j-1]+Cost_term[i]
			}
		}
	}  
  }


#### generate grouped qualty-adjuste life time (in days), each within a k interval #####
#quality of life=0 if dead; between (0,0.3)if after a heart failure event; between(0.3,1) if before a heart failure event

  # generate heart failure time T.HF, using similar coef for survival  
  T.HF1=rexp(n,1/exp(0.8*cbind(1,Z)%*%Coef.Surv.MainEff+3*(cbind(1,Z)%*%Coef.Surv.InterEff))) 
  T.HF0=T.HF1/exp(3*cbind(1,Z)%*%Coef.Surv.InterEff)
  T.HF=T.HF1
  T.HF[Trt==0]=T.HF0[Trt==0]

  q_fix=matrix(rep(runif(n)*0.15,length(t.int)-1),n,length(t.int)-1)	#fixed qualty for a person
  q_ran=matrix(runif((length(t.int)-1)*n)*0.15,n,length(t.int)-1)   #random qualty for a person
  q_ran.befHF=matrix(runif((length(t.int)-1)*n)*0.4+0.3,n,length(t.int)-1)   #random qualty for a person before HF

  QALY_obs=QALY_true=QALY_true1=QALY_true0=matrix(0,n,length(t.int)-1)
  for(i in 1:n)
  {
  	for(j in 2:length(t.int))
	{
  		# true uncensored grouped quality-adjusted life time 
		if(TL[i]>=t.int[j])
		{
			QALY_true[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
			if(T.HF[i]>=t.int[j])QALY_true[i,j-1]=QALY_true[i,j-1]+(t.int[j]-t.int[j-1])*q_ran.befHF[i,j-1]
			if(T.HF[i]<t.int[j] & T.HF[i]>=t.int[j-1])QALY_true[i,j-1]=QALY_true[i,j-1]+(T.HF[i]-t.int[j-1])*q_ran.befHF[i,j-1]
		}
		if(TL[i]<t.int[j])
		{
			if (TL[i]<t.int[j-1]) QALY_true[i,j-1]=0 
			else 
			{
				QALY_true[i,j-1]=(TL[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
				if(T.HF[i]>=TL[i])QALY_true[i,j-1]=QALY_true[i,j-1]+(TL[i]-t.int[j-1])*q_ran.befHF[i,j-1]
				if(T.HF[i]<TL[i] & T.HF[i]>=t.int[j-1])QALY_true[i,j-1]=QALY_true[i,j-1]+(T.HF[i]-t.int[j-1])*q_ran.befHF[i,j-1]
			}
		}

  		# true uncensored potential QALY in group 1  
		if(TL1[i]>=t.int[j])
		{
			QALY_true1[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
			if(T.HF1[i]>=t.int[j])QALY_true1[i,j-1]=QALY_true1[i,j-1]+(t.int[j]-t.int[j-1])*q_ran.befHF[i,j-1]
			if(T.HF1[i]<t.int[j] & T.HF1[i]>=t.int[j-1])QALY_true1[i,j-1]=QALY_true1[i,j-1]+(T.HF1[i]-t.int[j-1])*q_ran.befHF[i,j-1]
		}
		if(TL1[i]<t.int[j])
		{
			if (TL1[i]<t.int[j-1]) QALY_true1[i,j-1]=0 
			else 
			{
				QALY_true1[i,j-1]=(TL1[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
				if(T.HF1[i]>=TL1[i])QALY_true1[i,j-1]=QALY_true1[i,j-1]+(TL1[i]-t.int[j-1])*q_ran.befHF[i,j-1]
				if(T.HF1[i]<TL1[i] & T.HF1[i]>=t.int[j-1])QALY_true1[i,j-1]=QALY_true1[i,j-1]+(T.HF1[i]-t.int[j-1])*q_ran.befHF[i,j-1]
			}
		}

  		# true uncensored potential QALY in group 0 
		if(TL0[i]>=t.int[j])
		{
			QALY_true0[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
			if(T.HF0[i]>=t.int[j])QALY_true0[i,j-1]=QALY_true0[i,j-1]+(t.int[j]-t.int[j-1])*q_ran.befHF[i,j-1]
			if(T.HF0[i]<t.int[j] & T.HF0[i]>=t.int[j-1])QALY_true0[i,j-1]=QALY_true0[i,j-1]+(T.HF0[i]-t.int[j-1])*q_ran.befHF[i,j-1]
		}
		if(TL0[i]<t.int[j])
		{
			if (TL0[i]<t.int[j-1]) QALY_true0[i,j-1]=0 
			else 
			{
				QALY_true0[i,j-1]=(TL0[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
				if(T.HF0[i]>=TL0[i])QALY_true0[i,j-1]=QALY_true0[i,j-1]+(TL0[i]-t.int[j-1])*q_ran.befHF[i,j-1]
				if(T.HF0[i]<TL0[i] & T.HF0[i]>=t.int[j-1])QALY_true0[i,j-1]=QALY_true0[i,j-1]+(T.HF0[i]-t.int[j-1])*q_ran.befHF[i,j-1]
			}
		}

  		# observed grouped quality-adjusted life time until follow-up time X
		if(X[i]>=t.int[j])
		{
			QALY_obs[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
			if(T.HF[i]>=t.int[j])QALY_obs[i,j-1]=QALY_obs[i,j-1]+(t.int[j]-t.int[j-1])*q_ran.befHF[i,j-1]
			if(T.HF[i]<t.int[j] & T.HF[i]>=t.int[j-1])QALY_obs[i,j-1]=QALY_obs[i,j-1]+(T.HF[i]-t.int[j-1])*q_ran.befHF[i,j-1]
		}
		if(X[i]<t.int[j])
		{
			if (X[i]<t.int[j-1]) QALY_obs[i,j-1]=0 
			else 
			{
				QALY_obs[i,j-1]=(X[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
				if(T.HF[i]>=X[i])QALY_obs[i,j-1]=QALY_obs[i,j-1]+(X[i]-t.int[j-1])*q_ran.befHF[i,j-1]
				if(T.HF[i]<X[i] & T.HF[i]>=t.int[j-1])QALY_obs[i,j-1]=QALY_obs[i,j-1]+(T.HF[i]-t.int[j-1])*q_ran.befHF[i,j-1]
			}
		}
	}  
  }  

  Z=data.frame(Z)
  if(!is.null(cov.name))names(Z)=cov.name	
	
  return(list(Patition.times=t.int[-1],Follow.up=X,Covariate=Z,Death=delta,Trt=Trt,Cost.grp.obs=Cost_obs,QASurv.grp.obs=QALY_obs,	#observed data
	TL=TL,Cost_true=Cost_true,QALY_true=QALY_true,Trt.Prob=Trt.Prob,		#true uncensored values, used for evaluation 
	TL1=TL1,TL0=TL0,Cost_true1=Cost_true1,Cost_true0=Cost_true0,QALY_true1=QALY_true1,QALY_true0=QALY_true0))	# potential outcomes, used for evaluation 
}




