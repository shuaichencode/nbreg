# nbreg
------------------------------------------------------------------------------
DESCRIPTION:    
Sample R code to illustrate how to perform net-benefit regression methods using a simulated censored observational cost-effectiveness data

ARTICLE:    
A Tutorial on Net-benefit Regression for Real World Cost-effectiveness Analysis Using Administrative Data (Shuai Chen, Heejung Bang, and Jeffrey S. Hoch) 

  
------------------------------------------------------------------------------

FILES:    
Censored_CEdata.csv: The simulated censored data    
True_CEdata.csv: The true uncensored data, used to evaluate methods    
Main_CreateData.r: main program which simulates the data using Data_Gen.r and save the datasets into .csv files    
Main_NBR.r: main program which reads in the .csv datasets and performs net-benefit regression methods     
Data_Gen.r: program for data generation, used by Main_CreateData.r    
nbreg.r: program for net-benefit regressions and doubly robust methods for censored data, used by Main_NBR.r    
nbreg_help.pdf: help file about the technical details and explanation of sample R code and output
