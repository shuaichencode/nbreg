# nbreg
 ------------------------------------------------------------------------------
DESCRIPTION:    
Sample R code to illustrate how to perform net-benefit regression methods using a simulated censored observational data

REFERENCE:
Chen S, Hoch JS. Net-benefit regression with censored cost-effectiveness data from randomized or observational studies. Stat Med. Jun 03 2022;doi:10.1002/sim.9486
https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9486

CODE AUTHOR:  
Shuai Chen  <shschen@ucdavis.edu>          
Division of Biostatistics, Department of Public Health Sciences, University of California, Davis     
------------------------------------------------------------------------------

FILES:

Censored_CEdata.csv: The simulated censored data
True_CEdata.csv: The true uncensored data, used to evaluate methods
Main_CreateData.r: main program which simulates the data using Data_Gen.r and save the datasets into .csv files.
Main_NBR.r: main program which reads in the .csv datasets and performs net-benefit regression methods using Est_NBR.r.
Data_Gen.r: program for data generation, used by Main_CreateData.r
Est_NBR.r: program for net-benefit regressions and doubly robust methods for censored data, used by Main_NBR.r
