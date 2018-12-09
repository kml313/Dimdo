rm(list=ls()) #empty workspace if wanted ...
ls()
library(MASS)
library(AlgDesign)

#preparation
dimdoe_path<-"../"
script_path<-paste(dimdoe_path,"lib/",sep="")
ex_path<-paste(dimdoe_path,"data/",sep="")

#transferring functions
#source(paste(script_path,"dim_doe_functions_dataList_Vin.R",sep=""))
source(paste(script_path,"lib.R",sep=""))

data <-myRound()
data
plot(data)
