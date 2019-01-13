###########################################################################
#initialisation
###########################################################################

setPath<-function(code=0){
  rm(list=ls())    #empty workspace if wanted ...
  ls()
  print(code)
  library(MASS)
  library(AlgDesign)
  #preparation
  #dimdoe_path<-"../"
  dimdoe_path<-"/home/kamal/Desktop/R-Jenkins/"
  #script_path<-"/home/kamal/Desktop/R-Jenkins/lib/"
  script_path<-paste(dimdoe_path,"lib/",sep="")
  #ex_path<-"/home/kamal/Desktop/R-Jenkins/data/UE3_spruehtrock17/"
  switch(
    code,
    ex_path<-paste(dimdoe_path,"data/UE1_pendulum17/",sep=""),
    ex_path<-paste(dimdoe_path,"data/UE2_defoamer17/",sep=""),
    ex_path<-paste(dimdoe_path,"data/UE3_spruehtrock17/",sep=""),
    ex_path<-paste(dimdoe_path,"data/UE4_RelayExample/",sep="")
  )
 # ex_path<-paste(dimdoe_path,"data/"+exampleName+"/",sep="")
  #transferring functions
  #source(paste(script_path,"dim_doe_functions_dataList_Vin.R",sep=""))
  source(paste(script_path,"dim_doe_functions_TurnQual.R",sep=""))
}
