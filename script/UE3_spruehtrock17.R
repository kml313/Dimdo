###########################################################################
#initialisation
###########################################################################

rm(list=ls())    #empty workspace if wanted ...
ls()
library(MASS)
library(AlgDesign)

#preparation
dimdoe_path<-"../"
script_path<-paste(dimdoe_path,"lib/",sep="")
ex_path<-paste(dimdoe_path,"data/UE3_spruehtrock17/",sep="")

#transferring functions
#source(paste(script_path,"dim_doe_functions_dataList_Vin.R",sep=""))
source(paste(script_path,"dim_doe_functions_TurnQual.R",sep=""))



##########################################################################
#start example
##########################################################################
#fetch dataList
dataList<-readDimDoeData(ex_path)

#fetch u_infos
infosFromU<-fetchInfosFromU(dataList$U)
u_roles<-infosFromU$u_roles
u_roles
u_transform<-infosFromU$u_transform
u_transform
u_gradient<-infosFromU$u_gradient
u_gradient
u_offset<-infosFromU$u_offset
u_offset
u_low<-infosFromU$u_low
u_low
u_high<-infosFromU$u_high
D<-infosFromU$D
D
V<-infosFromU$V
V

#E2<-prepareProjection(infosFromU$u_roles,c("contr","scup","const"),transpose=TRUE)

Vok<-suggestVmatrix(u_roles,D,V)
Vok$flag
data.frame(u_roles,Vok$V)
V<-adjustVmatrix(Vok$V,2,0,"")
V

#establish connection between u and x
VsAndWs<-supplyVsAndWs(V,u_roles)
print(VsAndWs$flag)
VRES<-VsAndWs$VRES
VRES
WRES<-VsAndWs$WRES
WRES
x_roles<-VsAndWs$x_roles
x_roles
u_roles<-VsAndWs$u_roles
u_roles

#setup low and high for x - not needed when u-design is imported (variant 2)
DLFsettings<-fetchDLFSettings(VRES,u_low,u_high)  #x-levels from u-levels

test_eqpr<-FALSE
desParam<-list(numberCP=3)
desType="designFractionalFactorial"
desParam$reduce<-1
desParam$vorz<- -1
DLdesign<-fetchDesign(4,VsAndWs$VRES,VsAndWs$WRES,dataList,infosFromU$u_low,infosFromU$u_high,DLFsettings$x_low,DLFsettings$x_high,VsAndWs$u_roles,VsAndWs$x_roles,test_eqpr,desType,desParam)
RTdesign<-fetchDesign(1,VsAndWs$VRES,VsAndWs$WRES,dataList,infosFromU$u_low,infosFromU$u_high,x_D_levels$x_D_low,x_D_levels$x_D_high,VsAndWs$u_roles,VsAndWs$x_roles,test_eqpr,desType,desParam)

DLdesign$flag
RTdesign$flag
10^DLdesign$u_design
10^RTdesign$u_design
plot(DLdesign$u_design[,"MS"]~DLdesign$u_design[,"MG"])
dev.new()
plot(RTdesign$u_design[,"MS"]~RTdesign$u_design[,"MG"])
dev.new()
plot(DLdesign$u_design[,"T"]~DLdesign$u_design[,"p"])
dev.new()
plot(RTdesign$u_design[,"T"]~RTdesign$u_design[,"p"])

DLdesign$u_design%*%VRES
DLdesign$dimDoeMod
transDimDoeMod(DLdesign$dimDoeMod,x_roles,forward=TRUE) #should be the same as design$Dfrm
transDimDoeMod(DLdesign$Dfrm,x_roles,forward=FALSE)	#should be the same as design$$dimDoeMod

#Importing Response-names and values 

infosFromZ<-fetchInfosFromZ(dataList$Z)
V_resp<-infosFromZ$V_resp
z<-dataList$z_values

##########################################################################
# fit, diagnose and use the model - using lm() and Dfrm
##########################################################################

#define the linear (default for imported designs) or interaction model
Dfrml<-as.list(DLdesign$Dfrml)

#Fitting and Diagnosing the fitted model - this uses linear modelling lm( )
fittedModel<-analyseDesign(Dfrml,DLdesign$u_design,z,infosFromZ$z_transform,infosFromZ$z_gradient,infosFromZ$z_offset,V,V_resp,DLdesign$MD,debug=0)
fittedModel

#model coefficients
lapply(fittedModel,function(x) summary(x)$coeff)

#diagnostics using residual analysis y-response
lapply(fittedModel,plot)

#RSD
RSDY<-lapply(fittedModel,function(x) summary(x)$sigma)
RSDY
RSDinPercent<-lapply(fittedModel,function(x) (10^(summary(x)$sigma)-1)*100)
RSDinPercent


