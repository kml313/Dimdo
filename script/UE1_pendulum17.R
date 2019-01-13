###########################################################################
#initialisation
###########################################################################

rm(list=ls())    #empty workspace if wanted ...
ls()
library(MASS)
library(AlgDesign)

#preparation
dimdoe_path<-"../"
script_path<-"/home/kamal/Desktop/R-Jenkins/lib/"
#script_path<-paste(dimdoe_path,"lib/",sep="")
ex_path<-"/home/kamal/Desktop/R-Jenkins/data/UE1_pendulum17/"
#ex_path<-paste(dimdoe_path,"data/UE3_spruehtrock17/",sep="")

#transferring functions
#source(paste(script_path,"dim_doe_functions_dataList_Vin.R",sep=""))
#source(paste(script_path,"dim_doe_functions_TurnQual.R",sep=""))
source(paste(script_path,"lib_pendulum.R",sep=""))



##########################################################################
#start example
##########################################################################
#fetch dataList
dataList<-readDimDoeData(ex_path)
names(dataList)
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

Vok<-suggestVmatrix(u_roles,D,V,debug=1)
Vok$flag
Vok$V
V<-adjustVmatrix(Vok$V,1,1,"tanA")


#establish connection between u and x
VsAndWs<-supplyVsAndWs(V,u_roles)
print(VsAndWs$flag)
VSL<-VsAndWs$VSL
VSL
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
x_D_levels<-fetchXLevels(dataList$x_levels,DLFsettings$x_low,DLFsettings$x_high,eps=1E-6)

#variant 3: import a DL-design - given by the matrix MD
test_eqpr<-TRUE
design<-fetchDesign(3,VsAndWs$VRES,VsAndWs$WRES,dataList$MD,infosFromU$u_low,infosFromU$u_high,x_D_levels$x_D_low,x_D_levels$x_D_high,VsAndWs$u_roles,VsAndWs$x_roles,test_eqpr,desType,desParam)
design$flag
u2uu(design$u_design,u_transform,u_gradient,u_offset)
design$u_design%*%VRES
design$dimDoeMod

#variant 4: generate a RT-design
desParam<-list(numberCP=1)
desType="designFullFactorial"
design<-fetchDesign(4,VsAndWs$VRES,VsAndWs$WRES,NULL,infosFromU$u_low,infosFromU$u_high,x_D_levels$x_D_low,x_D_levels$x_D_high,VsAndWs$u_roles,VsAndWs$x_roles,test_eqpr,desType,desParam)
design$flag
u2uu(design$u_design,u_transform,u_gradient,u_offset)
design$u_design%*%VRES
design$dimDoeMod


#variant 1: generate a RT-design
desParam<-list(numberCP=1)
desType="designFullFactorial"
design<-fetchDesign(1,VsAndWs$VRES,VsAndWs$WRES,NULL,infosFromU$u_low,infosFromU$u_high,x_D_levels$x_D_low,x_D_levels$x_D_high,VsAndWs$u_roles,VsAndWs$x_roles,test_eqpr,desType,desParam)
design$flag
u2uu(design$u_design,u_transform,u_gradient,u_offset)
design$u_design%*%VRES
design$dimDoeMod


transDimDoeMod(design$dimDoeMod,x_roles,forward=TRUE) #should be the same as design$Dfrm
transDimDoeMod(design$Dfrm,x_roles,forward=FALSE)	#should be the same as design$$dimDoeMod

#Importing Response-names and values 
V_resp<-fetchV_resp(dataList$Z,V,import=1)
z<-dataList$z_values
infosFromZ<-fetchInfosFromZ(dataList$Z)


##########################################################################
# fit, diagnose and use the model - using lm() and Dfrm
##########################################################################

#define the linear (default for imported designs) or interaction model
Dfrml<-as.list(design$Dfrml)

#Fitting and Diagnosing the fitted model - this uses linear modelling lm( )
fittedModel<-analyseDesign(Dfrml,design$u_design,z,infosFromZ$z_transform,infosFromZ$z_gradient,infosFromZ$z_offset,V,V_resp,design$MD,debug=0)
fittedModel

#model coefficients
lapply(fittedModel,function(x) summary(x)$coeff)

#dianostics using residual analysis y-response
lapply(fittedModel,plot)

#RSD
RSDY<-lapply(fittedModel,function(x) summary(x)$sigma)
RSDY
RSDinPercent<-lapply(fittedModel,function(x) (10^(summary(x)$sigma)-1)*100)
RSDinPercent


#################################################################
#  Making predictions (e.g. for Scale-Up or for model validation)
#################################################################

u_new<-uu2u(dataList$uu_new,u_transform,u_gradient,u_offset)
scaleUp<-scaleUpProcess(u_new$u_vec,dataList$z_new,fittedModel,infosFromU$u_roles,infosFromU$u_low,infosFromU$u_high,VsAndWs$VRES,V_resp,design$x_D_low,design$x_D_high,validate=FALSE,debug=0)

10^scaleUp$u_new
scaleUp$z_new_pred


##########################################################################
# fit, diagnose and use the model - using matix multiplication and inversion
##########################################################################

fittedModel2<-fitAndXvalDesign(design$dimDoeMod,design$u_design,z,infosFromU$V,V_resp,design$MD,debug=0)
names(fittedModel2)
fittedModel2$b
coef(fittedModel[[1]])		#should be the same result


#################################################################
#  Making predictions (e.g. for Scale-Up or for model validation)
######################################################

