#              Initialisation and Import
rm(list=ls())    #empty workspace if wanted ...

ls()
library(MASS)
dimdoe_path<-"C:\\Dokumente und Einstellungen\\orth\\Desktop\\dim-doe\\"
#dimdoe_path<-"W:\\projekte\\dim-doe\\"
script_path<-paste(dimdoe_path,"R-skripte\\",sep="")
script_path
ex_path<-paste(dimdoe_path,"R-neue-Bspl\\",sep="")
ex_path<-paste(ex_path,"UE4_RelayExample\\",sep="")
ex_path
source(paste(script_path,"dim_doe_functions_TurnQual.R",sep=""))

#problem specification
dataList<-readDimDoeData(ex_path)
names(dataList)
head(dataList$uu_design)
infosFromU<-fetchInfosFromU(dataList$U)
#extracting results:
u_names<-infosFromU$u_names
u_names
u_roles<-infosFromU$u_roles
u_roles
u_transform<-infosFromU$u_transform
u_transform
u_gradient<-infosFromU$u_gradient
u_offset<-infosFromU$u_offset
u_low<-infosFromU$u_low
u_high<-infosFromU$u_high
D<-infosFromU$D
importedV<-infosFromU$V
V<-importedV
V
t(V)%*%D

Vok<-suggestVmatrix(u_roles,infosFromU$D,V)
Vok$flag
Vok$V
V<-Vok$V
data.frame(u_roles,Vok$V)

VsAndWs<-supplyVsAndWs(V,u_roles)
flag<-VsAndWs$flag
flag
VRES<-VsAndWs$VRES
x_roles<-VsAndWs$x_roles
data.frame(x_roles,t(VRES)%*%VRES)

#V<-adjustVmatrix(V,9,0,"",debug=0)
#V<-adjustVmatrix(V,8,0,"",debug=0)
#V<-adjustVmatrix(V,7,c(0,0,0,0,0,0,-1,0),"TVerh",debug=0) 
#V<-adjustVmatrix(V,,0,"",debug=0) 

data.frame(u_roles,V)
t(D)%*%V
#write.table(data.frame(u_roles,V),"clipboard",sep="\t",col.names=TRUE,dec=",",na="",row.names=TRUE)

R0<-makeRelayMatrix(dataList$relayDesign,u_roles, u_offset, u_gradient, u_transform,dataList$relayCoeffs,fromDesign=FALSE)$relayMatrix
R0

relayDesign1<-dataList$relayDesign1
relayDesign2<-dataList$relayDesign2
relayCoeffs<-dataList$relayCoeffs


# 0. relay
relay<-relayVmatrix(u_names,u_roles,u_transform,u_gradient,u_offset,u_low,u_high,V,R=R0)
#relay<-relayVmatrix(u_roles,u_transform,u_gradient,u_offset,u_low,u_high,V,R=dataList$R)
VR0<-relay$V
R0<-relay$R
u0_names<-relay$new_names
u0_roles<-relay$new_roles
u0_transform<-relay$new_transform
u0_gradient<-relay$new_gradient
u0_offset<-relay$new_offset
u0_low<-relay$new_low
u0_high<-relay$new_high
u0_mean<-(u0_low+u0_high)/2
10^cbind(u0_low,u0_high)
data.frame(u0_roles,VR0)

#write.table(data.frame(u0_roles,VR),"clipboard",sep="\t",col.names=TRUE,dec=",",na="",row.names=TRUE)


# 1. relay
R00<-makeRelayMatrix(dataList$relayDesign1,u0_roles, u0_offset, u0_gradient, u0_transform)$relayMatrix
qual<-makeRelayMatrix(dataList$relayDesign1,u0_roles, u0_offset, u0_gradient, u0_transform)$Rsquare
qual
diag(qual)

relay<-relayVmatrix(u0_names,u0_roles,u0_transform,u0_gradient,u0_offset,u0_low,u0_high,VR0,R=R00)
VR1<-relay$V
R1<-relay$R
u1_names<-relay$new_names
u1_roles<-relay$new_roles
u1_transform<-relay$new_transform
u1_gradient<-relay$new_gradient
u1_offset<-relay$new_offset
u1_low<-relay$new_low
u1_high<-relay$new_high
cbind(u1_low,u1_high)
data.frame(u1_roles,VR1)

# attempt at 2. relay
flag<-makeRelayMatrix(dataList$relayDesign2,u1_roles, u1_offset, u1_gradient, u1_transform)$flag
flag

#before the third relay qualRel has to be read, in order to have
#variation in W_dia and W_dis
u1_roles["H_w",1]<-"dep"
u1_roles["H_l",1]<-"dep"
R01<-makeRelayMatrix(dataList$relayDesign2,u1_roles, u1_offset, u1_gradient, u1_transform)$relayMatrix
R01
qual<-makeRelayMatrix(dataList$relayDesign2,u1_roles, u1_offset, u1_gradient, u1_transform)$Rsquare
diag(qual)

relay<-relayVmatrix(u1_names,u1_roles,u1_transform,u1_gradient,u1_offset,u1_low,u1_high,VR1,R=R01)
VR2<-relay$V
R2<-relay$R
u2_roles<-relay$new_roles
u2_names<-relay$new_names
u2_transform<-relay$new_transform
u2_gradient<-relay$new_gradient
u2_offset<-relay$new_offset
u2_low<-relay$new_low
u2_high<-relay$new_high
cbind(u2_low,u2_high)
data.frame(u2_roles,VR2)



#write.table(data.frame(u2_roles,R2%*%R1%*%R0),"clipboard",sep="\t",col.names=TRUE,dec=",",na="",row.names=TRUE)



VsAndWs<-supplyVsAndWs(VR2,u2_roles)
flag<-VsAndWs$flag
flag
VSL<-VsAndWs$VSL
VSL
VRES<-VsAndWs$VRES
cor(VRES)
#VRES
WRES<-VsAndWs$WRES
#WRES
x_roles<-VsAndWs$x_roles
u2_roles<-VsAndWs$u_roles
data.frame(x_roles,t(VRES)%*%VRES)
VRES%*%WRES
cosines<-VsAndWs$cosines     #for contr: the bigger the better, 
cbind(u1_roles[u1_roles=="contr"|u1_roles=="eqpr"|u1_roles=="quset"],cosines)     

###################################################
#          low and high levels for DL-factors
###################################################

DLFsettings<-fetchDLFSettings(VRES,u2_low,u2_high)
x_low<-DLFsettings$x_low      #careful: these are logs
x_high<-DLFsettings$x_high    #careful: these are logs
x_mean<-DLFsettings$x_mean
10^x_low    
10^x_high
10^x_mean
data.frame(x_roles,x_low,x_high)
data.frame(x_roles,10^x_low,10^x_high)


##################################################################
###### fetch and evaluate the design
##################################################################

test_eqpr<-FALSE
u_d<-provide_u_design(dataList,VR2,u2_names,u2_roles,u2_transform,u2_gradient,u2_offset,u2_low,u2_high,debug=0)
u_des<-as.matrix(u_d$u_design)
names(u_d)
rownames(u_d$VR)
colnames(u_des)
VR2<-u_d$VR
u2_roles<-u_d$u_roles
rownames(u2_roles)
u2_low<-u_d$u_low
u2_high<-u_d$u_high
rownames(u2_high)

VsAndWs<-supplyVsAndWs(VR2,u2_roles)
flag<-VsAndWs$flag
flag
VSL<-VsAndWs$VSL
VSL
VRES<-VsAndWs$VRES
cor(VRES)
#VRES
WRES<-VsAndWs$WRES
#WRES
x_roles<-VsAndWs$x_roles
u2_roles<-VsAndWs$u_roles
data.frame(x_roles,t(VRES)%*%VRES)
VRES%*%WRES


rownames(VRES)
design<-fetchDesign(2,VRES,WRES,u_des,u2_low,u2_high,x_low,x_high,u2_roles,x_roles,test_eqpr,desType,desParam)
flag<-design$flag
flag
design$Dfrml
u_design<-design$u_design
MD<-design$MD
head(MD)
uu_design<-u2uu(u_design,u1_transform,u1_gradient,u1_offset)
uu_design<-uu_design$uu_vec

cbind(colnames(uu_design),x_roles)
colnames(uu_design)[u1_roles=="contr"]

x_design<-u_design%*%VRES
x_design_sl<-u_design%*%as.matrix(VSL)
data.frame(10^x_design,10^x_design_sl)
#write.table(data.frame(x_design[,x_roles=="dldes"]),"clipboard",sep="\t",col.names=TRUE,dec=",",na="",row.names=TRUE)


