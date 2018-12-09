###########################################################################
#initialisation
###########################################################################

rm(list=ls())    #empty workspace if wanted ...
ls()
library(MASS)
library(AlgDesign)

#preparation
dimdoe_path<-"C:\\Dokumente und Einstellungen\\orth\\Desktop\\HIS-1718-wise\\"
script_path<-paste(dimdoe_path,"R-sequences\\",sep="")
ex_path<-paste(dimdoe_path,"examples\\UE2_defoamer17\\",sep="")

#transferring functions
source(paste(script_path,"dim_doe_functions_dataList_Vin.R",sep=""))

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
#establish connection between u and x
VsAndWs<-supplyVsAndWs(infosFromU$V,infosFromU$u_roles)
print(VsAndWs$flag)
if(VsAndWs$flag=="ok"){
  print(cbind(infosFromU$u_roles[infosFromU$u_roles!="const"],VsAndWs$cosines))
  print(VsAndWs$VRES)
  print(VsAndWs$WRES)
}

#setup low and high for x - not needed when u-design is imported (variant 2)
DLFsettings<-fetchDLFSettings(VsAndWs$VRES,infosFromU$u_low,infosFromU$u_high)  #x-levels from u-levels

#########################################################################
#generate design and model import response values
#########################################################################

#variant 2: importing an existing RT-design
u_design<-provide_u_design(dataList,u_roles,u_transform,u_gradient,u_offset)
design<-fetchDesign(2,VsAndWs$VRES,VsAndWs$WRES,u_design,infosFromU$u_low,infosFromU$u_high,DLFsettings$x_low,DLFsettings$x_high,infosFromU$u_roles,VsAndWs$x_roles,test_eqpr,desType,desParameters)
print(design$flag)
if(design$flag=="ok"){
  print(10^design$u_design)
  print(design$u_design%*%VsAndWs$VRES)
  print(cbind(design$x_D_low,design$x_D_high))
  print(design$MD)
  print(design$Dfrml)
  print(design$dimDoeMod)
}

#Importing Response-names and values 
V_resp<-fetchV_resp(dataList$Z,V,import=1)
z<-dataList$z_values 
infosFromZ<-fetchInfosFromZ(dataList$Z)


##########################################################################
# fit, diagnose and use the model
##########################################################################

#define the linear (default for imported designs) or interaction model
Dfrml<-design$Dfrml
Dfrml<-paste(design$Dfrml,"I(Fr*cT)",sep="+")
Dfrml

#Fitting and Diagnosing the fitted model - this uses linear modelling lm( )
fittedModel<-analyseDesign(Dfrml,design$u_design,z,infosFromZ$z_transform,infosFromZ$z_gradient,infosFromZ$z_offset,V,V_resp,design$MD,debug=0)

#model coefficients
lapply(fittedModel,function(x) summary(x)$coeff)

#diagnostics using residual analysis y-response
lapply(fittedModel,plot)

#R2 and RSD	 - regression measure and residual standard deviation for DL-response	
R2Y<-lapply(fittedModel,function(x) summary(x)$r.square)
R2Y
RSDY<-lapply(fittedModel,function(x) summary(x)$sigma)
RSDY
RSDinPercent<-lapply(fittedModel,function(x) (10^(summary(x)$sigma)-1)*100)
RSDinPercent

#observed vs predicted for DL-response
y_pred<-lapply(fittedModel, fitted)		#as a list of vectors
y_predM<-z		#this is a 2D-matrix with the correct size - it will be completely overwritten
for(colInd in (1:length(y_pred))) {
  y_predM[,colInd]<-y_pred[[colInd]]
  colnames(y_predM)<-names(y_pred)
}

y_res<-lapply(fittedModel, residuals)		#as a list of vectors
y_resM<-z		#this is a 2D-matrix with the correct size - it will be completely overwritten
for(colInd in (1:length(y_res))) y_resM[,colInd]<-y_res[[colInd]]
y_obsM=y_predM+y_resM
for(i in 1:length(z[1,])){
  dev.new()
  plot(y_predM[,i], y_obsM[,i],main=paste("Obs vs Pred for DL-resp: ",colnames(y_predM)[i],sep=""))
  abline(0,1)
  text(y_predM[,i], y_obsM[,i], row.names(y_obsM), cex=0.8, pos=4)
}

#observed vs predicted for RT-response
l_z_pred<-y_predM%*%solve(V_resp[(nrow(VsAndWs$VRES)+1):nrow(V_resp),])-design$u_design%*%V_resp[1:nrow(VsAndWs$VRES),]
z_pred<-10^l_z_pred
z_pred

R2Z<-diag(cor(z,z_pred))^2
R2Z
for(i in 1:length(z[1,])){
  dev.new()
  plot(z_pred[,i],z[,i],main=paste("Obs vs Pred for RT-resp: ",colnames(z)[i],sep=""))
  abline(0,1)
  text(z_pred[,i],z[,i], row.names(z), cex=0.8, pos=4)
}


#################################################################
#  Making predictions (e.g. for Scale-Up or for model validation)
#################################################################

u_new<-uu2u(dataList$uu_new,u_transform,u_gradient,u_offset)
scaleUp<-scaleUpProcess(u_new$u_vec,dataList$z_new,fittedModel,infosFromU$u_roles,infosFromU$u_low,infosFromU$u_high,VsAndWs$VRES,V_resp,design$x_D_low,design$x_D_high,validate=TRUE,debug=0)

#plot obs vs pred for y
for(i in 1:length(z[1,])){
  dev.new()
  plot(scaleUp$y_new[,i]~scaleUp$y_new_pred[,i])
  abline(0,1)
  text(scaleUp$y_new[,i]~scaleUp$y_new_pred[,i], row.names(scaleUp$y_new[,i]), cex=0.8, pos=4)
}

#plot y-res vs x-var Fr-sc
for(i in 1:length(z[1,])){
  d<-10^scaleUp$u_new[,"d"]
  dev.new()
  plot(scaleUp$y_new_res[,i]~scaleUp$x_new_sc[,"Fr"],cex=5*d)
  abline(0,0)
  abline(-0.05,0,col="red")
  abline(0.05,0,col="red")
  names<-as.list(10^scaleUp$u_new[,"d"])
  text(scaleUp$y_new_res[,i]~scaleUp$x_new_sc[,"Fr"], names, cex=0.8, pos=4)
           #labelling is just numbers 1 to 31 instead of d-values
}


#plot obs vs pred for z
for(i in 1:length(z[1,])){
  dev.new()
  plot(scaleUp$z_new[,i]~scaleUp$z_new_pred[,i])
  abline(0,1)
  text(scaleUp$z_new[,i]~scaleUp$z_new_pred[,i], row.names(scaleUp$y_new[,i]), cex=0.8, pos=4)
           #labelling is again just numbers 1 to 31 instead rownames
}
