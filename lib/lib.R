######################################################################
######################################################################
#     internal service functions - must be available for all functions
######################################################################
######################################################################

#To do's:
#warning in relay, if instead of __N _N is used, otherwise error message is uninteligible
#allow user defined slack variables - also include slack variables in output design
#output u_design as uu_design (i.e. in user-units)
#input u_design in user-units - reasonably awkward solution given - provide_u_design
#... check if this works well for relayed designs
#construct relay-matrix from relay-design - done (I think)
#relay z_transform, z_gradient and z_offset - done (I think)

#just plain picking cols out of a matrix (because R is too stupid to do it)
sCfM<-function(M,cols){
  Msel<-data.frame(M[,cols])
  if(nrow(M)==1)Msel<-t(Msel)
  colnames(Msel)<-colnames(M)[cols]
  rownames(Msel)<-rownames(M)
  return(as.matrix(Msel))
}
#just plain picking rows out of a matrix (because R is too stupid to do it)
sRfM<-function(M,rows){
  Msel<-cbind(M[rows,])
  if(nrow(M)==1 | ncol(M)!=ncol(Msel))Msel<-t(Msel)
  colnames(Msel)<-colnames(M)
  rownames(Msel)<-rownames(M)[rows]
  return(Msel)
}
#function to round to a desired presision
myRound<-function(x,prec){
  x<-x/prec
  x<-round(x,0)
  x<-x*prec
  return(x)
}


#function to calculate x_low/high for DL-factors
fetchDLFSettings<-function(V,u_low,u_high,debug=0){
  #extreme limits (from RF-limits) in DLF-space
  x_high<-matrix(0,ncol(V))
  x_low<-matrix(0,ncol(V))
  x_mean<-matrix(0,ncol(V))
  
  for(i in 1:ncol(V))
  {  
    x_low[i]<-sum(u_low[sign(V[,i])>0]*V[sign(V[,i])>0,i])+sum(u_high[sign(V[,i])<0]*V[sign(V[,i])<0,i])
    x_high[i]<-sum(u_low[sign(V[,i])<0]*V[sign(V[,i])<0,i])+sum(u_high[sign(V[,i])>0]*V[sign(V[,i])>0,i])
    x_mean[i]<-(x_low[i]+x_high[i])/2
  }
  rownames(x_low)<-colnames(V)
  colnames(x_low)<-"DL-calc-low"
  rownames(x_high)<-colnames(V)
  colnames(x_high)<-"DL-calc-high"
  rownames(x_mean)<-colnames(V)
  colnames(x_mean)<-"DL-calc-mean"
  
  DLFsettings<-list(x_low=x_low,x_high=x_high,x_mean=x_mean)
  return(DLFsettings)
}

#Function to suggest scup and eqpr factors, these are taken out of cntr-factors
suggestScupFactors<-function(u_roles,cosines,V,W,eps=0.00001,debug=0)
{ 
  flag="ok"
  u_roles[u_roles=="scup"]<-"contr"
  if(sum((W%*%V-(W%*%V)%*%(W%*%V))^2)>eps) flag<-"W is not inverse of V"
  if(debug&&sum((W%*%V-(W%*%V)%*%(W%*%V))^2)>eps) print("W not inverse of V: go back to the drawing board)")
  in_contr<-prepareProjection(u_roles,"contr")
  if(nrow(V)==ncol(V)) {    #dann gibt's keine Moeglichkeiten
    in_scup<-sCfM(in_contr,NULL)
    scupCosines<-cosines
  } else {
    #    in_scup<-sCfM(in_contr,sort(t(in_contr)%*%cosines,index.return=TRUE)$ix[1:(ncol(in_contr))])
    #    scupCosines<-sRfM(cosines,sort(t(in_contr)%*%cosines,index.return=TRUE)$ix[1:(ncol(in_contr))])
    in_scup<-sCfM(in_contr,sort(cosines,index.return=TRUE)$ix[1:(ncol(in_contr))])
    scupCosines<-sRfM(cosines,sort(cosines,index.return=TRUE)$ix[1:(ncol(in_contr))])
    in_scup<-sCfM(in_scup,scupCosines!=1)
    scupCosines<-sRfM(scupCosines,scupCosines!=1)
  }
  ScupFactors<-list(in_scup=in_scup,scupCosines=scupCosines,flag=flag)  
  return(ScupFactors)
}


#function ccW for calculating centres and weights
ccW<-function(low,high,debug)
{
  cp<-t((low+high)/2)    # we want these to be row-vectors
  rownames(cp)<-"l.ctr"
  if(debug) print(cp)
  weights<-as.matrix((high-low)/2)		#scaling weights
  SW<-diag(weights[,1],nrow=length(weights))
  rownames(SW)<-rownames(high)
  colnames(SW)<-rownames(high)
  if(debug) print(head(SW))
  
  centresWeights<-list(cp=cp,SW=SW)
  return(centresWeights)
}


#function   transferDLFcentreToR - probably deprecated
transferDLFcentreToR<-function(u_R_cp,x_D_cp,V,W,debug)
  #  W%*%V must be id in DL-space
{
  #Calculating the design centre point to be used in RF-space
  u_D_cp<-u_R_cp+(x_D_cp-u_R_cp%*%V)%*%W              # careful: V%*%W != id in R-space
  if(debug) print(u_D_cp)
  if(debug) print(u_D_cp%*%V)                         # check: should be x_D_cp
  
  return(u_D_cp)
}

#recalculating high low levels from any design (i.e. u- or x-design)
#to be used, when (and after) a design has been imported
recalcLevels<-function(design,debug=0)    #design should be in logs
{
  #redefine low, high
  D_low<-as.matrix(apply(design,2,min))   #rownames automatically set correctly!
  colnames(D_low)<-"l.low"
  D_high<-as.matrix(apply(design,2,max))  #rownames automatically set correctly!
  colnames(D_high)<-"l.high"
  
  newLevels<-list(D_low=D_low,D_high=D_high)
  return(newLevels)
}

#function to re-scale a design
scaleDesign<-function(x_design,x_D_cp,SW)
{
  MD<-(x_design-matrix(rep(x_D_cp,times=nrow(x_design)),ncol=length(x_D_cp),byrow=TRUE))%*% ginv(SW)
  colnames(MD)<-colnames(x_design)
  return(MD)
}

orScaleDesign<-function(design){
  newLevels<-recalcLevels(design)
  D_low<-newLevels$D_low
  D_high<-newLevels$D_high
  D_cp<-ccW(D_low,D_high,0)$cp
  SW<-ccW(D_low,D_high,0)$SW
  MDor<-round(scaleDesign(design,D_cp,SW),5)
  orSD<-list(D_low=D_low,D_cp=D_cp,D_high=D_high,SW=SW,MDor=MDor)
  return(orSD)
}

uvScaleDesign<-function(design){
  D_cp<-apply(design,2,mean)
  SW<-matrix(0,nrow=ncol(design),ncol=ncol(design))
  diag(SW)<-apply(design,2,sd)
  rownames(SW)<-colnames(design)  
  colnames(SW)<-colnames(design)
  MDuv<-round(scaleDesign(design,D_cp,SW),5)
  uvSD<-list(D_cp=D_cp,SW=SW,MDuv=MDuv)
  return(uvSD)
}

turnDesign<-function(design,preDet=0,scale="UV",prec=0.2){
  
  if(preDet[1]==0) diagOrdHead<-NULL else diagOrdHead<-preDet
  #print(preDet)
  if(scale=="UV") {
    ctr_des<-uvScaleDesign(design)$MDuv
    SW<-uvScaleDesign(design)$SW 
    diagOrdTail<-sort(diag(SW),decreasing=TRUE,index.return=TRUE)$ix
  } else if(scale=="OR") {
    ctr_des<-orScaleDesign(design)$MDor 
    SW<-orScaleDesign(design)$SW
    diagOrdTail<-sort(diag(SW),decreasing=TRUE,index.return=TRUE)$ix
  } else {
    ctr_des<-apply(design,2,function(x) x-mean(x))
    diagOrdTail<-sort(diag(t(ctr_des)%*%ctr_des),decreasing=TRUE,index.return=TRUE)$ix
  }
  #rules for using diagOrd: It must only be used for immediate display.
  #                         must not be used to reorder t_eigenM
  #                         must not be used to reorder x_design or t_design or r_design    
  diagOrdTail<-setdiff(diagOrdTail,preDet)
  diagOrd<-c(diagOrdHead,diagOrdTail) 
  #print(diagOrd)
  diagInfo<-diag(t(ctr_des)%*%ctr_des)[diagOrd]
  diagSW<-diag(SW)[diagOrd]
  mv<-apply(design,2,mean)[diagOrd]
  #PCA-turning starts here (could be done iteratively with NIPALS-algorithm)
  if(preDet[1]==0){
    ev<-eigen(t(ctr_des)%*%ctr_des)
    t_eigenV<-ev$values
    t_eigenM<-ev$vectors
  } else {	#add diag(1)-matrix as Eigen-M and diagInfo as Eigen-V
    ev<-eigen(t(ctr_des[,-preDet])%*%ctr_des[,-preDet])
    t_eigenV<-c(diagInfo[1:length(preDet)],ev$values)
    t_eigenM<-matrix(0,ncol=ncol(ctr_des),nrow=ncol(ctr_des))
    t_eigenM[-preDet,]<-cbind(matrix(0,nrow=nrow(ev$vectors),ncol=length(preDet)),ev$vectors)
    t_eigenM[preDet,]<-cbind(diag(1,length(preDet)),matrix(0,ncol=nrow(ev$vectors),nrow=length(preDet)))
  }
  rownames(t_eigenM)<-colnames(design)
  colnames(t_eigenM)<-paste("Comp",1:ncol(t_eigenM),sep="_")#temporary names for components
  for(cn in 1:ncol(t_eigenM)){					#components are named after variable with highest loading
    rn<-which(abs(t_eigenM[,cn])==max(abs(t_eigenM[,cn])))
    if(t_eigenM[rn,cn]<0) t_eigenM[,cn]<- -t_eigenM[,cn]   #so that the most important loading is positive
    colnames(t_eigenM)[cn] <- paste(rownames(t_eigenM)[rn],cn,sep="_") #to seperate comps with the same highest loading
  }
  r_eigenM<-myRound(t_eigenM,prec)
  t_design<-ctr_des%*%t_eigenM
  r_design<-ctr_des%*%r_eigenM
  t_des<-list(diagOrd=diagOrd,mv=mv,diagSW=diagSW,diagInfo=diagInfo,t_eigenV=t_eigenV,t_eigenM=t_eigenM,t_design=t_design,r_eigenM=r_eigenM,r_design=r_design)
  #Careful: mv, diagSW and diagInfo are in descending order for DISPLAY only: they should no be used for scaling designs!
  return(t_des)
}

y_turnDesign<-function(design,y,scale="UV",prec=0.2,preDet=0){
  #print(preDet)
  if(preDet[1]==0) diagOrdHead<-NULL else diagOrdHead<-preDet
  
  if(scale=="UV") {
    res_des<-ctr_des<-uvScaleDesign(design)$MDuv 
    res_y<-uvScaleDesign(y)$MDuv
    SW<-uvScaleDesign(design)$SW 
    diagOrdTail<-sort(diag(SW),decreasing=TRUE,index.return=TRUE)$ix 
  } else if(scale=="OR") {
    res_des<-ctr_des<-orScaleDesign(design)$MDor 
    res_y<-orScaleDesign(y)$MDor
    SW<-orScaleDesign(design)$SW 
    diagOrdTail<-sort(diag(SW),decreasing=TRUE,index.return=TRUE)$ix 
  } else {
    res_des<-ctr_des<-apply(design,2,function(x) x-mean(x))
    res_y<-y-mean(y)
    SW<-diag(1,ncol(design))
    diagOrdTail<-sort(diag(t(ctr_des)%*%ctr_des),decreasing=TRUE,index.return=TRUE)$ix
  }
  diagOrdTail<-setdiff(diagOrdTail,preDet)
  diagOrd<-c(diagOrdHead,diagOrdTail)
  if(preDet[1]!=0) res_des<-res_des[,-preDet]
  
  #rules for using diagOrd: It must only be used for immediate display.
  #                         must not be used to reorder t_eigenM
  #                         must not be used to reorder x_design or t_design or r_design  
  diagInfo<-diag(t(ctr_des)%*%ctr_des)[diagOrd]
  diagSW<-diag(SW)[diagOrd]
  mv<-apply(design,2,mean)[diagOrd]
  #PLS-turning starts here (for just one y no iteration is necessary)
  
  w<-p<-matrix(0,nrow=ncol(res_des),ncol=ncol(res_des))
  rownames(w)<-rownames(p)<-colnames(res_des)
  t<-matrix(0,nrow=nrow(res_des),ncol=ncol(res_des))
  c<-matrix(0,nrow=ncol(res_des),ncol=1)
  for(i in 1:(ncol(res_des))){
    w[,i]<-(t(res_des)%*%res_y)/sqrt(sum((t(res_des)%*%res_y)^2))
    t[,i]<-res_des%*%sCfM(w,i)
    p[,i]<-solve(t(sCfM(t,i))%*%sCfM(t,i))%*%t(sCfM(t,i))%*%res_des
    res_des<-res_des-sCfM(t,i)%*%t(sCfM(p,i))
    c[i]<-solve(t(sCfM(t,i))%*%sCfM(t,i))%*%t(sCfM(t,i))%*%res_y
    res_y<-res_y-sCfM(t,i)*c[i]
  }
  if(preDet[1]==0){
    t_eigenM<-w%*%solve(t(p)%*%w)
  } else {						#add diag(1)-matrix as Eigen-M and diagInfo as Eigen-V
    t_eigenM<-matrix(0,ncol=ncol(ctr_des),nrow=ncol(ctr_des))
    t_eigenM[-preDet,]<-cbind(matrix(0,nrow=ncol(res_des),ncol=length(preDet)),w%*%solve(t(p)%*%w))
    t_eigenM[preDet,]<-cbind(diag(1,length(preDet)),matrix(0,ncol=ncol(res_des),nrow=length(preDet)))
  }
  
  rownames(t_eigenM)<-colnames(design)
  colnames(t_eigenM)<-paste("Comp",1:ncol(t_eigenM),sep="_")
  for(cn in 1:ncol(t_eigenM)){
    rn<-which(abs(t_eigenM[,cn])==max(abs(t_eigenM[,cn])))
    if(t_eigenM[rn,cn]<0) t_eigenM[,cn]<- -t_eigenM[,cn]
    colnames(t_eigenM)[cn] <- paste(rownames(t_eigenM)[rn],cn,sep="_")
  }
  r_eigenM<-myRound(t_eigenM,prec)
  #colnames(t)<-colnames(t_eigenM)
  t_eigenV<-apply(t,2,function(x) sqrt(sum(x^2)))
  t_design<-ctr_des%*%t_eigenM
  r_design<-ctr_des%*%r_eigenM
  t_des<-list(diagOrd=diagOrd,mv=mv,diagSW=diagSW,diagInfo=diagInfo,t_eigenV=t_eigenV,t_eigenM=t_eigenM,t_design=t_design,r_eigenM=r_eigenM,r_design=r_design)
  #Careful: mv, diagSW and diagInfo are in descending order for DISPLAY only: they should no be used for scaling designs!
  return(t_des)
}

evalDesign<-function(design,leaveOut=0,compCut=0,prec=0.2,showPlots=TRUE,preDet=0){
  if(leaveOut[1]) design<-sRfM(design,-leaveOut)
  
  t_des<-turnDesign(design,prec=prec,preDet)
  
  x_mv<-t_des$mv
  diagSW<-t_des$diagSW
  diagOrd<-t_des$diagOrd
  diagInfo<-round(t_des$diagInfo,2)
  t_design<-t_des$t_design
  t_eigenV<-round(t_des$t_eigenV,2)
  t_eigenM<-round(t_des$t_eigenM,2)
  oblongness<-prod(t_eigenV[1:(ncol(design)-compCut)])/prod(diagInfo[1:(ncol(design)-compCut)])
  
  c_design<-sCfM(design[,diagOrd],1:(ncol(design)-compCut))
  t_design<-sCfM(t_des$t_design,1:(ncol(design)-compCut))
  r_design<-sCfM(t_des$r_design,1:(ncol(design)-compCut))
  if(showPlots){
    dev.new()
    pairs(c_design,lower.panel=OKPlot,upper.panel=OKPlot,main="Correlations of dl Design")  
    dev.new()
    pairs(t_design,lower.panel=OKPlot,upper.panel=OKPlot,main="Correlations of \"turned\" dl Design")
    dev.new()
    pairs(r_design,lower.panel=OKPlot,upper.panel=OKPlot,main="Correlations of \"rounded&turned\" dl Design")
  }
  r_eigenM<-t_des$r_eigenM
  t_design<-t_des$t_design
  r_design<-t_des$r_design
  
  evaluation<-list(diagInfo=diagInfo,diagSW=diagSW,t_eigenV=t_eigenV,oblongness=oblongness,t_eigenM=t_eigenM,r_eigenM=r_eigenM,t_design=t_design,r_design=r_design)
  return(evaluation)
}

evalDesignAndY<-function(design,y,leaveOut=0,compCut=0,preDet=0,prec=0.2,showPlots=TRUE){
  print(preDet)
  if(leaveOut[1]) design<-sRfM(design,-leaveOut)
  if(leaveOut[1]) y<-sRfM(y,-leaveOut)  
  t_des<-y_turnDesign(design,y,preDet=preDet,prec=prec)
  
  x_mv<-t_des$mv
  diagSW<-t_des$diagSW
  diagOrd<-t_des$diagOrd
  diagInfo<-round(t_des$diagInfo,2)
  t_design<-t_des$t_design
  t_eigenV<-round(t_des$t_eigenV,2)
  t_eigenM<-round(t_des$t_eigenM,2)
  oblongness<-prod(t_eigenV[1:(ncol(design)-compCut)])/prod(diagInfo[1:(ncol(design)-compCut)])
  
  c_design<-sCfM(design[,diagOrd],1:(ncol(design)-compCut))
  t_design<-sCfM(t_des$t_design,1:(ncol(design)-compCut))
  r_design<-sCfM(t_des$r_design,1:(ncol(design)-compCut))
  if(showPlots){
    dev.new()
    pairs(c_design,lower.panel=OKPlot,upper.panel=OKPlot,main="Correlations of dl Design")  
    dev.new()
    pairs(t_design,lower.panel=OKPlot,upper.panel=OKPlot,main="Correlations of \"turned\" dl Design")
    dev.new()
    pairs(r_design,lower.panel=OKPlot,upper.panel=OKPlot,main="Correlations of \"rounded&turned\" dl Design")
  }
  r_eigenM<-t_des$r_eigenM
  t_design<-t_des$t_design
  r_design<-t_des$r_design
  
  evaluation<-list(diagInfo=diagInfo,diagSW=diagSW,t_eigenV=t_eigenV,oblongness=oblongness,t_eigenM=t_eigenM,r_eigenM=r_eigenM,t_design=t_design,r_design=r_design,y=y)
  return(evaluation)
}


#internal design function for 2 level full fac
fullFactorial<-function(numfac){
  vorz_sp<-matrix(nrow=2^numfac,ncol=numfac)
  for(d in 1:numfac){
    vorz_sp[,d]<-t(t(rep(c(rep(-1,2^(d-1)),rep(1,2^(d-1))),2^(numfac-d))))
  }
  return(vorz_sp)
}
#internal design function for 3 level full fac
full3LevelFactorial<-function(numfac){
  vorz_sp<-matrix(nrow=3^numfac,ncol=numfac)
  for(d in 1:numfac){
    vorz_sp[,d]<-t(t(rep(c(rep(-1,3^(d-1)),rep(0,3^(d-1)),rep(1,3^(d-1))),3^(numfac-d))))
  }
  return(vorz_sp)
}
#internal design function for 4 level full fac - used to generate candidates
full4LevelFactorial<-function(numfac){
  vorz_sp<-matrix(nrow=4^numfac,ncol=numfac)
  for(d in 1:numfac){
    vorz_sp[,d]<-t(t(rep(c(rep(-1,4^(d-1)),rep(-1/3,4^(d-1)),rep(1/3,4^(d-1)),rep(1,4^(d-1))),4^(numfac-d))))
  }
  return(vorz_sp)
}

#internal design function for 2 level fractionalfac
reduceFactorial<-function(halvedFull,redcol,generator,vorz){
  halvedFull[,redcol]<-apply(halvedFull[,generator],1,function(x) vorz*prod(x))
  return(halvedFull)
}

#internal design function for adding centre points to a design
addCentrePoints<-function(design,numberCP){
  design<-rbind(design,rbind(matrix(0,ncol=ncol(design),nrow=numberCP)))
  return(design)
}

#internal design function for adding starpoints to a design
addStarPoints<-function(design,whFsq,starDistance){
  design<-rbind(design,rbind(matrix(0,ncol=ncol(design),nrow=2*length(whFsq))))
  for(i in (0:(2*length(whFsq)-1)) ) design[nrow(design)-i,whFsq[length(whFsq)-floor((i)/2)]]<-starDistance*(-1)^i
  return(design)
}

#internal quasi-design function to create design-matrix part for const and scup
constDepScupDesign<-function(u_low,u_high,nR,in_constDep,in_scup,useLowForScup=TRUE)
{
  if(ncol(in_constDep)) {
    constDepDes<-matrix(rep(t(u_low)%*%in_constDep,times=nR),ncol=ncol(in_constDep),byrow=TRUE)
    colnames(constDepDes)<-colnames(in_constDep)
  } else constDepDes<-matrix(nrow=nR,ncol=0)
  if(ncol(in_scup)) {  
    if(useLowForScup)
      scupWishes<-matrix(rep(t(u_low)%*%in_scup,times=nR),ncol=ncol(in_scup),byrow=TRUE) else
        scupWishes<-matrix(rep(t(u_high)%*%in_scup,times=nR),ncol=ncol(in_scup),byrow=TRUE)
      colnames(scupWishes)<-colnames(in_scup)
  } else scupWishes<-matrix(nrow=nR,ncol=0)
  constDepScupDes<-cbind(constDepDes,scupWishes)
  return(constDepScupDes)
}

#function to project candidates in a candidate set before calling d-optimal in order to fulfill linear constraints 
projectCandidates<-function(x_candidates,coeffs,lowLims,highLims,eps=0.00001,debug=0)
  #probably deprecated
  #how to feed projectCandidates with constraint equations:
  #coeffs<-SW%*%Wtop                              #ndes rows and nctr cols
  #lowLims<-U_D_cp[,new_roles=="contr"]-          #nruns cols and nctr cols
  #      matrix(rep(t(new_high)[,new_roles=="contr"],nrow(x_candidates)),nrow=nrow(x_candidates),byrow=TRUE)
  #lowLims<-U_D_cp[,new_roles=="contr"]-          #nruns cols and nctr cols
  #      matrix(rep(t(new_low)[,new_roles=="contr"],nrow(x_candidates)),nrow=nrow(x_candidates),byrow=TRUE)
{ 
  WtWt<-diag((t(coeffs)%*%coeffs))
  for( i in 1:length(WtWt))
  {
    tooHigh<-x_candidates%*%coeffs+highLims	#if positive outside of domain
    tooLow<-x_candidates%*%coeffs+lowLims       #if negative outside of domain
    
    x_candidates[sCfM(tooHigh,i)>eps,]<-x_candidates[sCfM(tooHigh,i)>eps,]-
      ginv(WtWt[i])[1,1]*cbind(sCfM(tooHigh,i)[sCfM(tooHigh,i)>eps,])%*%(coeffs)[,i]
    
    x_candidates[sCfM(tooLow,i)<(-eps),]<-x_candidates[sCfM(tooLow,i)<(-eps),]-
      ginv(WtWt[i])[1,1]*cbind(sCfM(tooLow,i)[sCfM(tooLow,i)<(-eps),])%*%(coeffs)[,i]
  }
  
  if(debug)
  {
    tooHigh<-x_candidates%*%coeffs+highLims
    print(tooHigh)
    tooLow<-x_candidates%*%coeffs+lowLims
    print(tooLow)
  }
  return(x_candidates)
}

#function to remove candidates from a candidate set before calling d-optimal in order to fulfill linear constraints 
removeCandidates<-function(x_candidates,coeffs,lowLims,highLims,eps=0.00001,debug=0)
  #probably also deprecated
  #how to feed removeCandidates with constraint equations:
  #coeffs<-SW%*%Wtop                              #ndes rows and nctr cols
  #lowLims<-U_D_cp[,new_roles=="contr"]-          #nruns cols and nctr cols
  #      matrix(rep(t(new_high)[,new_roles=="contr"],nrow(x_candidates)),nrow=nrow(x_candidates),byrow=TRUE)
  #lowLims<-U_D_cp[,new_roles=="contr"]-          #nruns cols and nctr cols
  #      matrix(rep(t(new_low)[,new_roles=="contr"],nrow(x_candidates)),nrow=nrow(x_candidates),byrow=TRUE)
{ 
  WtWt<-diag((t(coeffs)%*%coeffs))
  for( i in 1:length(WtWt))
  {
    tooHigh<-x_candidates%*%coeffs+highLims
    tooLow<-x_candidates%*%coeffs+lowLims
    
    x_candidates<-x_candidates[sCfM(tooHigh,i)<eps,]
    x_candidates<-x_candidates[sCfM(tooLow,i)>(-eps),]
  }
  
  if(debug)
  {
    tooHigh<-x_candidates%*%coeffs+highLims
    print(tooHigh)
    tooLow<-x_candidates%*%coeffs+lowLims
    print(tooLow)
  }
  return(x_candidates)
}


#function to extend a row for an observation or a whole design matrix
extend_x<-function(x,expo)
  #expo is a multi-index of the form "list(0,1,2)" (easy linear model for 2 factors) 
  #  "list(0,1,2,3,4,5,c(1,1),c(2,2),c(3,3),c(4,4),c(5,5),c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5))"
  #  (not quite so easy cross&square-model for 5 factors)
{
  x_ext<-matrix(1,nrow=nrow(x),ncol<-length(expo))
  #print(x_ext)
  for(i in 1:length(expo)) 
  {	
    for(j in expo[[i]]) 
    {
      #print(j)
      if(j!=0){x_ext[,i]<-x_ext[,i]*x[,j]}
    }
  }
  colnames(x_ext)<-expo
  return(x_ext)
}

#function to extend a u_design row for an observation or a whole u_design matrix
extend_u<-function(u,V,expo)		
  #expo is a multi-index of the form "list(0,1,2)" (easy linear model for 2 factors) 
  #  "list(0,1,2,3,4,5,c(1,1),c(2,2),c(3,3),c(4,4),c(5,5),c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5))"
  #  (not quite so easy cross&square-model for 5 factors)
{
  x_ext<-matrix(1,nrow=nrow(u),ncol<-length(expo))
  
  for(i in 1:length(expo)) 
  {	
    for(j in expo[[i]]) 
    {
      #print(j)
      if(j!=0){x_ext[,i]<-x_ext[,i]*u%*%V[,j]}
    }
  }
  colnames(x_ext)<-expo
  return(x_ext)
}

#function to create a list that represents a linear model e.g. list(0,1,2)
linearDimDoeMod<-function(dlInd,dleqInd,test_eqpr=TRUE,const=TRUE){
  if(const) dimDoeMod<-list(0) else
    dimDoeMod<-list()
  if(test_eqpr)
    dimDoeMod<-append(dimDoeMod,dleqInd) else
      dimDoeMod<-append(dimDoeMod,dlInd)
    return(dimDoeMod)
}

#function to create a list that represents a interaction model e.g. list(0,1,2,c(1,2))
interactionDimDoeMod<-function(dlInd,dleqInd,test_eqpr=TRUE,const=TRUE){
  if(const) dimDoeMod<-list(0) else
    dimDoeMod<-list()
  if(test_eqpr)
    dimDoeMod<-append(dimDoeMod,dleqInd) else
      dimDoeMod<-append(dimDoeMod,dlInd)
    if(length(dlInd)>1) for(i in 1:(length(dlInd)-1))
      for(j in (i+1):length(dlInd)) {
        term<-c(i,j)
        dimDoeMod[[length(dimDoeMod)+1]]<-term
      }
    return(dimDoeMod)
}

#function to create a list that represents a quadratic model e.g. list(0,1,2,c(1,2),c(1,1),c(2,2))
quadraticDimDoeMod<-function(dlInd,dleqInd,test_eqpr=TRUE,const=TRUE){
  if(const) dimDoeMod<-list(0) else
    dimDoeMod<-list()
  if(test_eqpr)
    dimDoeMod<-append(dimDoeMod,dleqInd) else
      dimDoeMod<-append(dimDoeMod,dlInd)
    if(length(dlInd)>0) for(i in 1:length(dlInd)) {
      term<-c(i,i)
      dimDoeMod[[length(dimDoeMod)+1]]<-term
      if(i<length(dlInd))
        for(j in (i+1):length(dlInd)) {
          term<-c(i,j)
          dimDoeMod[[length(dimDoeMod)+1]]<-term
        }
    }
    return(dimDoeMod)
}

#function to create a list that represents a hybrid model e.g. list(0,1,2,c(1,2),c(1,1))
hybridDimDoeMod<-function(hl,dlInd,dleqInd,test_eqpr=TRUE,const=TRUE){
  if(const) dimDoeMod<-list(0) else
    dimDoeMod<-list()
  if(test_eqpr)
    dimDoeMod<-append(dimDoeMod,dleqInd) else
      dimDoeMod<-append(dimDoeMod,dlInd)
    for(i in 1:length(hl)) {
      if(hl[i]>2) {
        term<-c(i,i)
        dimDoeMod[[length(dimDoeMod)+1]]<-term
      }
      if(i<length(dlInd))
        for(j in (i+1):length(dlInd)) {
          term<-c(i,j)
          dimDoeMod[[length(dimDoeMod)+1]]<-term
        }
    }
    return(dimDoeMod)
}

#####################################################
#              Variation 1: factorial in RF-space
#####################################################

createDesigninRFspace<-function(u_low,u_high,in_contr,in_eqpr,dlInd,dleqInd,desType,desParam,debug=0)
  #returns just a u_design with columns for contr and eqpr factors
  #desType contains the function name of the design to be generated - can be user defined!
{
  flag<-"ok"
  desParam$numfac<-ncol(in_contr)+ncol(in_eqpr)
  u_design<-0
  model<-0
  if(substr(try(substr(cccc<-try(match.fun(desType),silent=TRUE)[1],1,1),silent=TRUE)[1],1,5)!="Error") flag<-"design Type not found" else
  {
    modde<-match.fun(desType)(desParam,dlInd,dleqInd)
    MD<-as.matrix(modde$design)
    model<-modde$model
    u_low<-t(cbind(in_contr,in_eqpr))%*%u_low
    u_high<-t(cbind(in_contr,in_eqpr))%*%u_high
    u_design<-matrix(ccW(u_low,u_high,0)$cp,nrow=nrow(MD),ncol=desParam$numfac,byrow=TRUE)+MD%*%ccW(u_low,u_high,0)$SW
    if(debug) print(head(10^u_design))
  }
  
  modde<-list(design=u_design,model=model,flag=flag)
  return(modde)
}


#####################################################
#      Variation 2: importing a design in  RF-space (original coordinates)
#####################################################

importDesigninRFspace<-function(u_design,in_contr,in_eqpr,dlInd,dleqInd,debug=0)
{
  flag<-"ok"
  if(ncol(u_design)!=ncol(in_contr)+ncol(in_eqpr)) flag<-"Inconsistent col number in imported RF-design" else
    if(ncol(u_design)!=length(intersect(colnames(u_design),union(colnames(in_contr),colnames(in_eqpr))))) flag<-"Inconsistent column names in imported RF-design"
    #u_design<-log10(u_design)
    if(debug) print(head(u_design))
    model<-linearDimDoeMod(dlInd,dleqInd)		#i.e. just linear
    modde<-list(design=u_design,model=model,flag=flag)
    return(modde)
}


#####################################################################################
#      Variation 3: importing a scaled and centred design in  DLF-space (-1 +1 -array)
#      value: x_design for dldes and eqpr-factors - no cost or scup
#             dldes-design scaled to x_low/high
#####################################################################################


importSCDesigninDLFspace<-function(MD,x_low,x_high,in_eqpr,dlInd,dleqInd,test_eqpr,debug=0)
{
  flag<-"ok"
  x_D_cp<-ccW(x_low,x_high,0)$cp[dleqInd]			#scaling for both DLF and in_eqpr
  SW<-ccW(x_low,x_high,0)$SW[dleqInd,dleqInd]
  
  #scaled design in x-space: MD
  
  if(test_eqpr) {
    if(ncol(MD)!=length(dleqInd)) {
      flag<-"Inconsistent col number in imported sc DL-design, should include eqpr-factors"
    }
  } else {
    if(ncol(MD)!=length(dlInd)) {
      flag<-"Inconsistent col number in imported sc DL-design, should only include DL-factors"
    }
    eqprDes<-matrix(0,ncol=ncol(in_eqpr),nrow=nrow(MD))
    colnames(eqprDes)<-colnames(in_eqpr)
    MD<-cbind(MD,eqprDes)
    SW[colnames(in_eqpr),colnames(in_eqpr)]<-0			#scaling for eqpr set to 0 ?????????????????
  }
  if(debug) print(MD)
  
  if(flag=="ok"){
    X_D_cp<-t(matrix(rep(x_D_cp,times=nrow(MD)),nrow=length(x_D_cp)))
    colnames(X_D_cp)<-colnames(x_D_cp)
    if(debug)print(head(X_D_cp))
    x_design<-X_D_cp+MD%*%SW						
    model<-linearDimDoeMod(dlInd,dleqInd,test_eqpr)
    modde<-list(design=x_design,model=model,flag=flag)
  } else modde<-list(flag=flag)
  return(modde)
}



####################################################################
#      Variation 4: Design in  DLF-space
#      value: x_design for dldes and eqpr-factors - no const or scup
#             
####################################################################

createDesigninDLFspace<-function(x_low,x_high,in_eqpr,dlInd,dleqInd,test_eqpr,desType,desParam,debug=0)
{
  flag<-"ok"
  x_design<-0
  model<-0
  x_D_cp<-ccW(x_low,x_high,0)$cp[dleqInd]			#scaling for both DLF and in_eqpr
  SW<-ccW(x_low,x_high,0)$SW[dleqInd,dleqInd]
  
  if(test_eqpr){
    desParam$numfac<-length(dleqInd)
  } else {
    desParam$numfac<-length(dlInd)
  }
  
  #scaled design in x-space
  if(substr(try(substr(cccc<-try(match.fun(desType),silent=TRUE)[1],1,1),silent=TRUE)[1],1,5)!="Error") flag<-"design Type not found" else
  {
    modde<-match.fun(desType)(desParam,dlInd,dleqInd,test_eqpr)
    MD<-as.matrix(modde$design)
    model<-modde$model
    if(!test_eqpr){					#needs to be augmented
      eqprDes<-matrix(0,ncol=ncol(in_eqpr),nrow=nrow(MD))
      colnames(eqprDes)<-colnames(in_eqpr)
      MD<-cbind(MD,eqprDes)
      SW[colnames(in_eqpr),colnames(in_eqpr)]<-0			#scaling for eqpr set to 0
    }
    
    if(debug) print(MD)
    X_D_cp<-matrix(rep(x_D_cp,times=nrow(MD)),ncol=length(x_D_cp),byrow=TRUE)
    colnames(X_D_cp)<-colnames(x_D_cp)
    if(debug)print(head(X_D_cp))
    if(prod(diag(SW))==0) flag="warning: one of the factors has constant low/high setting"
    x_design<-X_D_cp+MD%*%SW						
  }
  modde<-list(design=x_design,model=model,flag=flag)
  return(modde)
}



#########################################################################################
#      Variation 5: D-optimal design in  DLF-space to satisfy low/high levels in RT-space
#      value: x_design for dldes and eqpr-factors - no cost or scup
#             dldes-design scaled to x_low/high
#########################################################################################


createUconstrainedDoptimalinDLFspace<-function(u_candidates,VRES,WRES,u_low,u_high,x_low,x_high,x_roles,in_constDep,in_scup,in_eqpr,dlInd,dleqInd,test_eqpr,desType,desParam,debug=0,eps=0.00001)
{
  
  flag<-"ok"
  #x_D_cp<-ccW(x_low,x_high,0)$cp[x_low!=x_high]			#scaling for both DLF and in_eqpr
  #SW<-ccW(x_low,x_high,0)$SW[x_low!=x_high,x_low!=x_high]
  
  if(test_eqpr){
    desParam$numfac<-length(dleqInd)      
  } else {
    desParam$numfac<-length(dlInd)
  }
  
  if(desParam$import){
    if(test_eqpr) {
      if(ncol(u_candidates)!=length(dleqInd)) flag<-"Inconsistent col number in imported u_candidates, should include eqpr-factors"
    } else {
      if(ncol(u_candidates)!=length(dlInd)) flag<-"Inconsistent col number in imported u_candidates, should only include DL-factors"
      eqprCtr<-sCfM(t(10^((u_low+u_high)/2)),colnames(in_eqpr))
      eqprDes<-matrix(eqprCtr,ncol=ncol(in_eqpr),nrow=nrow(u_candidates),byrow=TRUE)         #1 because logs taken later
      colnames(eqprDes)<-colnames(in_eqpr)
      u_candidates<-cbind(u_candidates,eqprDes)
    }
    if(flag=="ok"){
      u_candidates<-log10(u_candidates)
      
      constDepScupDes<-constDepScupDesign(u_low,u_high,nrow(u_candidates),in_constDep,in_scup,useLowForScup=TRUE)
      u_candidates<-cbind(u_candidates,constDepScupDes)[,rownames(VRES)]
      u_candidates<-u_candidates[ncol(u_candidates)==apply(u_candidates<=eps+matrix(rep(t(u_high),times=nrow(u_candidates)),ncol=ncol(u_candidates),byrow=TRUE),1,sum),]
      u_candidates<-u_candidates[ncol(u_candidates)==apply(u_candidates+eps>=matrix(rep(t(u_low),times=nrow(u_candidates)),ncol=ncol(u_candidates),byrow=TRUE),1,sum),]
      x_cp<-ccW(x_low,x_high,0)$cp
      SW<-ccW(x_low,x_high,0)$SW
    } #else error on import, but flag is still set
  } else { #(!desParam$import)
    if(length(desParam$hyblev))
    {
      MD_candidates<-designFullHybLevFactorial(desParam,dlInd,dleqInd,test_eqpr)$design 
      #u_low<-t(cbind(in_contr,in_eqpr))%*%u_low
      #u_high<-t(cbind(in_contr,in_eqpr))%*%u_high
      u_candidates<-matrix(ccW(u_low,u_high,0)$cp,nrow=nrow(MD_candidates),ncol=nrow(u_low),byrow=TRUE)+MD_candidates%*%ccW(u_low,u_high,0)$SW
      x_cp<-ccW(x_low,x_high,0)$cp
      SW<-ccW(x_low,x_high,0)$SW
      
    } else {
      MD_candidates<-full4LevelFactorial(desParam$numfac)
      MD_candidates<-rbind(MD_candidates,full3LevelFactorial(desParam$numfac))
      MD_candidates<-addCentrePoints(MD_candidates,1)
      MD_candidates<-addStarPoints(MD_candidates,(1:ncol(MD_candidates)),1)
      
      if(!test_eqpr){
        MD_candidates<-cbind(MD_candidates,matrix(0,ncol=ncol(in_eqpr),nrow=nrow(MD_candidates)))
      } 
      colnames(MD_candidates)<-colnames(VRES)[dleqInd]
      constDepScupDes<-constDepScupDesign(u_low,u_high,nrow(MD_candidates),in_constDep,in_scup,useLowForScup=TRUE)
      
      MD_candidates<-cbind(MD_candidates,constDepScupDes)[,colnames(VRES)]
      
      x_cp<-ccW(x_low,x_high,0)$cp
      SW<-ccW(x_low,x_high,0)$SW
      X_cp<-matrix(rep(x_cp,times=nrow(MD_candidates)),ncol=length(x_cp),byrow=TRUE)
      x_candidates<-X_cp+MD_candidates%*%SW
      u_candidates<-x_candidates%*%WRES
    }
    if(length(desParam$projectCand)!=0)
    { #project
      u_candidates<-pmin(u_candidates,matrix(rep(t(u_high),times=nrow(u_candidates)),ncol=ncol(u_candidates),byrow=TRUE))
      u_candidates<-pmax(u_candidates,matrix(rep(t(u_low),times=nrow(u_candidates)),ncol=ncol(u_candidates),byrow=TRUE))
    } else {
      #cut off
      #print(head(u_candidates))
      u_candidates<-u_candidates[ncol(u_candidates)==apply(u_candidates<=eps+matrix(rep(t(u_high),times=nrow(u_candidates)),ncol=ncol(u_candidates),byrow=TRUE),1,sum),]
      #print(u_candidates)
      #print(u_high)
      u_candidates<-u_candidates[ncol(u_candidates)==apply(u_candidates+eps>=matrix(rep(t(u_low),times=nrow(u_candidates)),ncol=ncol(u_candidates),byrow=TRUE),1,sum),]
    }
  }
  if(flag=="ok"){
    x_candidates<-u_candidates%*%VRES
    x_candidates<-x_candidates[ncol(x_candidates)==apply(x_candidates<=eps+matrix(rep(t(x_high),times=nrow(x_candidates)),ncol=ncol(x_candidates),byrow=TRUE),1,sum),]
    x_candidates<-x_candidates[ncol(x_candidates)==apply(x_candidates+eps>=matrix(rep(t(x_low),times=nrow(x_candidates)),ncol=ncol(x_candidates),byrow=TRUE),1,sum),]
    
    MD_candidates<-scaleDesign(x_candidates,x_cp,SW)
    if(test_eqpr){
      MD_desCandidates<-sCfM(MD_candidates,dleqInd)
    } else {
      MD_desCandidates<-sCfM(MD_candidates,dlInd)
    }
    dimDoeMod<-desParam$dimDoeMod
    if(length(dimDoeMod)==0) dimDoeMod<-linearDimDoeMod(dlInd,dleqInd,test_eqpr)
    DmodelTemplate<-transDimDoeMod(dimDoeMod,x_roles)
    MD<-MD_candidates[optFederov(DmodelTemplate,data.frame(MD_desCandidates))$rows,]
    MD<-addCentrePoints(MD,desParam$numberCP)
    X_D_cp<-matrix(rep(x_cp,times=nrow(MD)),ncol=length(x_cp),byrow=TRUE)
    colnames(X_D_cp)<-colnames(x_cp)
    if(debug)print(head(X_D_cp))
    x_design<-X_D_cp+MD%*%SW						#cols for in_eqpr should be centred and scaled
  } else x_design<-0
  modde<-list(design=x_design,dimDoeMod=dimDoeMod,flag=flag)
  return(modde)
  #print(modde)
}

#function to calculate y values (logs) from z in SI-units
fetchDL_resp<-function(u_matrix,z_val,nRF,V_resp)
  #z_val ist jetzt berets logarithmiert - muss ggfs. zuvor in uu2u passiert
{
  flag="ok"
  #l_zv<-log10(z_val)
  if(nrow(z_val)!=nrow(u_matrix)) {
    flag<-"number of response values not equal to the number of observations" 
    yv<-0
  } else if(ncol(z_val)!=nrow(V_resp)-nRF) {
    flag<-"number of response values not equal to the number of observations" 
    yv<-0
  } else  {
    yv<-u_matrix%*%sRfM(V_resp,1:nRF)+z_val%*%sRfM(V_resp,(nRF+1):nrow(V_resp))
    colnames(yv)<-colnames(V_resp)
  }
  DL_resp<-list(flag=flag,yv=yv)
  return(DL_resp)
}



###################################################################
###################################################################
#     external design functions   (wrappers useful for both x- and u-designs)
#     user can modify or add design functions
#     rule 1: name of function must start with design
#     rule 2: input parameter must be a named list called desParam
#     rule 3: return-object must be a list, containing model AND sc-design, as well as an ok-flag
###################################################################
###################################################################

#wrapper for a 2-level full factorial design
designFullFactorial<-function(desParam,dlInd,dleqInd,test_eqpr=TRUE){
  #  	desParam$numfac: number of factors to be varied
  #  	desParam$numberCP: number of centre points to be added)
  design<-fullFactorial(desParam$numfac)
  if(desParam$numberCP) design<-addCentrePoints(design,desParam$numberCP)
  model<-interactionDimDoeMod(dlInd,dleqInd,test_eqpr)
  modde<-list(model=model,design=design)
  return(modde)
}

#wrapper for a N-level full factorial design
designFullNLevelFactorial<-function(desParam,dlInd,dleqInd,test_eqpr=TRUE){
  #  	desParam$numfac: number of factors to be varied
  #     desParam$numlev: number of levels
  #  	desParam$numberCP: number of centre points to be added)
  if(desParam$numlev==2)design<-fullFactorial(desParam$numfac) else
    if(desParam$numlev==3)design<-full3LevelFactorial(desParam$numfac) else
      design<-full4LevelFactorial(desParam$numfac)
    if(desParam$numberCP) design<-addCentrePoints(design,desParam$numberCP)
    model<-quadraticDimDoeMod(dlInd,dleqInd,test_eqpr)
    modde<-list(model=model,design=design)
    return(modde)
}


#wrapper for a hybrid-level full factorial design
designFullHybLevFactorial<-function(desParam,dlInd,dleqInd,test_eqpr=TRUE){
  #  	desParam$numfac: number of factors to be varied
  #     desParam$hyblev: vector with number of levels
  #  	desParam$numberCP: number of centre points to be added)
  hl<-desParam$hyblev
  vorz_sp<-matrix(0,nrow=prod(hl[1:length(hl)]),ncol=length(hl))
  
  if(length(hl)>1)
    #d=1
    for(j in (1:prod(hl[2:length(hl)]))){
      if(hl[1]>1) for(i in 1:hl[1]){
        vorz_sp[(1+(j-1)*hl[1]+(i-1)):((j-1)*hl[1]+i),1]<- -1+2*(i-1)/(hl[1]-1)
      }
    }
  
  if(length(hl)>2)
    for(d in 2:(length(hl)-1)){
      if(hl[d]>1) for(j in (1:prod(hl[(d+1):length(hl)]))){
        for(i in 1:hl[d]){
          vorz_sp[(1+(j-1)*prod(hl[1:d])+(i-1)*prod(hl[1:(d-1)])):((j-1)*prod(hl[1:d])+i*prod(hl[1:(d-1)])),d]<-rep(-1+2*(i-1)/(hl[d]-1),prod(hl[1:(d-1)]))
        }
      }
    }
  #d=length(hl), j=1
  if(hl[length(hl)]>1) for(i in 1:hl[length(hl)]){
    vorz_sp[(1+(i-1)*prod(hl[1:(length(hl)-1)])):(i*prod(hl[1:(length(hl)-1)])),length(hl)]<-rep(-1+2*(i-1)/(hl[length(hl)]-1),prod(hl[1:(length(hl)-1)]))
  }
  
  
  design<-vorz_sp
  if(desParam$numberCP) design<-addCentrePoints(design,desParam$numberCP)
  model<-hybridDimDoeMod(hl,dlInd,dleqInd,test_eqpr)
  modde<-list(model=model,design=design)
  return(modde)
}


#wrapper for fractional factorial design
designFractionalFactorial<-function(desParam,dlInd,dleqInd,test_eqpr=TRUE){
  dim<-desParam$numfac				#number of factors, n
  vorz<-desParam$vorz				#list of signs with length reduce 
  reduce<-desParam$reduce			#number of generators, k = 1 or 2
  addCentre<-desParam$numberCP		#number of centre points
  
  design<-fullFactorial(dim)
  
  if(length(vorz)!=reduce) vorz[1:reduce]<-1        #set default
  if(reduce==1&(dim>=3)) design<-reduceFactorial(design[1:(nrow(design)/2),], dim, 1:(dim-1),vorz[1])
  if(reduce==2&(dim>=5)){
    design<-reduceFactorial(design[1:(nrow(design)/2),], dim, 1:3,vorz[1])
    design<-reduceFactorial(design[1:(nrow(design)/2),], dim-1, 2:(dim-2),vorz[2])
  }
  if(addCentre)design<-addCentrePoints(design,addCentre)
  if(dim-reduce>3)
    model<-interactionDimDoeMod(dlInd,dleqInd,test_eqpr) else
      model<-linearDimDoeMod(dlInd,dleqInd,test_eqpr)
  modde<-list(model=model,design=design)
  return(modde)
}

#wrapper for a central composite design
designCCD<-function(desParam,dlInd,dleqInd,test_eqpr=TRUE){
  dim<-desParam$numfac				#number of factors, n
  whichFactorssquare<-desParam$whichFactorssquare
  starDistance<-desParam$starDistance
  addCentre<-desParam$numberCP
  
  design<-fullFactorial(dim)
  design<-addStarPoints(design,whichFactorssquare,starDistance)
  if(addCentre)design<-addCentrePoints(design,addCentre)
  
  model<-quadraticDimDoeMod(dlInd,dleqInd,test_eqpr)
  modde<-list(model=model,design=design)
  return(modde)
}


#wrapper for d-optimal design using Federov's exchange algorithm in AlgDesign package
designDoptimal<-function(desParam,dlInd,dleqInd,test_eqpr=TRUE){
  
  candidates<-designFullNLevelFactorial(desParam,dlInd,dleqInd,test_eqpr)
  if(length(desParam$DoptModel)&&!desParam$DoptModel[[1]]) {
    model<-desParam$DoptModel
    ext_candidates<-extend_x(candidates,model)
  } else {						#if no model is given linear will be used
    model<-linearDimDoeMod(dlInd,dleqInd,test_eqpr)
    ext_candidates<-candidates
  }
  design<-optFederov(~.,ext_candidates)$design[,1:ncol(candidates)]
  modde<-list(model=model,design=design)
  return(modde)
}


##########################################################################
##########################################################################
#     functions that communicate with the user
##########################################################################
##########################################################################

#function to fetch all available data into the list dataList
readDimDoeData<-function(ex_path){
  dataList<-list()
  tryCatch(dataList$U<-(read.csv(paste(ex_path,"U.csv",sep=""), sep=";", dec=",", row.names=1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$qualRel<-as.matrix(read.csv(paste(ex_path,"qualRel.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$relayDesign1<-as.matrix(read.csv(paste(ex_path,"relayDesign1.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$relayDesign2<-as.matrix(read.csv(paste(ex_path,"relayDesign2.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$relayDesign3<-as.matrix(read.csv(paste(ex_path,"relayDesign3.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$relayDesign4<-as.matrix(read.csv(paste(ex_path,"relayDesign4.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$relayCoeffs<-as.matrix(read.csv(paste(ex_path,"relayCoeffs.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$x_levels<-as.matrix(read.csv(paste(ex_path,"x_levels.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$u_design<-as.matrix(read.csv(paste(ex_path,"u_design.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$uu_design<-data.frame(read.csv(paste(ex_path,"uu_design.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$MD<-as.matrix(read.csv(paste(ex_path,"MD.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$u_candidates<-as.matrix(read.csv(paste(ex_path,"u_candidates.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$uu_candidates<-data.frame(read.csv(paste(ex_path,"uu_candidates.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$Z<-read.csv(paste(ex_path,"Z.csv",sep=""), sep=";", dec=",", row.names = 1),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$z_values<-as.matrix(read.csv(paste(ex_path,"z_values.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$uz_values<-dataframe(read.csv(paste(ex_path,"uz_values.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$u_new<-as.matrix(read.csv(paste(ex_path,"u_new.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$uu_new<-data.frame(read.csv(paste(ex_path,"uu_new.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$z_new<-as.matrix(read.csv(paste(ex_path,"z_new.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  tryCatch(dataList$uz_new<-data.frame(read.csv(paste(ex_path,"uz_new.csv",sep=""), sep=";", dec=",", row.names = 1)),error=function(e) e, warning = function(w) w<-w)
  
  return(dataList)
}

#importing factor information from U.csv
fetchInfosFromU<-function(U,debug=0)
{ 
  
  rownames(U)<-U[,2]; head(U)
  #falls hierbei ein Fehler auftritt: 
  # - "kann Verbindung nicht oeffnen": dann Pfad und Name genau ueberpruefen
  # - "doppelte rownames nicht zulaessig": bitte nach der letzten Zeile in Excel die naechsten Zeilen mit "Strg--" entfernen
  #    denn wahrscheinlich gibt es leere Zeilen, mit denen R nichts anfangen kann
  
  #names
  u_names<-as.matrix(U[,1])
  colnames(u_names)<-"names"
  rownames(u_names)<-U[,2]
  if(debug) print(head(u_names))
  
  u_roles<-as.matrix(U[,3])
  colnames(u_roles)<-"roles"
  rownames(u_roles)<-U[,2]
  if(debug) print(head(u_roles))
  
  #userunits
  u_userunit<-as.matrix(U[,4])
  colnames(u_userunit)<-"userunit"
  rownames(u_userunit)<-U[,2]
  if(debug) print(head(u_userunit))
  
  #offsets
  u_offset<-as.matrix(U[,8])
  colnames(u_offset)<-"offset" 
  rownames(u_offset)<-U[,2] 
  if(debug) print(head(u_offset))
  
  #gradients
  u_gradient<-as.matrix(U[,9])
  colnames(u_gradient)<-"gradient" 
  rownames(u_gradient)<-U[,2] 
  if(debug) print(head(u_gradient))
  
  #transforms
  u_transform<-as.matrix(U[,5])
  colnames(u_transform)<-"transform" 
  rownames(u_transform)<-U[,2] 
  if(debug) print(head(u_transform))
  
  #low limits 
  u_low<-as.matrix(U[,6])
  rownames(u_low)<-U[,2] 
  u_low<-u_offset+u_gradient*u_low      #so to invert: uu_low<-(u_low-offset)/gradient
  u_low<-apply(data.frame(u_low,u_transform),1, function(x) {if(x[2]=="none") as.numeric(x[1]) else log10(as.numeric(x[1]))} )
  #vielleicht besser mit for-Schleife und dann
  #if(u_transform[i]=="log") u_low[i,]<-log10(u_low[i,]) dann kann man weitere Transf. zulassen
  #vergl. Funktionen u2uu und uu2u (s.  unten)
  u_low<-as.matrix(data.frame(u_low))
  colnames(u_low)<-"l.low"
  if(debug) print(head(u_low))
  
  #high limits 
  u_high<-as.matrix(U[,7])
  rownames(u_high)<-U[,2] 
  u_high<-u_offset+u_gradient*u_high
  u_high<-apply(data.frame(u_high,u_transform),1, function(x) {if(x[2]=="none") as.numeric(x[1]) else log10(as.numeric(x[1]))} )
  u_high<-as.matrix(data.frame(u_high))
  colnames(u_high)<-"l.high" 
  if(debug) print(head(u_high))
  
  #dimensions D:
  D<-as.matrix(U[,10:16])
  colnames(D)<-c("m","k","s","Kel","mol","amp","cand") 
  rownames(D)<-U[,2]
  if(debug) print(head(D))
  
  #imported V:
  if(ncol(U)>16) {
    V<-as.matrix(U[,17:ncol(U)])
    colnames(V)<-colnames(U)[17:ncol(U)]
  } else V<-matrix(0,nrow=nrow(U),ncol=0)
  rownames(V)<-U[,2]
  
  infosFromU<-list(u_names=u_names,u_roles=u_roles,u_low=u_low,u_high=u_high,u_userunit=u_userunit,u_offset=u_offset,u_gradient=u_gradient,u_transform=u_transform,D=D,V=V)
  return(infosFromU)
}

#importing response information from Z.csv
fetchInfosFromZ<-function(Z,debug=0)
{ 
  
  rownames(Z)<-Z[,2]
  #falls hierbei ein Fehler auftritt: 
  # - "kann Verbindung nicht oeffnen": dann Pfad und Name genau ueberpruefen
  # - "doppelte rownames nicht zulaessig": bitte nach der letzten Zeile in Excel die naechsten Zeilen mit "Strg--" entfernen
  #    denn wahrscheinlich gibt es leere Zeilen, mit denen R nichts anfangen kann
  
  #userunits
  z_userunit<-as.matrix(Z[,3])
  colnames(z_userunit)<-"userunit"
  rownames(z_userunit)<-Z[,2]
  if(debug) print(head(z_userunit))
  
  #offsets
  z_offset<-as.matrix(Z[,5])
  colnames(z_offset)<-"offset" 
  rownames(z_offset)<-Z[,2] 
  if(debug) print(head(z_offset))
  
  #gradients
  z_gradient<-as.matrix(Z[,6])
  colnames(z_gradient)<-"gradient" 
  rownames(z_gradient)<-Z[,2] 
  if(debug) print(head(z_gradient))
  
  #transforms
  z_transform<-as.matrix(Z[,4])
  colnames(z_transform)<-"transform" 
  rownames(z_transform)<-Z[,2] 
  if(debug) print(head(z_transform))
  
  #dimensions D_resp:
  D_resp<-as.matrix(Z[,7:13])
  colnames(D_resp)<-c("m","k","s","Kel","mol","amp","cand") 
  rownames(D_resp)<-Z[,2]
  if(debug) print(head(D))
  
  #imported V_resp:
  if(ncol(Z)>13) {
    V_resp<-as.matrix(Z[,14:ncol(Z)])
    colnames(V_resp)<-colnames(Z)[14:ncol(Z)]
  } else V_resp<-matrix(0,nrow=nrow(Z),ncol=0)
  rownames(V_resp)<-Z[,2]
  
  infosFromZ<-list(z_userunit=z_userunit,z_offset=z_offset,z_gradient=z_gradient,z_transform=z_transform,D_resp=D_resp,V_resp=V_resp)
  return(infosFromZ)
}

provide_u_design<-function(dataList,VR,u_names,u_roles,u_transform,u_gradient,u_offset,u_low,u_high,debug=0,newData=FALSE){
  flag<-"ok"
  factorNames<-u_names
  if(exists("uu_design",where=dataList)&&!newData|exists("uu_new",where=dataList)&&newData) {
    if(!newData) uu_design<-dataList$uu_design else uu_design<-dataList$uu_new
    
    if(qRex<-exists("qualRel",where=dataList)) qualRel<-dataList$qualRel
    #identify how many qual-factors
    qusetInd<-(1:nrow(u_roles))[u_roles=="quset"]		
    
    
    for(qsI in sort(qusetInd,decreasing=TRUE)) {				#qualitative factors are handled in reverse order
      if(debug) print(qsI)
      if(qRex) inqR<-rownames(u_roles)[qsI] %in% rownames(qualRel)
      if(qRex&inqR){			#transfer levels from qalRel
        
        qRelInd<-which(rownames(qualRel)==rownames(u_roles)[qsI])
        #colI<-which(colnames(uu_design)==factorNames[rownames(u_roles)[qsI],])			#column-ind in uu_design, where new columns go
        
        for(cN in names(qualRel[qRelInd,]))  {   #[length(names(qualRel[qRelInd,])):1] ){	#here too: columns are added in reverse order
          if(!cN %in% colnames(uu_design)) {					#then column has not already been generated
            u_roles[cN,1]<-"quset"
            newcol<-rep(0,nrow(uu_design))
            uu_design<-data.frame(uu_design,newcol)
            #uu_design<-data.frame(sCfM(uu_design,1:colI),newcol,sCfM(uu_design,(colI+1):ncol(uu_design))) 
            #colnames(uu_design)[colI]<-factorNames[rownames(u_roles)[qsI],]
            colnames(uu_design)[colnames(uu_design)=="newcol"]<-cN
          }
          for (rowI in 1:nrow(uu_design)) {
            if(debug) print(qsI)
            if(debug) print(paste(substring(factorNames[rownames(u_roles)[qsI],],1,1),"__",as.character(uu_design[,factorNames[rownames(u_roles)[qsI],]])[rowI],sep=""))
            if(paste(substring(factorNames[rownames(u_roles)[qsI],],1,1),"__",as.character(uu_design[,factorNames[rownames(u_roles)[qsI],]])[rowI],sep="") == rownames(u_roles)[qsI])
              uu_design[rowI,cN]<-qualRel[rownames(u_roles)[qsI],cN]
          }
        }
        if(debug) print(qualRel[qRelInd,])
        if(debug) print(qRelInd)
      } else {			#make a column of 0's and 1's for each setting
        if(debug) print(rownames(u_roles)[qsI])
        if(debug) print(factorNames[rownames(u_roles)[qsI],])
        #fnI<-(1:ncol(uu_design))[colnames(uu_design)==factorNames[rownames(u_roles)[qsI],]]
        newcol<-rep(0,nrow(uu_design))
        #names(newcol)<-rownames(u_roles)[qsI]
        for (rowI in 1:nrow(uu_design)) 
          if(paste(substring(factorNames[rownames(u_roles)[qsI],],1,1),"__",as.character(uu_design[,factorNames[rownames(u_roles)[qsI],]])[rowI],sep="") == rownames(u_roles)[qsI])
            newcol[rowI]<-1
        if(debug) print(sum(newcol))
        uu_design<-data.frame(uu_design,newcol) 
        #colnames(uu_design)[fnI]<-factorNames[rownames(u_roles)[qsI],]
        colnames(uu_design)[colnames(uu_design)=="newcol"]<-rownames(u_roles)[qsI]
      }
    }
    #Now the original char-cols can be deleted
    if(length(which(!(colnames(uu_design) %in% rownames(u_roles))))) 
      uu_design<-uu_design[,-which(!(colnames(uu_design) %in% rownames(u_roles)))]
    
    #if levels were transferred from qualRel, u_roles etc are deleted:
    if(qRex){
      remRI<-which(rownames(u_roles) %in% rownames(qualRel))
      VR<-sRfM(VR,-remRI)
      u_names<-sRfM(u_names,-remRI)
      u_roles<-sRfM(u_roles,-remRI)
      u_transform<-sRfM(u_transform,-remRI)
      u_gradient<-sRfM(u_gradient,-remRI)
      u_offset<-sRfM(u_offset,-remRI)
      u_low<-sRfM(u_low,-remRI)
      u_high<-sRfM(u_high,-remRI)
      
      if(!newData){	#bei uu_design werden contr, eqpr eingelesen
        ctrInd<-(1:nrow(u_roles))[u_roles=="contr"|u_roles=="eqpr"|rownames(u_roles) %in% colnames(qualRel)|u_roles=="quset"]
      } else {
        ctrInd<-(1:nrow(u_roles))[u_roles=="contr"|u_roles=="eqpr"|u_roles=="scup"|rownames(u_roles) %in% colnames(qualRel)|u_roles=="quset"]
      }
      u<-uu2u(uu_design[,rownames(u_roles)[ctrInd]],sRfM(u_transform,ctrInd),sRfM(u_gradient,ctrInd),sRfM(u_offset,ctrInd))
    } else { 
      if(!newData){	#bei uu_design werden contr, eqpr eingelesen
        ctrInd<-(1:nrow(u_roles))[u_roles=="contr"|u_roles=="eqpr"]
      } else {
        ctrInd<-(1:nrow(u_roles))[u_roles=="contr"|u_roles=="eqpr"|u_roles=="scup"]
      }
      if(debug) print(colnames(uu_design))
      if(debug) print(colnames(u_roles))
      if(debug) print(colnames(VR))
      if(debug) print(head(uu_design))
      if(debug) print(ctrInd)
      uu_design<-uu_design[,rownames(u_roles)[ctrInd]]								####change20181119
      u<-uu2u(uu_design[,rownames(u_roles)[ctrInd]],sRfM(u_transform,ctrInd),sRfM(u_gradient,ctrInd),sRfM(u_offset,ctrInd))     
    }
    if(u$flag=="ok") u_design<-list(flag=u$flag,u_design=as.matrix(u$u_vec),VR=VR,u_names=u_names,u_roles=u_roles,u_transform=u_transform,
                                    u_gradient=u_gradient,u_offset=u_offset,u_low=u_low,u_high=u_high ) else u_design<-list(flag=u$flag)
  } else u_design<-list(flag="error: no design found to import")
  return(u_design)
}  


#following function is deprecated: covered by provide_u_design with option newData=TRUE
provide_u_new<-function(dataList,u_roles,u_transform,u_gradient,u_offset){
  if(exists("uu_new",where=dataList)) {
    ctrInd<-(1:nrow(u_roles))[u_roles=="contr"|u_roles=="eqpr"|u_roles=="scup"|u_roles=="const"]
    u<-uu2u(dataList$uu_new,sRfM(u_transform,ctrInd),sRfM(u_gradient,ctrInd),sRfM(u_offset,ctrInd))
    if(u$flag=="ok") u_new<-u$u_vec else u_new<-u$flag
  } else if(exists("u_new",where=dataList)) u_new<-dataList$u_new else u_new<-"no new factor settings found to import"
  return(u_new)
}

#function to check if scale-up is possible
checkPotential<-function(D,V,u_roles)
{
  ndlf<-ncol(V)
  Vtop<-sRfM(V,u_roles=="contr"|u_roles=="scup")
  ndes<-qr(Vtop)$rank
  rkD<-qr(D)$rank
  nctr<-sum(u_roles=="contr")
  ncnst<-sum(u_roles=="const")+sum(u_roles=="dep")
  ncnstd<-sum((u_roles=="const"|u_roles=="dep")&apply(D,1,function(x) sum(x^2))!=0)
  
  nctr+ncnstd-rkD
  
  infText<-"Information: number of constants equals dimensional gain."
  if(ncnstd>rkD) infText<-"Warning: const-dominated design."
  if(ncnstd<rkD) infText<-paste("Information: ",rkD-ncnstd," factor(s) available for Scale Up.",sep="")
  if(ndes<ndlf) infText<-paste(infText,"Warning: DL-factors are not independent.")
  
  return(infText)
}


#function prepareProjection for constructing so called "initialization" vectors
#    ... also useful for permuting D and V-matrices
prepareProjection<-function(u_roles,roleList,transpose=FALSE,debug=0,origOrd=FALSE)
  #creates a matrix with a 1 in each row and col according to the role list
{
  colnum<-Reduce("+",lapply(roleList,function(x) sum(u_roles==x)))
  inx<-matrix(0,nrow=nrow(u_roles),ncol=colnum)
  rownames(inx)<-rownames(u_roles)
  if(colnum)colnames(inx)<-1:colnum
  
  if(origOrd){
    sp<-Reduce("|",lapply(roleList,function(x) {u_roles==x}))
    x_ind<-as.array(1:length(sp))[sp]
    colnames(inx)<-rownames(u_roles)[x_ind]
    inx[x_ind,]<-diag(1,nrow=length(x_ind))
    if(debug) print(inx[,])      #has a col for each roleList-factor
  } else {
    begin<-1
    for(x in roleList){
      if(debug) print(x)
      sp<-u_roles==x
      x_ind<-as.array(1:length(sp))[sp]
      if(length(x_ind)) {
        colnames(inx)[begin:(begin-1+length(x_ind))]<-rownames(u_roles)[x_ind]
        inx[x_ind,begin:(begin-1+length(x_ind))]<-diag(1,nrow=length(x_ind))
        if(debug) print(inx[,begin:(begin-1+length(x_ind))])#has a col for each x-factor
      }
      begin<-begin+length(x_ind)
    }
  }
  if(transpose) inx<-t(inx)
  return(inx)
}


#semi-automatic generation of V-matrix based on so called PI-theorem
#i.e particular solution of the homogeneous linear set of equations

suggestVmatrix<-function(u_roles,D,Vin=as.matrix(0),responses=0,debug=0)
{
  if(debug) print(u_roles)  
  if(debug) print(D)  
  if(debug) print(Vin)  
  if(debug) print(responses)  
  if(debug) print(debug)
  DD<-D
  flag<-"ok"
  if(!responses)
  {
    E2<-prepareProjection(u_roles,c("const","quset","scup","contr"),transpose=TRUE) ###changed 180415
  } else { 
    E2<-matrix(0,nrow=nrow(D),ncol=nrow(D))
    rownames(E2)<-rownames(D)
    colnames(E2)<-rownames(D)
    diag(E2)<-1
  }
  
  if(sum(Vin^2)) {
    if(nrow(D)!=nrow(Vin)) flag="error: number of rows in D and Vin is not the same" else
      if(!identical(rownames(D),rownames(Vin))) flag="error: rownames in D and Vin are not the same"
  }
  
  if(debug) print(flag)
  
  r<-qr(t(D)%*%D,tol=1e-07)$rank
  n<-nrow(D)
  
  D<-E2%*%D
  if(debug) print(n)
  if(debug) print(r)
  if(r==n){   # impossible to make things dimensionless: return unity matrix
    flag<-"warning: it is impossible to make things dimensionless: return unity matrix" 
    V<-matrix(0,nrow=n,ncol=n)
    diag(V)<-1
    colnames(V)<-rownames(D) ####changed 180415
    rownames(V)<-rownames(D)
    if(debug) print(flag)
  } else {  
    if(debug) print("restructering rows of D")
    E<-matrix(0,nrow=n,ncol=n)
    diag(E)<-1	#E is a permutation matrix (at the moment just identity)
    
    i<-1
    j<-0
    while(i<=r){
      if(qr(D[(1:i),],tol=1e-07)$rank==i) {
        i<-i+1	#full rank up to now: no need to permute the row away    
      } else {    #row is permuted (not switched!) with row n-j
        cache<-sRfM(D,i)
        D[i:(n-j-1),]<-sRfM(D,(i+1):(n-j))
        rownames(D)[i:(n-j-1)]<-rownames(sRfM(D,(i+1):(n-j)))
        D[n-j,]<-cache
        rownames(D)[n-j]<-rownames(cache)				
        cache<-sRfM(E,i)
        E[i:(n-j-1),]<-sRfM(E,(i+1):(n-j))
        rownames(E)[i:(n-j-1)]<-rownames(sRfM(E,(i+1):(n-j)))
        E[n-j,]<-cache
        rownames(E)[n-j]<-rownames(cache)
        j<-j+1
      }
    }  
    if(debug) print("restructering cols of D")
    
    k<-1
    L<-0
    while(k+L<=7){
      if(qr(D[,(1:k)],tol=1e-07)$rank==k) {
        k<-k+1	#full rank up to now: no need to permute the col away    
      } else {    #col is just skipped
        D<-sCfM(D,-k)
        L<-L+1
      }
    } 
    
    #Now D[1:r,1:r] is invertible and V can be constructed
    #if t(D)%*%V wants to be 0, then t(D_upleft)%*%V_up = -t(D_botleft)%*%V_bot
    if(debug) print("constructing V")
    V<-matrix(0,nrow=n,ncol=n-r)
    colnames(V)<-paste("PI",1:ncol(V),sep="")
    rownames(V)<-rownames(D)
    if(!responses)
      if((r+1)==n) V[(r+1),]<- 1 else diag(V[(r+1):n,])<- 1 else #this tends to put constants into the denominator 
        if((r+1)==n) V[(r+1),]<-1 else diag(V[(r+1):n,])<- 1       #this tends to put responses into the numerator
    V<-rbind(solve(t(sCfM(sRfM(D,1:r),1:r)),-t(sRfM(D,(r+1):n))%*%as.matrix(sCfM(sRfM(V,(r+1):n),1:(n-r)))),sRfM(V,(r+1):n))
    #if(r==1) G<-ginv(cbind(D[(1:r),]))%*%t(D) else G<-ginv(t(D[(1:r),]))%*%t(D)
    #V[1:r,]<--G[,(r+1):n]
    dlCols<-(1:ncol(V))[apply(V,2,function(x) sum(x^2)==1)]	 	####change180415
    if(debug) print(dlCols)
    if(any(dlCols))							####change181119
    {
      dlRows<-sapply(dlCols,function(x) which(V[,x]==1)) 		####change180415
      if(debug) print(dlRows)
      if(length(dlRows)) colnames(V)[dlCols]<-rownames(V)[dlRows]				####change180415
    }
    V<-t(E)%*%V
    V<-t(E2)%*%V
  }
  if(qr(sRfM(V,u_roles=="contr"))$rank < sum(!!apply(sRfM(V,u_roles=="contr"),2,function(x) sum(x^2)))){    #control factors degenerate
    flag="warning: There are more DL-factors than controlled factors: Try to find a DL-factor that is constant"
  }
  #New code that replaces columns in V by Vin and finds constant DL to be put at the end
  if(sum(Vin^2)) {
    if(sum((t(Vin)%*%DD)^2)) { flag<-"warning: Vin is not dimensionless and is ignored" 
    } else {
      repeat{
        Vnew<-cbind(Vin,sCfM(V,-sample(1:ncol(V),ncol(Vin),replace=FALSE)))
        if(qr(Vnew,tol=1e-07)$rank==qr(V,tol=1e-07)$rank) break
      }
      V<-Vnew
    }
  }
  if(debug) print(flag)
  if(debug) print(V)
  Vok<-list(flag=flag,V=round(V,7))
  return(Vok)
} 

#user can redefine cols of the V-matrix by linear combinations
adjustVmatrix<-function(V,columnNo,weights,DLname,debug=0){
  if(debug&&ncol(V)!=length(weights)) print("modify weights")
  if(length(weights)==1&&weights==0) V<-sCfM(V,-columnNo) else {
    V[,columnNo]<-apply(V,1,function(x) sum(x*weights))
    colnames(V)[columnNo]<-DLname
  }
  return(V)
}

############################################################################################################
#three function to handle so called "exponent relaying", when material "constants" depend on state variables
############################################################################################################
#creates  relay matrix from a relay design
makeRelayMatrix<-function(relayDesign,u_roles, u_offset, u_gradient, u_transform,relayCoeffs=0,fromDesign=TRUE)
{
  flag<-"ok"			####change180417
  Rsquare<-0
  if(fromDesign){
    rDcols<-match(colnames(relayDesign),rownames(u_roles))
    rD_roles<-u_roles[rDcols]
    rD_offset<-u_offset[rDcols]
    rD_gradient<-u_gradient[rDcols]
    rD_transform<-u_transform[rDcols]
    
    relayDesign<-t(rD_offset+rD_gradient*t(relayDesign))
    for(i in 1:ncol(relayDesign))
      relayDesign[,i]<-apply(data.frame(relayDesign[,i],rD_transform[i]),1, function(x) {if(x[2]=="none") as.numeric(x[1]) else log10(as.numeric(x[1]))} )
    
    rDcolmeans<-apply(relayDesign,2,mean)
    relayDesign<-t(t(relayDesign)-rDcolmeans)
    
    
    X<-as.matrix(sCfM(relayDesign,rD_roles=="contr"|rD_roles=="dep"))		####change180417
    
    Y<-as.matrix(sCfM(relayDesign,rD_roles!="contr"&rD_roles!="dep")) ## i.e. const ####change180417
    if(!ncol(X)) flag<-"error: no valid columns in the relay design" else   ####change180417
    {												####change180417
      relayCoeffs<-solve(t(X)%*%X)%*%t(X)%*%Y					####change180417
      Rsquare<-apply(Y,2,function(x) {cor(X%*%relayCoeffs,x)^2})		####change180417
    } 											####change180417
  }
  
  if(!flag=="ok")   relayInfo<-list(flag=flag) else   		####change180417
  {											####change180417 - all subsequent rows shifted 2 slots to the right
    relayMatrix<-diag(length(u_roles))
    colnames(relayMatrix)<-rownames(u_roles)
    rownames(relayMatrix)<-rownames(u_roles)
    insertPosV<-unlist(lapply(rownames(relayCoeffs),function(x) {which(x==rownames(relayMatrix))}))
    coeffPosV<-unlist(lapply(colnames(relayCoeffs),function(x) {which(x==rownames(relayMatrix))}))
    for(insertPos in sort(insertPosV,decreasing=TRUE)){		#sorting backward makes things easier!
      if(insertPos<nrow(relayMatrix)) { 
        if(rownames(relayMatrix)[insertPos+1]==paste(colnames(relayMatrix)[insertPos],"__N",sep="")){	#name is already there, coeff is just added
          for(coeffPos in coeffPosV){					#no need to sort this time
            relayMatrix[insertPos,coeffPos]<-relayMatrix[insertPos,coeffPos]+relayCoeffs[rownames(relayMatrix)[insertPos] ,colnames(relayMatrix)[coeffPos]]
            relayMatrix[insertPos+1,coeffPos]<- relayMatrix[insertPos+1,coeffPos]-relayMatrix[insertPos,coeffPos]
          }
        } else {												#row is inserted, coeff entered
          relayMatrix<-rbind(sRfM(relayMatrix,1:insertPos),rep(0,ncol(relayMatrix)),sRfM(relayMatrix,(insertPos+1):nrow(relayMatrix)))
          rownames(relayMatrix)[insertPos+1]<-paste(colnames(relayMatrix)[insertPos],"__N",sep="")
          for(coeffPos in coeffPosV){					#no need to sort this time
            relayMatrix[insertPos,coeffPos]<-relayCoeffs[rownames(relayMatrix)[insertPos] ,colnames(relayMatrix)[coeffPos]]
            relayMatrix[insertPos+1,coeffPos]<- -relayMatrix[insertPos,coeffPos]
          }
        }
      } else {												#row is also inserted, coeff entered
        relayMatrix<-rbind(sRfM(relayMatrix,1:insertPos),rep(0,ncol(relayMatrix)),sRfM(relayMatrix,(insertPos+1):nrow(relayMatrix)))
        rownames(relayMatrix)[insertPos+1]<-paste(colnames(relayMatrix)[insertPos],"__N",sep="")
        for(coeffPos in coeffPosV){					#no need to sort this time
          relayMatrix[insertPos,coeffPos]<-relayCoeffs[rownames(relayMatrix)[insertPos] ,colnames(relayMatrix)[coeffPos]]
          relayMatrix[insertPos+1,coeffPos]<- -relayMatrix[insertPos,coeffPos]
        }
      }
    }
    coeffPosV<-unlist(lapply(colnames(relayCoeffs),function(x) {which(x==rownames(relayMatrix))}))
    for(depPos in coeffPosV){		
      rownames(relayMatrix)[depPos]<-paste(rownames(relayMatrix)[depPos],"__D",sep="")
    }
    relayInfo<-list(flag=flag,relayMatrix=relayMatrix,Rsquare=Rsquare)   		####change180417
  }														####change180417
  return (relayInfo)
}



#imports first and second relay matrices and handles all the naming of columns
relayVmatrix<-function(u_names,u_roles,u_transform,u_gradient,u_offset,u_low,u_high,V,R=0,debug=0){
  flag="ok"
  
  if(length(R)==1) if(R==0) {			#if accidentally called
    R<-diag(nrow(V)) 
    rownames(R)<-rownames(V)
    colnames(R)<-rownames(V)
  }  
  
  if (ncol(R)!=nrow(V)) flag<-"incompatible R and V"
  #apply(R,1,function(x) sum(x^2))-1
  
  V<-R%*%V
  
  new_names<-R%*%cbind(rep(0,ncol(R)))
  new_names[intersect(colnames(R),rownames(R)),]<-u_names[intersect(colnames(R),rownames(R)),]
  new_names[setdiff(rownames(R),colnames(R)),1]<-
    u_names[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_names)<-colnames(u_names)
  
  new_roles<-R%*%cbind(rep(0,ncol(R)))
  new_roles[intersect(colnames(R),rownames(R)),]<-u_roles[intersect(colnames(R),rownames(R)),]
  new_roles[setdiff(rownames(R),colnames(R)),1]<-
    unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {if(unlist(strsplit(x,'__',fixed=TRUE))[length(unlist(strsplit(x,'__',fixed=TRUE)))]=="N") "const" else "dep"}))
  colnames(new_roles)<-colnames(u_roles)
  
  new_transform<-R%*%cbind(rep(0,ncol(R)))
  new_transform[intersect(colnames(R),rownames(R)),]<-u_transform[intersect(colnames(R),rownames(R)),]
  new_transform[setdiff(rownames(R),colnames(R)),1]<-
    u_transform[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_transform)<-colnames(u_transform)
  
  new_gradient<-R%*%cbind(rep(0,ncol(R)))
  new_gradient[intersect(colnames(R),rownames(R)),]<-u_gradient[intersect(colnames(R),rownames(R)),]
  new_gradient[setdiff(rownames(R),colnames(R)),1]<-
    u_gradient[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_gradient)<-colnames(u_gradient)
  
  new_offset<-R%*%cbind(rep(0,ncol(R)))
  new_offset[intersect(colnames(R),rownames(R)),]<-u_offset[intersect(colnames(R),rownames(R)),]
  new_offset[setdiff(rownames(R),colnames(R)),1]<-
    u_offset[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_offset)<-colnames(u_offset)
  
  new_low<-R%*%as.matrix(rep(0,ncol(R)))
  colnames(new_low)<-colnames(u_low)
  new_low[intersect(colnames(R),rownames(R)),]<-u_low[intersect(colnames(R),rownames(R)),]
  new_low[setdiff(rownames(R),colnames(R)),1]<-
    (u_low[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1] + 
       u_high[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1])/2  
  
  new_high<-R%*%as.matrix(rep(0,ncol(R)))
  colnames(new_high)<-colnames(u_high)
  new_high[intersect(colnames(R),rownames(R)),]<-u_high[intersect(colnames(R),rownames(R)),]
  new_high[setdiff(rownames(R),colnames(R)),1]<-
    (u_low[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1] + 
       u_high[unlist(lapply(setdiff(rownames(R),colnames(R)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1])/2  
  
  if(debug) print(cbind(new_low,new_high))
  
  relay<-list(flag=flag,new_names=new_names,new_roles=new_roles,new_transform=new_transform,new_gradient=new_gradient,new_offset=new_offset,new_low=new_low,new_high=new_high,V=V,R=R)
  return(relay)
}

relayV_resp<-function(z_transform,z_gradient,z_offset,V_resp,R_resp,debug=0){
  
  topRight<-matrix(0,nrow<-nrow(R_resp),ncol<-nrow(V_resp)-ncol(R_resp))
  bottomLeft<-matrix(0,nrow<-nrow(V_resp)-ncol(R_resp),ncol<-ncol(R_resp))
  rownames(bottomLeft)<-rownames(V_resp)[(ncol(R_resp)+1):nrow(V_resp)]
  Right<-rbind(topRight,diag(nrow(V_resp)-ncol(R_resp)))
  colnames(Right)<-rownames(V_resp)[(ncol(R_resp)+1):nrow(V_resp)]
  R_resp<-cbind(rbind(R_resp,bottomLeft),Right)
  V_resp<-R_resp%*%V_resp  
  
  new_names<-R_resp%*%cbind(rep(0,ncol(R_resp)))
  new_names[intersect(colnames(R_resp),rownames(R_resp)),]<-z_names[intersect(colnames(R_resp),rownames(R_resp)),]
  new_names[setdiff(rownames(R_resp),colnames(R_resp)),1]<-
    z_names[unlist(lapply(setdiff(rownames(R_resp),colnames(R_resp)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_names)<-colnames(z_names)
  
  new_transform<-R_resp%*%cbind(rep(0,ncol(R_resp)))
  new_transform[intersect(colnames(R_resp),rownames(R_resp)),]<-z_transform[intersect(colnames(R_resp),rownames(R_resp)),]
  new_transform[setdiff(rownames(R_resp),colnames(R_resp)),1]<-
    z_transform[unlist(lapply(setdiff(rownames(R_resp),colnames(R_resp)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_transform)<-colnames(z_transform)
  
  new_gradient<-R_resp%*%cbind(rep(0,ncol(R_resp)))
  new_gradient[intersect(colnames(R_resp),rownames(R_resp)),]<-z_gradient[intersect(colnames(R_resp),rownames(R_resp)),]
  new_gradient[setdiff(rownames(R_resp),colnames(R_resp)),1]<-
    z_gradient[unlist(lapply(setdiff(rownames(R_resp),colnames(R_resp)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_gradient)<-colnames(z_gradient)
  
  new_offset<-R_resp%*%cbind(rep(0,ncol(R_resp)))
  new_offset[intersect(colnames(R_resp),rownames(R_resp)),]<-z_offset[intersect(colnames(R_resp),rownames(R_resp)),]
  new_offset[setdiff(rownames(R_resp),colnames(R_resp)),1]<-
    z_offset[unlist(lapply(setdiff(rownames(R_resp),colnames(R_resp)),function(x) {unlist(strsplit(x,'__',fixed=TRUE))[1]})),1]
  colnames(new_offset)<-colnames(z_offset)
  
  relay_resp<-list(new_names=new_names,new_transform=new_transform,new_gradient=new_gradient,new_offset=new_offset,V_resp=V_resp,R_resp=R_resp)
  return(relay_resp)
}




#function to fetch response specification
fetchV_resp<-function(Z,V,u_roles,import=1,debug=0) 
{
  #Z<-dataList$Z
  rownames(Z)<-Z[,2]
  if(import) {
    V_resp<-as.matrix(sCfM(Z,(14:ncol(Z))))
  } else {
    D<-as.matrix(Z[,7:13])
    rownames(D)<-rownames(Z)
    V_resp<-suggestVmatrix(u_roles,D,responses=nrow(Z),debug=debug)
    respCols<-apply(sRfM(V_resp$V,(nrow(V)+1):nrow(V_resp$V)),2,function(x) sum(x^2))
    V_resp<-sCfM(V_resp$V,(1:ncol(V_resp$V))[respCols==1])
    rownames(V_resp)<-rownames(Z)
    for(i in ((nrow(V)+1):nrow(V_resp)))
      for(j in (1:ncol(V_resp)))
        if(V_resp[i,j]==1) colnames(V_resp)[j]<-paste(colnames(V_resp)[j],rownames(V_resp)[i],sep="_")    
  }
  
  if(debug){print(t(Z[,7:13])%*%V_resp)}             #dimension check
  return(V_resp)
}

#deprecated
#function to fetch RT-response values (untransformed)
#fetchZ<-function(dataList)    
#{
#  z_values<-dataList$z_values
#  return(z_values)
#}


#######################################################################################
#######################################################################################
#
#     Vs and Ws are the link between u-designs for RT-factors and x_design for DL-factors
#
#######################################################################################
#######################################################################################

supplyVsAndWs<-function(VR,u_roles,keepCols=FALSE,whichScup=0,eps=0.00001,debug=0){
  #there may be three problems with inverting V (or VR):
  #(1) one is that due to relaying there may be more DL-factors needed than RT-factors available;
  #     then either the most important (DLF-importance) or a user specific (keepCols) are kept
  #(2) another one is that the resulting V (or VD matrix is not of full rank
  #     the the user has badly chosen his DL-factors and must redefine 
  #(3) another is that the scup factor is collinear with the design-DL-factors
  #     then either the most suitable (cosines) or a user alternative (which) are to be used
  
  eqprIndex<-(1:length(u_roles))[u_roles=="eqpr"]  #remember which factor is eqpr
  u_roles[eqprIndex]<-"contr"    #for now: don't differentiate between "contr" and "eqpr" here
  scupIndex<-(1:length(u_roles))[u_roles=="scup"]  #remember which factor is scup - needed later
  u_roles[scupIndex]<-"contr"    #for now: don't differentiate between "contr" and "scup" here
  nctrR<-sum(u_roles=="contr"|u_roles=="quset")			 	####change180415
  
  if(ncol(VR)>nctrR)  {                    #if too many DL-factors then slack factor(s) introduced
    if(keepCols[1]) {                      #either user wish
      VDS<-sCfM(VR,keepCols) 
      VSL<-sCfM(VR,setdiff(1:ncol(VR),keepCols))  
    } else {                               #least important
      DLFimportance<-apply(VR[u_roles=="contr"|u_roles=="quset",],2,function(x) sum(x^2))  	 	####change180415
      
      VDS<-sCfM(VR,sort(DLFimportance,decreasing=TRUE,index.return=TRUE)$ix[1:nctrR])
      VSL<-sCfM(VR,sort(DLFimportance,decreasing=TRUE,index.return=TRUE)$ix[(nctrR+1):ncol(VR)])
    } 
  } else {
    VDS<-VR
    VSL<-VR[,NULL]                         #no slack factor needed
  }
  
  if(!ncol(VR)||qr(sRfM(VDS,u_roles=="contr"|u_roles=="quset"))$rank < ncol(VDS)){    #control factors degenerate   	 	####change180415
    
    flag="Please reduce, change or supply DL-factors, the exponent matrix, V, is singular"
    
    VsAndWs<-list(flag=flag)
  } else {                   #the rest of the function is in this else-branch
    flag<-"ok"
    x_roles<-matrix(0,ncol=1,nrow=nrow(u_roles))
    colnames(x_roles)<-"roles"
    x_roles[1:ncol(VDS)]<-"dldes"
    in_const<-prepareProjection(u_roles,"const")
    in_dep<-prepareProjection(u_roles,"dep")
    in_constDep<-prepareProjection(u_roles,c("const","dep"),origOrd=TRUE)
    in_contr<-prepareProjection(u_roles,c("contr","quset"),origOrd=TRUE)    	 	####change180415
    
    #VD
    VDtop<-as.matrix(sRfM(VDS,u_roles=="contr"|u_roles=="quset"))    	 	####change180415
    
    VDbot<-as.matrix(sRfM(VDS,u_roles=="const"|u_roles=="dep"))
    
    VDB<-cbind(VDS,in_constDep)
    if(sum(u_roles=="const")+sum(u_roles=="dep")) x_roles[(ncol(VDS)+1):ncol(VDB)]<-"const"
    #WD
    WDtop<-solve(t(VDtop)%*%VDtop)%*%t(VDtop)        #upper left corner
    WDbot<-(-VDbot)%*%WDtop                          #lower left corner
    WDS<-rbind(WDtop,WDbot)                          #join them
    if(sum(u_roles=="const")+sum(u_roles=="dep")){
      WDB<-cbind(WDS,t(VDB)[,u_roles=="const"|u_roles=="dep"])        #add right side
      WDB[1:nrow(WDtop),(ncol(WDtop)+1):ncol(WDB)]<-0  #top right to 0
      if(debug){
        print(rownames(rbind(t(in_contr),t(in_constDep))))
        print(colnames(WDB))
        print("muessten gleich sein")
      }
      WDB<-WDB%*%rbind(t(in_contr),t(in_constDep))        #sort its columns correctly
    } else WDB<-WDS
    WDW<-sRfM(WDB,1:nrow(WDtop))                     #select top rows
    
    if(debug)
    {
      print(VDS%*%WDW)       
      print(round(WDW%*%VDS,7))    #muesste kleine ID-Matrix sein!
      print("muesste kleine ID-Matrix sein!")
      print(VDB%*%WDB)
      print(round(WDB%*%VDB,7))    #muesste groessere ID-Matrix sein!
      print("muesste groessere ID-Matrix sein!")
    }
    
    cosines<-matrix(0,nrow(VDtop))    #at this point projection is on VR-space
    cosines<-round(as.matrix(diag(VDtop%*%solve(t(VDtop)%*%VDtop)%*%t(VDtop))),7)
    rownames(cosines)<-rownames(VDtop)
    colnames(cosines)<-"cosines"
    
    
    
    u_roles[scupIndex]<-"scup" 
    scupFactors<-suggestScupFactors(u_roles,cosines,VDtop,WDtop,debug=debug)
    in_pot_scup<-scupFactors$in_scup
    scupCosines<-scupFactors$scupCosines
    if(debug){ 
      print("scupFactorsFlag says: ")
      print(scupFactors$flag)
    }
    if(nrow(VDB)==ncol(VDB)) {    #then RES is empty and scale up is not possible
      RES<-sCfM(in_contr,NULL)
      RES_scup<-RES
      in_in_scup<-in_pot_scup
      if(length(u_roles[u_roles=="scup"])) {
        u_roles[u_roles=="scup"]<-"contr" 
        flag<-"warning: No scale up is possible at all"
        #in_in_scup<-in_scup
      }
    } else {                                      # RES needs to be constructed          
      if(length(u_roles[u_roles=="scup"])==1) {	#if a scup-factor has been given, ...
        if(cosines[u_roles=="scup"] > 1-eps ) {   # ... with cosine 1, it is rejected
          u_roles[u_roles=="scup"]<-"contr"   
          flag<-paste("Scale up is impossibe for", rownames(u_roles)[u_roles=="scup"])
        } else {                                  # ...with cosine < 1 ist, it is accepted
          in_in_scup<-sCfM(in_pot_scup,(colnames(in_pot_scup)==rownames(u_roles)[u_roles=="scup"]))
          whichScup<-(1:ncol(in_pot_scup))[(colnames(in_pot_scup)==rownames(u_roles)[u_roles=="scup"])]
        }
      } else {
        if (whichScup){
          in_in_scup<-sCfM(in_pot_scup,whichScup)
          if(debug) print(in_in_scup)
          u_roles[(1:length(u_roles))[in_in_scup[,1]==1]]<-"scup"
        } else {
          in_in_scup<-sCfM(in_pot_scup,whichScup)      
        }
      } 
      #jetzt sind in_in_scup, whichScup und u_roles konsistent gesetzt, aber RES ist noch nicht berechnet            
      nscup<-ncol(in_in_scup)
      if (nscup){ 
        RES_scup<-in_in_scup-(VDB%*%WDB)%*%in_in_scup
        #       RES_scup<-in_in_scup-(VDS%*%WDW)%*%in_in_scup
        P_scup<-RES_scup%*%solve(t(RES_scup)%*%RES_scup)%*%t(RES_scup)
        x_roles[ncol(VDB)+1]<-"scup"
        if(debug) print(qr(t(RES_scup)%*%RES_scup,eps)$rank)
      } else { 
        RES_scup<-in_in_scup
        P_scup<-0
      } 
      RES<-RES_scup
    }
    
    #eqpr-factors with corresponding RES: if possible user suggested eqprInd-factors are used, otherwise highest cosines  
    if(ncol(RES_scup)==nrow(VDB)-ncol(VDB)) {    #no eqpr possible, RES is commplete
      flag<-paste(flag,"info: equivalence check is not possible", sep=" & ")
      in_eqpr<-sCfM(in_pot_scup,NULL)
      RES<-RES_scup
    } else {
      in_eqpr<-sCfM(in_pot_scup,(1:ncol(in_pot_scup))!=whichScup)
      in_in_eqpr<-sCfM(in_eqpr,intersect(rownames(u_roles)[eqprIndex],colnames(in_eqpr)))
      colnames(in_in_eqpr)<-intersect(rownames(u_roles)[eqprIndex],colnames(in_eqpr))
      nAddeqpr<-nrow(VDB)-ncol(VDB)-ncol(RES_scup)-ncol(in_in_eqpr)
      if(nAddeqpr>0) {
        i<-0
        for(Ind in 1:nAddeqpr) {   #look for a good RES-direction
          in_in_eqpr<-cbind(in_in_eqpr,sCfM(in_eqpr,(Ind+i)))
          RES_eqpr<-in_in_eqpr-(VDB%*%WDB+P_scup)%*%in_in_eqpr
          while(qr(t(RES_eqpr)%*%RES_eqpr,tol=eps)$rank!=(ncol(in_in_eqpr))||sum(RES_eqpr^2)<eps^2)
          { 
            i<-i+1
            in_in_eqpr<-sCfM(in_in_eqpr,-ncol(in_in_eqpr))
            in_in_eqpr<-cbind(in_in_eqpr,sCfM(in_eqpr,Ind+i))
            RES_eqpr<-in_in_eqpr-(VDB%*%WDB+P_scup)%*%in_in_eqpr
          } #endwhile
        } #endfor
        
        in_eqpr<-in_in_eqpr
        RES_eqpr<-in_eqpr-(VDB%*%WDB+P_scup)%*%in_eqpr
        x_roles[(nrow(VDB)-ncol(in_eqpr)+1):nrow(VDB)]<-"eqpr"
      }  #endif
      
      #constructing RES-matrix 
      RES<-cbind(RES_scup,RES_eqpr) 
    } #end-else                
    
    if(ncol(RES)!=nrow(VDB)-ncol(VDB)) {
      if(debug) print(c(ncol(RES),nrow(VDB),ncol(VDB)))
      flag<-"internal error: RES-matrix has wrong dimension"
    }
    
    if(ncol(in_eqpr)) for(col in (1:ncol(in_eqpr))) u_roles[(1:length(u_roles))[in_eqpr[,col]==1] ]<-"eqpr"
    in_contr<-prepareProjection(u_roles,c("contr","quset"))    	 	####change180415
    
    
    VRES<-as.matrix(cbind(VDB,RES))
    rownames(x_roles)<-colnames(VRES)
    
    if(ncol(RES)) WRES<-as.matrix(rbind(WDB,solve(t(RES)%*%RES)%*%t(RES))) else
      WRES<-WDB
    if(debug) print(round(WRES%*%VRES,7))
    VsAndWs<-list(u_roles=u_roles,x_roles=x_roles,VRES=VRES,VSL=VSL,WRES=WRES,in_contr=in_contr,in_const=in_const,in_dep=in_dep,in_constDep=in_constDep,in_scup=in_in_scup,in_eqpr=in_eqpr,cosines=cosines,flag=flag)
  }
  return(VsAndWs)
}


#function fetchXLevels
fetchXLevels<-function(x_levels,x_low,x_high,import=1,debug=0,eps=1e-6)
{
  #importing logs of (!!) low and high levels for DLF-factors
  
  x_D_high<-x_high
  colnames(x_D_high)<-"DL-des-high"
  x_D_low<-x_low
  colnames(x_D_low)<-"DL-des-low"
  
  if(import) {
    #x_levels<-x_levels
    x_D_low[1:nrow(x_levels)]<-apply(cbind(as.matrix(x_levels[,1]),x_low),1,max)
    x_D_high[1:nrow(x_levels)]<-apply(cbind(as.matrix(x_levels[,2]),x_high),1,min)
    colnames(x_D_high)<-"DL-des-high"
    if(debug)print(cbind(x_high,x_D_high,(x_high+x_low)/2,x_D_low,x_low))
  }
  
  sumOk<-sum(x_high>=x_D_high&x_D_high>=(x_high+x_low)/2-eps&(x_high+x_low)/2+eps>=x_D_low&x_D_low>=x_low)==length(x_low)  #TRUE would be cool!
  if(debug) if(!sumOk) print(cbind(x_high>=x_D_high,x_D_high>(x_high+x_low)/2,(x_high+x_low)/2>x_D_low,x_D_low>=x_low))
  
  XLevels<-list(sumOk=sumOk,x_D_low=x_D_low,x_D_high=x_D_high)
  return(XLevels)
}


######################################################################################
######################################################################################
#              FetchDesign: wrapper for all design variants 
######################################################################################
#      5 variants
#      Variation 1: generates a u_design i.e. design in RT-space
#      Variation 2: imports a u_design in  RT-space; uu_design's it should be preceeded by a call to provide_u_design() 
#      Variation 3: imports a scaled and centred x_design in  DLF-space (-1 +1 -array)
#      Variation 4: generates x_design i.e. design in  DLF-space
#      Variation 5: generates a D_optimal x_design s.t. u_design satisfies u_low/high levels
######################################################################################
######################################################################################

fetchDesign<-function(variant,VRES,WRES,toBeImported,u_low,u_high,x_low,x_high,u_roles,x_roles,test_eqpr,desType,desParameters,debug=0)
  #for variant 1: sc design generated for contr- and eqpr-factors, positioned, const and scup factors added and returned
  #for variant 2: complete design imported, checked for correctness, returned
  #for variant 3: complete DLF and eqpr-factor design imported, positioned, const and scup factors added, slider technique applied returned
  #for variant 4: design for DLF-factors and eqpr-factors generated, positioned, const and scup factors added, slider technique applied, returned
  #for variant 5: as for variant 4
  #return values of variants are u_design and flag for var 1&2, u_design and flas for 3&4
  #       const-columns have not been added!
  #return values of this wrapper are x_design, u_design and flag for error messages or problems
  #order of columns in designs:
  #u_designs: order as in U.csv, in_xxxx determin positioning
  #x_designs: design DL-factors at left, then const-factors, then scup-factors, then eqpr-factors - this can be checked in x_roles
{
  in_scup<-prepareProjection(u_roles,"scup")
  in_eqpr<-prepareProjection(u_roles,"eqpr")
  in_const<-prepareProjection(u_roles,"const")
  in_constDep<-prepareProjection(u_roles,c("const","dep"),origOrd=TRUE)
  in_dep<-prepareProjection(u_roles,"dep")
  in_contr<-prepareProjection(u_roles,c("contr","quset"),origOrd=TRUE)  ####change180424
  dleqInd<-(1:ncol(VRES))[x_roles=="dldes"|x_roles=="eqpr"]
  dlInd<-(1:ncol(VRES))[x_roles=="dldes"]
  if(debug) print(ncol(in_scup))
  if(debug) print(ncol(in_eqpr))
  if(debug) print(ncol(in_const))
  if(debug) print(ncol(in_constDep))
  if(debug) print(ncol(in_contr))
  
  if(variant==1){
    modde<-createDesigninRFspace(u_low,u_high,in_contr,in_eqpr,dlInd,dleqInd,desType,desParameters,debug)
    
  } else if(variant==2){
    modde<-importDesigninRFspace(toBeImported,in_contr,in_eqpr,dlInd,dleqInd,debug)
    
  } else if(variant==3){
    modde<-importSCDesigninDLFspace(toBeImported,x_low,x_high,in_eqpr,dlInd,dleqInd,test_eqpr,debug)
    
  } else if(variant==4){
    modde<-createDesigninDLFspace(x_low,x_high,in_eqpr,dlInd,dleqInd,test_eqpr,desType,desParameters,debug)
    
  } else if(variant==5){
    modde<-createUconstrainedDoptimalinDLFspace(toBeImported,VRES,WRES,u_low,u_high,x_low,x_high,x_roles,in_constDep,in_scup,in_eqpr,dlInd,dleqInd,test_eqpr,desType,desParameters,debug)
    
  } else flag<-"variant does not exist"
  flag<-modde$flag
  if(debug) print(flag)
  if(flag=="ok") {
    constDepScupDes<-constDepScupDesign(u_low,u_high,nrow(modde$design),in_constDep,in_scup,useLowForScup=TRUE)
    if(variant==1 | variant==2){			  #i.e. u_design
      u_design<-cbind(modde$design,constDepScupDes)[,rownames(VRES)]
      dimDoeMod<-modde$model
      Dfrml<-transDimDoeMod(dimDoeMod,x_roles)
    } else if(variant==3 | variant==4){ 	  #i.e. x_design - scup factors treated using slider technique
      x_design<-cbind(modde$design,constDepScupDes)[,colnames(VRES)]   #careful col for scup factor not yet correct!!
      x_cp<-ccW(x_low,x_high,0)$cp
      SW<-ccW(x_low,x_high,0)$SW
      dimDoeMod<-modde$model
      Dfrml<-transDimDoeMod(dimDoeMod,x_roles)
      if(ncol(in_scup)){
        selectInd<-ncol(VRES)-ncol(in_eqpr) 
        RES<-sCfM(VRES,selectInd)
        RESred<-t(in_scup)%*%RES     #use slider technique here
        if(debug) print(RESred)
        MR<-x_design%*%WRES
        SLI<-(sCfM(constDepScupDes,ncol(constDepScupDes))-MR%*%in_scup)%*%solve(t(RESred))
        u_design<-MR+SLI%*%t(RES)
      } else { 
        u_design<-x_design%*%WRES
      }
    } else if(variant==5) {
      u_design<-modde$design%*%WRES
      dimDoeMod<-modde$dimDoeMod
      Dfrml<-transDimDoeMod(dimDoeMod,x_roles)
    }
    x_design<-u_design%*%VRES
    if(debug) print(head(x_design))
    orSC<-orScaleDesign(x_design)
    #newXLevels<-recalcLevels(x_design)
    #x_D_low<-newXLevels$D_low
    #x_D_high<-newXLevels$D_high
    
    #x_D_cp<-ccW(x_D_low,x_D_high,0)$cp
    #SW<-ccW(x_D_low,x_D_high,0)$SW
    #MD<-round(scaleDesign(x_design,x_D_cp,SW),5)
    designInfos<-list(u_design=u_design, x_D_low=orSC$D_low, x_D_high=orSC$D_high, MD=orSC$MD, Dfrml=Dfrml,dimDoeMod=dimDoeMod,flag=flag)
  } else designInfos<-list(flag=flag)         #i.e. flag is not "ok"
  
  return(designInfos)
}


u2uu<-function(u_vec,u_transform,u_gradient,u_offset){
  flag<-"ok"
  uu_vec<-u_vec
  if(length(colnames(u_vec))!=length(rownames(u_transform)))  flag<-"inconsistent RF-number" else
    if(sum(colnames(u_vec)==rownames(u_transform))!=ncol(u_vec)) flag<-"inconsistent RF-names" else
      for(i in 1:ncol(u_vec))  {
        if(u_transform[i]=="log") uu_vec[,i]<-10^u_vec[,i]
        uu_vec[,i]<-(uu_vec[,i]-u_offset[i])/u_gradient[i]
      }
  uu<-list(flag=flag,uu_vec=uu_vec)
  return(uu)
}

uu2u<-function(uu_vec,u_transform,u_gradient,u_offset,debug=0){
  flag<-"ok"
  u_vec<-uu_vec
  if(length(colnames(uu_vec))!=length(rownames(u_transform)))  flag<-"inconsistent RF-number - check input column number" else
    if(sum(colnames(uu_vec)==rownames(u_transform))!=ncol(u_vec)) flag<-"inconsistent RF-names - check input column names" else
      for(i in 1:ncol(uu_vec))  {
        u_vec[,i]<-u_offset[i]+u_gradient[i]*uu_vec[,i]   #so to invert: uu<-(u-offset)/gradient
        if(u_transform[i]=="log") u_vec[,i]<-log10(u_vec[,i])
        # s. auch Umrechnung von u_low und u_high - dort ist es (noch) anders:
        # u_low<-apply(data.frame(u_low,u_transform),1, function(x) {if(x[2]=="none") as.numeric(x[1]) else log10(as.numeric(x[1]))} )
      }
  if(debug) print(flag)
  if(debug) print(colnames(u_vec))
  if(debug) print(rownames(u_transform))
  if(debug) print(colnames(uu_vec)==rownames(u_transform))
  u<-list(flag=flag,u_vec=as.matrix(u_vec))			#### as.matrix added 20181119
  return(u)
}


#####################################################################################
#####################################################################################
#
#	Get responses and fit and diagnose model using lm() and R-formula-model
#
#####################################################################################
#####################################################################################

#function to transfer a DimDoe-model to a R-formula-model and vice versa
transDimDoeMod<-function(ModelExpr,x_roles,forward=TRUE)
{
  if(forward) # transcribe from DimDoeMod to R-formula
  {
    frml<-"~1" 
    for(term in ModelExpr)
    { 
      if(length(term)==1) {
        if(term!=0&&x_roles[term]!="const")            #linear term (not model-const and not const DL-factor)
          frml<-paste(frml,rownames(x_roles)[term],sep="+")
      } else if (length(term)==2&&term[1]==term[2]&&x_roles[term[1]]!="const") {     #quadratic term
        frml<-paste(frml,"+I(",rownames(x_roles)[term[1]],"^2)",sep="")
      } else if (length(term)==2&&term[1]<term[2]&&x_roles[term[1]]!="const"&&x_roles[term[2]]!="const")  {
        frml<-paste(frml,"+I(",rownames(x_roles)[term[1]],"*",rownames(x_roles)[term[2]],")",sep="")
      }
    }
    ret<-frml 
  } else {     # transcribe from R-formula to DimDoeMod
    ModTermList<-unlist(strsplit(ModelExpr,"[+]"))
    dimDoeMod<-list(0)	                         #set the constant
    for(term in ModTermList[-1])
    {
      if(term%in%rownames(x_roles))          #linear terms
        dimDoeMod<-append(dimDoeMod,(1:nrow(x_roles))[rownames(x_roles)==term])
      else {                                  #interaction or quadratic
        termp<-substr(term,3,nchar(term)-1)
        terms<-gsub(" +","",unlist(strsplit(termp,"[*]")))
        if(substr(terms[1],(nchar(terms[1])-1),nchar(terms[1]))!="^2") {
          dimDoeMod[[length(dimDoeMod)+1]]<-c((1:nrow(x_roles))[rownames(x_roles)==terms[1]],(1:nrow(x_roles))[rownames(x_roles)==terms[2]])
        } else {
          termsa<-substr(terms[1],1,nchar(terms[1])-2)
          dimDoeMod[[length(dimDoeMod)+1]]<-c((1:nrow(x_roles))[rownames(x_roles)==termsa],(1:nrow(x_roles))[rownames(x_roles)==termsa])
        } 
      }
    }
    ret<-dimDoeMod     
  }
  return(ret)
}

prepareAnalysis<-function(u_design,z,z_transform,z_gradient,z_offset,VR,V_resp,debug=0)  
{
  flag="ok"
  nr<-nrow(VR)
  nt<-nrow(V_resp)
  if((nt-nr)!=ncol(z)) {
    flag<-"inconsistency in the number of responses"
    preparedAnalysis<-flag
  } else {
    zList<-uu2u(z,sRfM(z_transform,(nr+1):nt),sRfM(z_gradient,(nr+1):nt),sRfM(z_offset,(nr+1):nt),debug=debug)                          	#z in user units gets transformed to SI-units and logs
    z<-zList$u_vec
    flag<-zList$flag
    if(debug) print(flag)
    if(flag!="ok") {
      fittedModel<-flag
      if(debug) print(flag)
    } else {
      DL_resp<-fetchDL_resp(u_design,z,nr,V_resp)   			
      flag<-DL_resp$flag
      if(debug) print(flag)
    }
    if(flag!="ok") {
      preparedAnalysis<-flag
    } else {
      y<-DL_resp$y
      preparedAnalysis<-list(flag=flag,y=y)
    }
  }
  return(preparedAnalysis)
}


#function that calls lm( ) to fit the models 
#returns a list containing an lm-object for each response(i.e. list length is no of responses)
#DmodelTempateList is a list of R-formula-models, one for each response
doAnalysis<-function(DmodelTemplateList,MD,y,debug=0){
  if(nrow(y)!=nrow(MD)){
    flag="design runs don't match response values"
    modelList<-list(flag=flag)
  } else {
    YMD<-data.frame(y,MD)
    if(length(colnames(y))!=length(DmodelTemplateList)) #then all models as the first
      for(Ind in (1:length(colnames(y)))) DmodelTemplateList[[Ind]]<-DmodelTemplateList[[1]]
      modelList<-list()
      for(DL_resp in as.list(colnames(y))) {
        FmodelTemplate<-paste(DL_resp,DmodelTemplateList[[length(modelList)+1]],sep="")
        modelList[[length(modelList)+1]]<-lm(FmodelTemplate,data=YMD)
      }
      names(modelList)<-colnames(y)
  }
  if(debug) print(flag) 
  return(modelList)
}

#the following function is deprecated 31.03.2018
analyseDesign<-function(DmodelTemplateList,u_design,z,z_transform,z_gradient,z_offset,VR,V_resp,MD,debug=0)  
  #DmodelTempateList is a list of R-formula-models, one for each response
{
  print("please replace call to analyseDesign by 2 calls, 
        one to prepareAnalysis and on to doAnalysis")
  flag="ok"
  nr<-nrow(VR)
  nt<-nrow(V_resp)
  if((nt-nr)!=ncol(z)) {
    flag<-"inconsistency in the number of responses"
    modelList<-flag
  } else {
    zList<-uu2u(z,sRfM(z_transform,(nr+1):nt),sRfM(z_gradient,(nr+1):nt),sRfM(z_offset,(nr+1):nt))                          	#z in user units gets transformed to SI-units and logs
    z<-zList$u_vec
    flag<-zList$flag
    if(debug) print(flag)
    if(flag!="ok") {
      fittedModel<-flag
      if(debug) print(flag)
    } else {
      DL_resp<-fetchDL_resp(u_design,z,nr,V_resp)   			
      flag<-DL_resp$flag
      if(debug) print(flag)
    }
    if(flag!="ok") {
      modelList<-flag
    } else {
      y<-DL_resp$y
      YMD<-data.frame(y,MD)
      if(length(colnames(y))!=length(DmodelTemplateList)) #then all models as the first
        for(Ind in (1:length(colnames(y)))) DmodelTemplateList[[Ind]]<-DmodelTemplateList[[1]]
      modelList<-list()
      for(DL_resp in as.list(colnames(y))) {
        FmodelTemplate<-paste(DL_resp,DmodelTemplateList[[length(modelList)+1]],sep="")
        modelList[[length(modelList)+1]]<-lm(FmodelTemplate,data=YMD)
      }
      names(modelList)<-colnames(y)
    }
  }
  if(debug) print(flag) 
  return(modelList)
}


#this function allows prediction with the fitted model to a new dataset
#if validate is TRUE new y-values are imported
scaleUpProcess<-function(u_new,z_new,fittedModel,u_roles,u_low,u_high,VRES,V_resp,x_low,x_high,validate=TRUE,debug=0)
{
  flag="ok"
  #u_new<-dataList$u_new
  #u_new<-log10(u_new)		#careful from now on u is in log
  in_const<-prepareProjection(u_roles,"const")
  in_dep<-prepareProjection(u_roles,"dep") 
  in_constDep<-prepareProjection(u_roles,c("const","dep"),origOrd=TRUE)   
  in_scup<-prepareProjection(u_roles,"scup")
  constScupNew<-constDepScupDesign(u_low,u_high,nrow(u_new),in_constDep,in_scup,useLowForScup=FALSE)
  u_new<-cbind(u_new,constScupNew)[,rownames(VRES)][,rownames(VRES)]
  x_new<-u_new%*%VRES
  if(debug) print(head(x_new))
  x_D_cp<-ccW(x_low,x_high,0)$cp
  SW<-ccW(x_low,x_high,0)$SW
  x_new_sc<-data.frame(scaleDesign(x_new,x_D_cp,SW))
  if(debug) print(paste("stelle 0:",head(x_new_sc),sep=" "))
  y_new_pred<-matrix(0,nrow=nrow(x_new),ncol=length(fittedModel))
  colnames(y_new_pred)<-colnames(V_resp)
  for(modInd in (1:length(fittedModel))) {
    y_new_pred[,modInd]<-predict(fittedModel[[modInd]],x_new_sc) 
  }
  if(debug) print(paste("stelle1:",head(y_new_pred),sep=" "))
  if(validate) {
    #z_new<-dataList$z_new
    if(debug) print(paste("stelle2:",head(z_new),sep=" "))
    l_z<-log10(z_new)
    if(debug) print(paste("stelle3:",head(l_z),sep=" "))
    DL_resp<-fetchDL_resp(u_new,l_z,nrow(VRES),V_resp)
    flag=DL_resp$flag
    if(flag!="ok") {
      validatedModel<-list(flag=flag)
    } else { #flag=="ok"
      y_new<-DL_resp$y
      if(debug) print(head(y_new))
      y_new_res<-y_new-y_new_pred
      RSD_new<-sqrt(apply(y_new_res^2,2,sum)/nrow(y_new))
      RSD_new_PCNT<-(10^RSD_new-1)*100
      R2Y_new<-1-apply(y_new_res^2,2,sum)/apply(sweep(y_new,2,apply(y_new,2,mean))^2,2,sum)
      l_z_new_pred<-y_new_pred%*%ginv(V_resp[(nrow(VRES)+1):nrow(V_resp),])-u_new%*%V_resp[1:nrow(VRES),]
      z_new_pred<-10^l_z_new_pred
      validatedModel<-list(u_new=u_new,x_new=x_new,x_new_sc=x_new_sc,y_new=y_new,y_new_pred=y_new_pred,y_new_res=y_new_res,R2Y_new=R2Y_new,RSD_new=RSD_new,RSD_new_PCNT=RSD_new_PCNT,z_new=z_new,z_new_pred=z_new_pred,flag=flag)
    } 
  } else {  #!validate
    l_z_new_pred<-y_new_pred%*%ginv(V_resp[(nrow(VRES)+1):nrow(V_resp),])-u_new%*%V_resp[1:nrow(VRES),]
    z_new_pred<-10^l_z_new_pred
    if(debug) print(head(y_new))
    flag="ok for prediction, no validation"
    validatedModel<-list(u_new=u_new,x_new=x_new,x_new_sc=x_new_sc,y_new_pred=y_new_pred,z_new_pred=z_new_pred,flag=flag)
  }
  return(validatedModel)
}    

#########################################################################################
#########################################################################################
#
#	Get responses and fit and diagnose model using just motrix-operations and dimDoeMod
#
#	must be improved to accomodate different models for different responses
#	... and to accomodate PLS, resp. OPLS
#
#########################################################################################
#########################################################################################


fitAndXvalDesign<-function(model,u_design,z,VR,V_resp,MD,debug=0){
  DL_resp<-fetchDL_resp(u_design,z,nrow(VR),V_resp)   #log!!
  flag=DL_resp$flag
  if(flag!="ok") {
    fittedModel<-flag
  } else {
    y<-DL_resp$y
    EMD<-extend_x(MD,model)
    I <- solve(t(EMD)%*%EMD)   #Informationsmatrix XTX^-1
    g <- diag(I)           #Diagonale der I-Matrix
    b <- I%*%t(EMD)%*%y      #Modellkoeffizienten b
    H <- EMD%*%I%*%t(EMD)      #Hat-Matrix zur Berechnung der Prognosen
    h <- diag(H)	         #Diagonale der H-Matrix (leverages)
    p <- H%*%y             #Prognosen p
    e <- y - p             #Residuen    e
    r <- e/(1-h)           #Studentisierte Residuen r
    R2b<-1-apply(e^2,2,sum)/apply(sweep(y,2,apply(y,2,mean))^2,2,sum)
    R2b		#by 1 - residual/total corrected
    Q2b<-1-apply(r^2,2,sum)/apply(sweep(y,2,apply(y,2,mean))^2,2,sum)
    Q2b		#by 1 - press/total corrected
    
    RSD<-sqrt(apply(r^2,2,sum)/(nrow(MD)-length(b[,1])))
    est<-e/RSD
    RSDPCNT<-100*(10^RSD-1)
    fittedModel<-list(flag=flag,b=b,p=p,R2b=R2b,Q2b=Q2b,RSD=RSD,RSDPCNT=RSDPCNT)
  }
}
