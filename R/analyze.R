.First.lib <- function(lib, pkg) {
  library.dynam("LN3GV", pkg, lib)
#  library.dynam("product", pkg, lib)
}

analyze<-function(dat,pats,method,log.scale=TRUE,MaxIt=30,parm.conv=.005,prob.conv=.0005,parm.init=NULL,prob.init=NULL){

if(sum(method==c("LNNMV*","LNNGV","LN3","LN3MV*","LN3GV"))==0) stop("Error: method argument must be specified as one of
LNNMV*, LNNGV, LN3, LN3MV*, or LN3GV")


if(log.scale==FALSE){
if(min(dat)<=0) stop("Error: Data contains non-positive values. Using log.scale=FALSE produces NaN's")
 dat<-log(dat)
}

if(sum(is.nan(dat))!=0) stop("Error: Data contains NaN's")
if(sum(is.na(dat))!=0) stop("Error: Data contains NA's")
if(ncol(dat)!=ncol(pats)) stop("Error: Number of columns in data matrix does not match number of columns in expression pattern matrix")

Gn<-apply(pats,1,max)


## Reformat pats to use smallest integers possible  
## (speeds up processes and fixes pattern matrices that aren't specified with integers)
for(i in 1:nrow(pats)){
x<-1:length(unique(pats[i,]));names(x)<-as.character(unique(pats[i,]))
pats[i,]<-x[as.character(pats[i,])]
}

## Infer treatment structure from columns of pattern matrix 
## (Experimental units that have identical columns are assumed to be replicates when finding MSE for each gene)
#if(is.null(trt)){
x<-1:length(unique(apply(pats,2,paste,collapse="")));names(x)<-unique(apply(pats,2,paste,collapse=""))
trt<-x[apply(pats,2,paste,collapse="")]; names(trt)<-NULL
print(trt)
#}

patmat<-vector("list",nrow(pats))
for(i in 1:nrow(pats)){
patmat[[i]]<-matrix(0,ncol(pats),ncol(pats))
for(j in 1:ncol(pats)){
for(k in 1:j){
if(pats[i,j]==pats[i,k]){
patmat[[i]][j,k]<-1
patmat[[i]][k,j]<-1
} } } }



if(method=="LN3MV*"|method=="LN3GV"){

##### LN3MV ANALYSIS 
results<-LN3MVem(dat,pats,MaxIt,parm.conv,prob.conv,trt,parm.init)

if(method=="LN3MV*") post.probs<-LN3MVjoint(dat,results[[1]],results[[2]],results[[3]],pats)[[3]]  ### LN3MV Method

else{ ####LN3GV Method

# empirically integrate out vr
SUB.IT<-1000; jnt<-0;  est.parms<-results[[2]]; est.vg<-est.parms[1]; est.vt<-est.parms[2]; 
est.mu<-est.parms[3]; d0<-results[[4]][1]; s02<-results[[4]][2]

for(sub.it in 1:SUB.IT){
vr.sim<-d0*s02/qchisq(sub.it/(SUB.IT+1),d0)    ### Plug in quantiles from prior distribution for vr using estimated nu and phi. 
jnt<-jnt+LN3joint(dat,results[[1]],c(est.vg,est.vt,vr.sim,est.mu),patmat)[[1]] ### Using estimated parameters
  }

### Convert joint densities to posterior probs (and average over vr simulation iterations)
post.probs<-jnt/rowSums(jnt)

}

parms<-c(results[[1]],results[[2]],results[[4]]); names(parms)<-c(paste("pi",1:length(results[[1]]),sep=""),"gn.var","trt.var","mu","nu","phi")

}

### End of LN3MV/LN3GV analyses


if(method=="LN3"){

results<-LN3em(dat,pats,MaxIt,parm.conv,prob.conv,trt,parm.init)
post.probs<-LN3joint(dat,results[[1]],results[[2]],patmat)[[3]]
parms<-c(results[[1]],results[[2]]);names(parms)<-c(paste("pi",1:length(results[[1]]),sep=""),"gn.var","trt.var","err.var","mu")
}

if(method=="LNNMV*"|method=="LNNGV"){

##### LNNMV ANALYSIS 
results<-LNNMVem(dat,pats,MaxIt,parm.conv,prob.conv,trt,parm.init)

if(method=="LNNMV*") post.probs<-LNNMVjoint(dat,results[[1]],results[[2]],results[[3]],pats)[[3]]  ### LNNMV Method

else{ ####LNNGV Method

#empirically integrate out vr
SUB.IT<-1000; jnt<-0

est.parms<-results[[2]]; est.vt<-est.parms[1];  est.mu<-est.parms[2]
d0<-results[[4]][1];s02<-results[[4]][2]
for(sub.it in 1:SUB.IT){
vr.sim<-d0*s02/qchisq(sub.it/(SUB.IT+1),d0)    ### Using estimated parameters of vr dist'n
jnt<-jnt+LNNjoint(dat,results[[1]],c(est.vt,vr.sim,est.mu),patmat)[[1]] ### Using estimated parameters
  }
### Convert joint densities to posterior probs (and average over vr simulation iterations)
post.probs<-jnt/rowSums(jnt)

}
parms<-c(results[[1]],results[[2]],results[[4]]); names(parms)<-c(paste("pi",1:length(results[[1]]),sep=""),"trt.var","mu","nu","phi")
}

return(list("Pattern Probs"=post.probs,"Parameter Estimates"=parms))
}





