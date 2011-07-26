LN3joint<-function(dat,prob,parm,patmat){
vg<-parm[1];vt<-parm[2];vr<-parm[3];mu<-parm[4]
neu<-nrow(patmat[[1]])
npat<-length(patmat)
pmat<-matrix(NA,nrow(dat),length(prob))
inv<-array(0,c(neu,neu,npat))
dt<-rep(0,npat)
for(u in 1:npat){
	cmat<-matrix(vg,neu,neu)+vt*patmat[[u]]+diag(rep(vr,neu))
	dt[u]<-abs(det(cmat))
	inv[,,u]<-solve(cmat)
	}
prod<-matrix(.C("product",d=as.numeric(dat-mu),dims=as.integer(c(dim(dat),
npat)),inv=as.numeric(inv),
prod=as.numeric(rep(0,nrow(dat)*npat)),PACKAGE="LN3GV")[[4]],nrow(dat),npat)

prod<--.5*t(t(prod)+neu*log(2*pi)+log(dt))
prod<-t(t(prod)+log(prob))
pmat<-exp(prod)+1e-200

h<-rowSums(pmat)  #### marginal densities for data
k<-pmat/h    ##### conditional densities for pattern given data (by gene)
return(list(pmat,h,k))
}


#### find expectation for E part of EM algorithm
LN3q<-function(parm,dat,pmat,prob,patmat){  
Q<-log(LN3joint(dat,prob,parm,patmat)[[1]])*pmat
qq<--sum(Q)
return(qq)
}


LN3em<-function(dat,pats,MaxIt,parm.conv,prob.conv,trt,parm.init=NULL,prob.init=NULL){

patmat<-vector("list",nrow(pats))
for(i in 1:nrow(pats)){
patmat[[i]]<-matrix(0,ncol(pats),ncol(pats))
for(j in 1:ncol(pats)){
for(k in 1:j){
if(pats[i,j]==pats[i,k]){
patmat[[i]][j,k]<-1
patmat[[i]][k,j]<-1
} } } }

if(is.null(prob.init)) pnew<-rep(1/nrow(pats),nrow(pats))
else pnew<-prob.init

if(is.null(parm.init)){
mu0<-mean(as.vector(dat))

tmn<-matrix(NA,nrow(dat),max(trt))
for(grp in 1:max(trt)){
tmn[,grp]<-apply(dat[,trt==grp],1,mean) }

## residual variance component estimate
vr0<-mean(as.vector((dat-tmn[,trt])^2))*length(trt)/(length(trt)-max(trt))

gmn<-apply(dat,1,mean)

xx<-apply(pats,1,table)
trt.var.coef<-0
for(i in 1:length(xx)){
trt.var.coef<-trt.var.coef+sum(xx[[i]]^2)*pnew[i]
}
trt.var.coef<-trt.var.coef/length(trt)^2

gn.var.wt<-c(.01,.1*(1:9),.99)

gn.var.try<-(var(gmn)-vr0/length(trt))*gn.var.wt
trt.var.try<-(var(gmn)-vr0/length(trt))*(1-gn.var.wt)/trt.var.coef

init.var<-gn.var.try
for(i in 1:length(init.var)){
init.var[i]<-sum(log(LN3joint(dat,pnew,c(gn.var.try[i],trt.var.try[i],vr0,mu0),patmat)[[2]]+1e-200))}

i<-order(-init.var)[1]
vg0<-max(gn.var.try[i],.001); vt0<-max(trt.var.try[i],.001)
parm<-c(vg0,vt0,vr0,mu0)
}

else parm<-parm.init
 
p<-0*pnew
parmnew<-parm

it<-1
###### Set convergence criteria below ##########

while((sum(abs(p-pnew)>prob.conv)+sum(abs((parm-parmnew)/parmnew)>parm.conv))>0&it<MaxIt){

#### Update pattern probabilities ####
p<-pnew
jnt<-LN3joint(dat,p,parmnew,patmat)

prob<-rowSums(is.na(jnt[[3]])+is.nan(jnt[[3]]))

pnew <- apply(jnt[[3]][prob==0,],2,mean)   ### average of conditional probabilities across genes

#### Update vg, vt, vr and mu ####
parm<-parmnew; 
#print(paste("LN3 vg, vt, vr and mu:",parm))

optpar<- optim(parm,fn=LN3q,method = c("L-BFGS-B"),lower=c(.3*parm[-4],
-2*abs(parm[4])),upper=c(2*parm[-4],2*abs(parm[4])),dat=dat,pmat=jnt[[3]],prob=pnew,patmat=patmat)

parmnew<-optpar[[1]]; names(parmnew)<-c("gn.var","trt.var","err.var","mu")

names(pnew)<-paste("pi",1:length(pnew))
print(paste("Iteration Number:",it))
print(list("Pattern Weights"=pnew,"Model Parameters"=parmnew))

it<-it+1
}

return(list(patprob=pnew,parameters=parmnew))
}


