LN3MVjoint<-function(dat,prob,parm,vr,pats){
vg<-parm[1];vt<-parm[2];mu<-parm[3]

Gn<-apply(pats,1,max)

pmat<-matrix(.C("like",dat=as.numeric(dat-mu),dimd=as.integer(dim(dat)),
varr=as.numeric(vr),vart=as.numeric(vt),varg=as.numeric(vg),
as.integer(pats),as.integer(nrow(pats)),as.integer(Gn),
as.numeric(rep(0,nrow(dat)*nrow(pats))),PACKAGE="LN3GV")[[9]],nrow(dat),nrow(pats))


for(u in 1:nrow(pats)){
pmat[,u]=prob[u]*pmat[,u]+1e-200
}

h<-rowSums(pmat)  #### marginal densities for data
k<-pmat/h    ##### conditional densities for pattern given data (by gene)


return(list(pmat,h,k))
}

#### find expectation for E part of EM algorithm
LN3MVq<-function(parm,dat,pmat,prob,vr,pats){  

Q<-log(LN3MVjoint(dat,prob,parm,vr,pats)[[1]])*pmat
qq<--sum(Q)

return(qq)
}


LN3MVem<-function(dat,pats,MaxIt,parm.conv,prob.conv,trt,parm.init=NULL,prob.init=NULL){
neu<-length(trt)
if(is.null(prob.init)) pnew<-rep(1/nrow(pats),nrow(pats))
else pnew<-prob.init

tmn<-matrix(NA,nrow(dat),max(trt))
for(grp in 1:max(trt)){
tmn[,grp]<-apply(dat[,trt==grp],1,mean)
}

### Employ Smyth's method to compute gene specific residual variance estimates "vr"
MSE<-apply((dat-tmn[,trt])^2,1,sum)/(neu-max(trt))

z<-log(MSE)
mnz<-mean(z)

## solve for nu and phi
nuarg<-sum((z-mnz)^2)/(length(z)-1)-trigamma((neu-max(trt))/2)
dif<-function(x,y) abs(trigamma(x)-y)
inverse.trigamma<-function(y) optimize(dif,interval=c(0,100),y=y)$minimum
nu<-2*inverse.trigamma(nuarg)
phi<-exp(mnz-digamma((neu-max(trt))/2)+digamma(nu/2)- log(nu/(neu-max(trt))))

## compute vr
vr<- ((neu-max(trt))*MSE+nu*phi)/(neu-max(trt)+nu-2)

### initialize other model parameters

if(is.null(parm.init)) {
## compute Means
mu0<-mean(as.vector(as.matrix(dat)))


gmn<-apply(dat,1,mean)

xx<-apply(pats,1,table)
trt.var.coef<-0
for(i in 1:length(xx)){
trt.var.coef<-trt.var.coef+sum(xx[[i]]^2)*pnew[i]
}
trt.var.coef<-trt.var.coef/length(trt)^2

gn.var.wt<-c(.01,.1*(1:9),.99)

gn.var.try<-(var(gmn)-mean(vr)/length(trt))*gn.var.wt
trt.var.try<-(var(gmn)-mean(vr)/length(trt))*(1-gn.var.wt)/trt.var.coef

init.var<-gn.var.try
for(i in 1:length(init.var)){
init.var[i]<-sum(log(LN3MVjoint(dat,pnew,c(gn.var.try[i],trt.var.try[i],mu0),vr,pats)[[2]]+1e-200))}

i<-order(-init.var)[1]
vg0<-max(gn.var.try[i],.001); vt0<-max(trt.var.try[i],.001)
parm<-c(vg0,vt0,mu0); 
}
else parm<-parm.init

p<-0*pnew
parmnew<-parm

###### Set convergence criteria below ##########
it<-1

while((sum(abs(p-pnew)>prob.conv)+sum(abs((parm-parmnew)/parmnew)>parm.conv))>0&it<MaxIt){

#### Update pattern probabilities ####
p<-pnew

jnt<-LN3MVjoint(dat,p,parmnew,vr,pats)
pnew <- apply(jnt[[3]],2,mean)   ### average of conditional probabilities across genes

#### Update vg, vt and mu ####
parm<-parmnew
#print(paste("LN3GV vg, vt and mu:",parm))

optpar<- optim(parm,fn=LN3MVq,method = c("L-BFGS-B"),lower=c(.4*parm[-3],-2*abs(parm[3])),
upper=c(2.5*parm[-3],2*abs(parm[3])),dat=dat,pmat=jnt[[3]],prob=pnew,vr=vr,pats=pats)

parmnew<-optpar[[1]]; names(parmnew)<-c("gn.var","trt.var","mu")

names(pnew)<-paste("pi",1:length(pnew))
print(paste("Iteration Number:",it))
print(list("Pattern Weights"=pnew,"Model Parameters"=parmnew))
it<-it+1
}
return(list(pnew,parmnew,vr,c(nu,phi)))
}

