LNNjoint<-function(dat,prob,parm,patmat){
vt<-parm[1];vr<-parm[2];mu<-parm[3]
neu<-nrow(patmat[[1]])
npat<-length(patmat)

pmat<-matrix(NA,nrow(dat),length(prob))
inv<-array(0,c(neu,neu,npat))
dt<-rep(0,npat)
for(u in 1:npat){
	cmat<-vt*patmat[[u]]+diag(rep(vr,neu))
	dt[u]<-abs(det(cmat))
	inv[,,u]<-solve(cmat)
	}

prod<-matrix(.C("product",d=as.numeric(dat-mu),dims=as.integer(c(dim(dat),npat)),inv=as.numeric(inv),
prod=as.numeric(rep(0,nrow(dat)*npat)),PACKAGE="LN3GV")[[4]],nrow(dat),npat)

prod<--.5*t(t(prod)+neu*log(2*pi)+log(dt))
prod<-t(t(prod)+log(prob))
pmat<-exp(prod)+1e-200


h<-rowSums(pmat)  #### marginal densities for data
k<-pmat/h    ##### conditional densities for pattern given data (by gene)

return(list(pmat,h,k))

}



