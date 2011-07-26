
//send

void product(double *d,int *dims,double *inv,double *prod)
{
    //dims=c(nrow(dat),ncol(dat),#of patterns)
    //inv is a 3-dim array where (in R) last index provides inverse matrices... e.g. (in R) inv[,,1] is inverse of covariance matrix under pattern 1
double m[dims[1]],s;
int row,col,pat,i;

for (pat=0; pat<dims[2]; pat++){

for (row=0; row<dims[0]; row++){
s=0;
for (col=0; col<dims[1]; col++){
m[col]=0;
for (i=0; i<dims[1]; i++){
m[col]=m[col]+d[row+i*dims[0]]*inv[dims[1]*dims[1]*pat+dims[1]*i+col];
}
m[col]=m[col]*d[row+col*dims[0]];
s=s+m[col];
}
prod[row+dims[0]*pat]=s;
}
}

}

//(from R)inv[c,b,a]=inv[a][b][c]=inv[dims[1]*dims[1]*a+dims[1]*b+c]
//(from R)dat[a,b]=d[a+b*dims[0]]
