#include <math.h>

void like(double *d,int *dimd,double *varr,double *vart, double *varg,int *pattern,int *numpat,int *Gn,double *like)
{
    //d is gene expression data in vector form (by column, NOT row)
    //dimd = c(nrow(dat),ncol(dat))
    //varr is a vector of the gene specific residual variances
    //vart and varg are scalars for treatment and gene variance components, respectively
    //pattern is a vector containing grouping structures under patterns of interest; more naturally thought of as a matrix w\ #of rows=#of patterns and #of columns=#of E.U.'s
    //numpat is a scalar indicating the number of patterns of interest
    //Gn is a vector containing the number of groups for each pattern
    //like will store the multivariate normal density evaluated for each combination of gene and pattern

    int gene,pat,i,j,p,g,obs,npat=*numpat;
    double vt=*vart, vg=*varg;
    double det, inv[dimd[1]][dimd[1]],sum,pi=3.14159265;

    for (gene=0;gene<dimd[0];gene++){
        double vr=varr[gene];
        for (pat=0;pat<npat;pat++){

//obtain determinant and inverse matrix
// start with residual and gene effect
            det=pow(vr,dimd[1])*(1+dimd[1]*vg/vr);
            for (i=0;i<dimd[1];i++){
                for(j=0;j<dimd[1];j++){
                    if(i==j)
                    inv[i][j]=1/vr-vg/((dimd[1]*vg+vr)*vr);
                    else
                    inv[i][j]=-vg/((dimd[1]*vg+vr)*vr);
                }
            }
//Add in treatment effect
            double a[dimd[1]],iab[dimd[1]][dimd[1]],iabi[dimd[1]][dimd[1]];
            for(g=1;g<Gn[pat]+1;g++){
// GET vector a
                for (i=0;i<dimd[1];i++){
                    if (pattern[pat+i*npat]==g)
                    a[i]=sqrt(vt);
                    else
                    a[i]=0;
                }

//GET inv*ab

                for (i=0;i<dimd[1];i++){
                    for(j=0;j<dimd[1];j++){
                        sum=0;
                        for (p=0;p<dimd[1];p++){
                            sum=sum+inv[i][p]*a[p]*a[j]; //ab[i][j]=a[i]*a[j];
                        }
                        iab[i][j]=sum;
                    }
                }

//GET inv*ab*inv

                for (i=0;i<dimd[1];i++){
                    for(j=0;j<dimd[1];j++){
                        sum=0;
                        for (p=0;p<dimd[1];p++){
                            sum=sum+iab[i][p]*inv[p][j];
                        }
                        iabi[i][j]=sum;
                    }
                }

//Get b'*inv*a
                double bia=0;
                for (i=0;i<dimd[1];i++){
                    for(j=0;j<dimd[1];j++){
                        if(pattern[pat+i*npat]==g&&pattern[pat+j*npat]==g)
                        bia=bia+inv[i][j];
                    }
                }
                bia=bia*vt;

//Update det
                det=det*(1+bia);

//Updata inv

                for (i=0;i<dimd[1];i++){
                    for(j=0;j<dimd[1];j++){
                        inv[i][j]= inv[i][j]-iabi[i][j]/(1+bia);
                    }
                }
            }

//Obtain product for exponent
            double prd=0;
            for (obs=0; obs<dimd[1]; obs++){
                sum=0;
                for (i=0; i<dimd[1]; i++){
                    sum=sum+d[gene+i*dimd[0]]*inv[i][obs];
                }
                sum=sum*d[gene+obs*dimd[0]];

                prd=prd+sum;
            }
//Evaluate likelihood
            like[gene+pat*dimd[0]]=exp(-.5*prd)/(pow(sqrt(2*pi),dimd[1])*sqrt(det));
        }
    }
}


