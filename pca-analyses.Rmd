
Early , Intermediate and Late Images of Transverse,saggetal and coronal data views

```{r}
  setwd("D:\\BRAIN-TUMOR-FDG")
  fi='fdg.dyn' ; fa='Attn'; faq='FDG.acqtimes90'
  tx  = matrix (scan(faq,skip=2),ncol=2,byrow = T)
  x=readBin(fa,n=128*128*35*31,numeric(),endian='big',size=4)
  x=array(x,c(128,128,35,31))
  z=readBin(fi,n=128*128*35*31,numeric(),endian='big',size=4)
  z=array(z,c(128,128,35,31))
  X = matrix(c(x),ncol = 31)
  gr=grey(c(1:128)/128)
  xs= apply(x,c(1,2,3),sum)
  n1=128;n2=128;n3=35
  #attenuation scan 
  xa = readBin(fa,n=128*128*35,numeric(),endian='big',size=4)
  xa=array(xa,c(128,128,35))
  xs=xa
  xe= apply(z[,,,2:9],c(1,2,3),sum) #early
  xi= apply(z[,,,15:19],c(1,2,3),sum) #intermediate
  xl= apply(z[,,,26:29],c(1,2,3),sum)  #late
  par(mfrow=c(3,3))
  #transverse
  image(xe[,rev(1:n2),n3/2],col=gr,axes=F,main="Transverse",
        xlab="Early")
  #coronal slice
  image(xe[,n2/2,],col=gr,axes=F,main="Coronal", xlab="Early")
  #saggatial slice
  image(xe[n1/2,rev(1:n2),],col=gr,axes=F,main="Sagittal",
        xlab="Early")
  #transverse
  image(xi[,rev(1:n2),n3/2],col=gr,axes=F,main="Transverse", 
        xlab="intermediate")
  #coronal slice
  image(xi[,n2/2,],col=gr,axes=F,main="Coronal", xlab="Intermediate")
  #saggatial slice
  image(xi[n1/2,rev(1:n2),],col=gr,axes=F,main="Sagittal",
        xlab="intermediate")
  #transverse
  image(xl[,rev(1:n2),n3/2],col=gr,axes=F,main="Transverse",
        xlab="late")
  #coronal slice
  image(xl[,n2/2,],col=gr,axes=F,main="Coronal",xlab="late")
  #saggatial slice
  image(xl[n1/2,rev(1:n2),],col=gr,axes=F,main="Sagittal",
        xlab="Late")
```

Scree Plot and Plots of Loadings of first 4 Principal Components

```{r}
  # Compute PCA of the result
  o=prcomp(X,scale=TRUE)
  #scree plot 
  lambda = o$sdev^2
  percentExp=100 * cumsum(lambda)/sum(lambda)   
  p=ncol(X)
  percentUnExp=100 - (100 * cumsum(lambda)/sum(lambda))
  matplot(1:p,cbind(percentExp,percentUnExp),main="scree",
          xlab="Number of Factors",ylab="% Variance",pch=16)
  legend(8,90, legend=c("variance unexplained", 
                        "variance explained"),col=c("red", "black"),
         text.font=1,pch=20:20, cex=1)
  #loading for 4 principle components
  matplot(tx[,1],o$rotation[,1:4],type="l",xlab="Time",ylab="",lty=5)
  legend('topright', legend=c("PCA 1","PCA 2","PCA 3","PCA 4"),
         col=c( "black","red","green","blue"), 
         text.font=1,lty=5, cex=0.9)
```

Images corresponding to Tranverse, Saggital and coronal views of first 4 PCAs

```{r}
  o=prcomp(X,scale=TRUE)
  z=array(c(o$x),c(128,128,35,31))  
  pca1=z[,,,1]
  pca2=z[,,,2]
  pca3=z[,,,3]
  par(mfrow=c(3,3))
  
  #Transverse
  image(pca1[,rev(1:n2),n3/2],col=gr,axes=F,
        main="Transverse",ylab="pca1")
  #coronal slice
  image(pca1[,n2/2,],col=gr,axes=F,
        main="coronal",ylab="pca1")
  #saggatial slice
  image(pca1[n1/2,rev(1:n2),],col=gr,axes=F,
        main="sagittal",ylab="pca1")
  #Transverse
  image(pca2[,rev(1:n2),n3/2],col=gr,axes=F,
        main="Transverse",ylab="pca2")
  #coronal slice
  image(pca2[,n2/2,],col=gr,axes=F,
        main="coronal",ylab="pca2")
  #saggatial slice
  image(pca2[n1/2,rev(1:n2),],col=gr,axes=F,
        main="sagittal",ylab="pca2")
  #Transverse
  image(pca3[,rev(1:n2),n3/2],col=gr,axes=F,
        main="Transverse",ylab="pca3")
  #coronal slice
  image(pca3[,n2/2,],col=gr,axes=F,
        main="coronal",ylab="pca3")
  #saggatial slice
  image(pca3[n1/2,rev(1:n2),],col=gr,axes=F,
        main="sagittal",ylab="pca3")
  
```


