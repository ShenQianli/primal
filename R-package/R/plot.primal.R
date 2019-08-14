plot.primal<-function(str, n=NULL){
  tt<-str$iterN
  if(is.null(n)){
    opar<-par(no.readonly=TRUE)
    par(mfrow=c(1,3),family="serif")
    matplot(str$lambda[1:tt],t(str$beta[,1:tt]),type="l",main="Regularization Path",
          xlab="Regularization Parameter", ylab="Coefficient",cex.main=2,cex.lab=1.6)
    matplot(1:tt, t(str$beta),type="l",main="Regularization Path",
            xlab="Iteration", ylab="Coefficient",cex.main=2,cex.lab=1.6)
    plot(1:tt, str$lambda[1:tt], type="l", main="Value of Lambda along the Path",
         xlab="Iteration", ylab="Lambda",cex.main=2,cex.lab=1.6)
    par(opar)
  }else{
    opar<-par(no.readonly=TRUE)
  par(family="serif")
  switch(n,
         matplot(str$lambda[1:tt],t(str$beta[,1:tt]),type="l",main="Regularization Path",
                   xlab="Regularization Parameter", ylab="Coefficient",cex.main=2,cex.lab=1.6),
         matplot(1:tt, t(str$beta[,1:tt]),type="l",main="Regularization Path",
                 xlab="Iteration", ylab="Coefficient",cex.main=2,cex.lab=1.6),
         plot(1:tt, str$lambda[1:tt], type="l", main="Value of Lambda along the Path",
              xlab="Iteration", ylab="Lambda",cex.main=2,cex.lab=1.6))
  par(opar)

  }
}


print.psm<-function(str,...){
  cat("\n Parametric Simplex Method solving ")
  cat(str$type,"  problem\n")
  cat("iteration times = ",str$t)
  cat("lambda list:\n")
  print(signif(str$lambda,digits=3))
  cat("Degree of freedom:",min(str$df),"----->",max(str$df),"\n")
  if (units.difftime(str$runtime)=="secs") unit="secs"
  if (units.difftime(str$runtime)=="mins") unit="mins"
  if (units.difftime(str$runtime)=="hours") unit="hours"
  cat("Runtime:",str$runtime," ",unit,"\n")
}
