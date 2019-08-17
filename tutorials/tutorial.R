library(PRIMAL)

## Dantzig selector
## Generate the design matrix and coefficient vector
n = 100 # sample number
d = 250 # sample dimension
c = 0.5 # correlation parameter
s = 20  # support size of coefficient
set.seed(1024)
X = scale(matrix(rnorm(n*d),n,d)+c*rnorm(n))/sqrt(n-1)*sqrt(n)
flag = runif(s,-1,1)
beta1<-c()
for(i in 1:s){
  if(flag[i]>=0) beta1[i]=rnorm(1,1,1)
  if(flag[i]<0)  beta1[i]=rnorm(1,-1,1)
}
beta = c(beta1, rep(0, d-s))

## Generate response using Gaussian noise, and solve the solution path
noise = rnorm(n)
Y = X%*%beta + noise

## Dantzig selection solved with parametric simplex method
fit.dantzig = Dantzig_solver(X, Y, max_it = 100, lambda_threshold = 0.01)

## lambdas used
print(fit.dantzig$lambda)

## number of nonzero coefficients for each lambda
print(fit.dantzig$df)

## coefficients and value for the i-th lambda
i = 30
print(fit.dantzig$lambda[i])
print(fit.dantzig$beta[,i])
print(fit.dantzig$value[i])

## Visualize the solution path
plot.primal(fit.dantzig)
plot.primal(fit.dantzig, n=1)


## Compressed Sensing
## Generate the design matrix and coefficient vector
n = 100 # sample number
d = 250 # sample dimension
c = 0.5 # correlation parameter
s = 20  # support size of coefficient
set.seed(1024)
X = scale(matrix(rnorm(n*d),n,d)+c*rnorm(n))/sqrt(n-1)*sqrt(n)
flag = runif(s,-1,1)
beta1<-c()
for(i in 1:s){
  if(flag[i]>=0) beta1[i]=rnorm(1,1,1)
  if(flag[i]<0)  beta1[i]=rnorm(1,-1,1)
}
beta = c(beta1, rep(0, d-s))

## Generate response using Gaussian noise, and solve the solution path
noise = rnorm(n)
Y = X%*%beta + noise

## Compressed Sensing solved with parametric simplex method
fit.compressed = CompressedSensing_solver(X, Y, max_it = 100, lambda_threshold = 0.01)

## lambdas used
print(fit.compressed$lambda)

## number of nonzero coefficients for each lambda
print(fit.compressed$df)

## coefficients and value for the i-th lambda
i = 30
print(fit.compressed$lambda[i])
print(fit.compressed$beta[,i])
print(fit.compressed$value[i])

## Visualize the solution path
plot.primal(fit.compressed)
plot.primal(fit.compressed, n=1)



## Quantile Regression
## Generate the design matrix and coefficient vector
n = 200 # sample number
d = 100 # sample dimension
c = 0.5 # correlation parameter
s = 20  # support size of coefficient
set.seed(1024)
X = scale(matrix(rnorm(n*d),n,d)+c*rnorm(n))/sqrt(n-1)*sqrt(n)
beta1<-rnorm(s)
beta = c(beta1, rep(0, d-s))

## Generate response using Gaussian noise, and solve the solution path
noise = rnorm(n)
Y = X%*%beta +noise

## Compressed Sensing solved with parametric simplex method
fit.quantile = QuantileRegression_solver(X, Y, max_it = 100, lambda_threshold = 0.01, tau = 0.3)

## lambdas used
print(fit.quantile$lambda)

## number of nonzero coefficients for each lambda
print(fit.quantile$df)

## coefficients and value for the i-th lambda
i = 30
print(fit.quantile$lambda[i])
print(fit.quantile$beta[,i])
print(fit.quantile$value[i])

## Visualize the solution path
plot.primal(fit.quantile)
plot.primal(fit.quantile, n=1)


## Sparse SVM
## Generate the design matrix and coefficient vector
n = 200 # sample number
d = 100 # sample dimension
c = 0.5 # correlation parameter
s = 20  # support size of coefficient
set.seed(1024)
X = matrix(rnorm(n*d),n,d)+c*rnorm(n)

## Generate response and solve the solution path
Y <- sample(c(-1,1),n,replace = T)

## Compressed Sensing solved with parametric simplex method
fit.SVM = SparseSVM_solver(X, Y, max_it = 1000, lambda_threshold = 0.01)

## lambdas used
print(fit.SVM$lambda)

## number of nonzero coefficients for each lambda
print(fit.SVM$df)

## coefficients and value for the i-th lambda
i = 500
print(fit.SVM$lambda[i])
print(fit.SVM$beta[,i])
print(fit.SVM$value[i])

## Visualize the solution path
plot.primal(fit.SVM)
plot.primal(fit.SVM, n=1)


