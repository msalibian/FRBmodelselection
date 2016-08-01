Model Selection with the Fast and Robust Bootstrap
================
Matias Salibian
2016-08-01

Robust model selection with the FRB
-----------------------------------

This repository contins `R` code implementing a robust model selection procedure for linear regression models based on MM-estimators and the [Fast and Robust Bootstrap](http://dx.doi.org/10.1214/aos/1021379865) as described in [Salibian-Barrera, M. and Van Aelst, S. (2008)](http://dx.doi.org/10.1016/j.csda.2008.05.007).

You will need to create a dynamic library from the code in `FRB-model-selection.c` using, for example, the following command in your shell:

``` r
R CMD SHLIB FRB-model-selection.c
```

Note that if you are running Windows, you will need to have installed the [RTools package](https://cran.r-project.org/bin/windows/Rtools/) from CRAN.

The following script illustrates this method when applied to the well-known Boston Housing data set.

``` r
library(robustbase)
source('FRB-model-selection-functions.R')

# this needs to be changed, when I build a proper R package
if(length(grep('linux', version$os)) > 0) dyn.load('FRB-model-selection.so')
if(length(grep('windows', version$os)) > 0) dyn.load('FRB-model-selection.dll')

# setup a simple example

data(Boston, package='MASS')
xx <- model.matrix(medv ~ ., data=Boston )
y <- as.vector(Boston$medv)

n <- length(y)

# Model selection using stepwise + AIC - BIC - Cp
library(MASS)
a.lm <- lm(medv ~ ., data = Boston)
a.aic <- stepAIC(a.lm, direction='both', k=2, trace=0)
a.bic <- stepAIC(a.lm, direction='both', k=log(n), trace=0)
a.cp <- stepAIC(a.lm, direction='both', k=2, scale=summary(a.lm)$sigma, trace=0)


# Model selection using backward stepwise + RFPE
p <- 14; model <- 1:14 ;
criterias.rfpe <- rep(NA, length(model))
sigma <- my.rlm(xx,y)$scale # full model scale
models <- vector('list', length(model))
models[[length(models)]] <- model
criterias.rfpe[length(models)] <- RFPE(x=xx, y=y, model=model, sigma.full=sigma)
while( p > 1 ) {
  best.c <- 10^10
  for(i in 2:p) {
    model.c <- model[-i]
    a <- RFPE(x=xx, y=y, model=model.c, sigma.full=sigma)
    if( a < best.c ) {
      new.model <- model.c
      best.c <- a
      criterias.rfpe[p - 1] <- a
    }
  }
  model <- models[[ p - 1]] <- new.model
  p <- p - 1
}
model.rfpe <- models[[ which.min(criterias.rfpe) ]]
print(c(min(criterias.rfpe),
  RFPE(x=xx, y=y, model=model.rfpe, sigma.full=sigma) ))



# Model selection using backward stepwise + Shao's and Welsh's criterion + FRB
b <- 1000; nboot <- 150; p <- 14; set.seed(123)
# Generate bootstrap samples
boot.samp <- matrix(sample(n, b*nboot, repl=TRUE), b, nboot)
model.s <- model.w <- 1:14
models.s <- models.w <- vector('list', length(model.s))
# Parameters for the criterion function
k <- 1 ; rho.type <- 1 ; tr <- 2
# control options for the MM estimator
my.ctrl <- rlm.control(M=b, calc.full=0)

criterias.s <- criterias.w <- rep(NA,14)
x0 <- xx[, model.s] # same as model.w -- should be full model
tmp <- roboot(x=x0, y=y, nboot=nboot, boot.samp=boot.samp,
              control=my.ctrl)
beta <- tmp$coef
betas <- tmp$ours
sigma <- tmp$scale
a <- criteria(beta=beta, sigma=sigma,
              y=y, x=x0, betas=betas, p=p, k=k, rho.type=rho.type, tr=tr)

criterias.w[length(models.w)] <- a$w
criterias.s[length(models.s)] <- a$s

while( p > 1 ) {
  # Welsh's criterion
  best.c <- 10^10
  for(i in 2:p) {
    model.c <- model.w[-i]
    x0 <- xx[, model.c, drop=FALSE]
    tmp <- roboot(x=x0, y=y, nboot=nboot, boot.samp=boot.samp,
                  control=my.ctrl)
    beta <- tmp$coef
    betas <- tmp$ours
    a <- criteria(beta=beta, sigma=sigma,
                  y=y, x=x0, betas=betas, p=p, k=k, rho.type=rho.type, tr=tr)
    if( a$w < best.c ) {
      new.model <- model.c
      best.c <- a$w
      criterias.w[p-1] <- a$w
    }
  }
  models.w[[ p - 1 ]] <- model.w <- new.model

  # Shao's criterion
  best.c <- 10^10
  for(i in 2:p) {
    model.c <- model.s[-i]
    x0 <- xx[, model.c, drop=FALSE]
    tmp <- roboot(x=x0, y=y, nboot=nboot, boot.samp=boot.samp,
                  control=my.ctrl)
    beta <- tmp$coef
    betas <- tmp$ours
    a <- criteria(beta=beta, sigma=sigma,
                  y=y, x=x0, betas=betas, p=p, k=k, rho.type=rho.type, tr=tr)
    if( a$s < best.c ) {
      new.model <- model.c
      best.c <- a$s
      criterias.s[p-1] <- a$s
    }
  }
  models.s[[ p - 1 ]] <- model.s <- new.model
  p <- p - 1
}


( model.shao <- models.s[[ which.min(criterias.s) ]] )
( model.welsh <- models.w[[ which.min(criterias.w) ]] )
( model.rfpe )
```
