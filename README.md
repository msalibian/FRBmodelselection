Model Selection with the Fast and Robust Bootstrap
================
Matias Salibian
2016-08-03

Robust model selection with the FRB
-----------------------------------

This package implements a robust model selection procedure for linear regression models based on MM-estimators and the [Fast and Robust Bootstrap](http://dx.doi.org/10.1214/aos/1021379865) as described in [Salibian-Barrera, M. and Van Aelst, S. (2008)](http://dx.doi.org/10.1016/j.csda.2008.05.007).

To install it use the following commands (assuming that you have the `devtools` package from [CRAN](https://cran.r-project.org) already installed):

``` r
library(devtools)
install_github("msalibian/FRBmodelselection")
```

The following script illustrates this method when applied to the well-known Boston Housing data set.

``` r
library(FRBmodelselection)
data(Boston, package='MASS')
xx <- model.matrix(medv ~ ., data=Boston )
y <- as.vector(Boston$medv)
```

The "classical" stepwise methods all yield the same optimal submodel with 11 features:

``` r
library(MASS)
a.lm <- lm(medv ~ ., data = Boston)
a.aic <- stepAIC(a.lm, direction='both', k=2, trace=0)
n <- length(y)
a.bic <- stepAIC(a.lm, direction='both', k=log(n), trace=0)
a.cp <- stepAIC(a.lm, direction='both', k=2, scale=summary(a.lm)$sigma, trace=0)
dimnames(a.aic$model)[[2]]
```

    ##  [1] "medv"    "crim"    "zn"      "chas"    "nox"     "rm"      "dis"    
    ##  [8] "rad"     "tax"     "ptratio" "black"   "lstat"

``` r
dimnames(a.bic$model)[[2]]
```

    ##  [1] "medv"    "crim"    "zn"      "chas"    "nox"     "rm"      "dis"    
    ##  [8] "rad"     "tax"     "ptratio" "black"   "lstat"

``` r
dimnames(a.cp$model)[[2]]
```

    ##  [1] "medv"    "crim"    "zn"      "chas"    "nox"     "rm"      "dis"    
    ##  [8] "rad"     "tax"     "ptratio" "black"   "lstat"

Using a backwards stepwise approach with the Robust Future Prediction Error criterion yields the following model, also with 11 variables:

``` r
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
(model.rfpe <- models[[ which.min(criterias.rfpe) ]])
```

    ##  [1]  1  2  3  5  6  7  8  9 10 11 12 13

This model is slighlty different from the one found by the AIC-based approaches:

``` r
dimnames(xx)[[2]][model.rfpe]
```

    ##  [1] "(Intercept)" "crim"        "zn"          "chas"        "nox"        
    ##  [6] "rm"          "age"         "dis"         "rad"         "tax"        
    ## [11] "ptratio"     "black"

``` r
dimnames(a.aic$model)[[2]]
```

    ##  [1] "medv"    "crim"    "zn"      "chas"    "nox"     "rm"      "dis"    
    ##  [8] "rad"     "tax"     "ptratio" "black"   "lstat"

Now we use backward stepwise using Shao's and Welsh's criterion combined with the fast and robust boostrap. We use `b = 1000` bootstrap samples of size `nboot = 150`.

``` r
b <- 1000; nboot <- 150; p <- 14
```

We generate the bootstrap samples, and set up lists to store the optimal models to be found for each size 1, 2, ..., p.

``` r
set.seed(123)
boot.samp <- matrix(sample(n, b*nboot, repl=TRUE), b, nboot)
model.s <- model.w <- 1:14
models.s <- models.w <- vector('list', length(model.s))
```

We will use a loss function in Tukey's bi-square family, and we do not want to compute the fully-bootstrapped estimators:

``` r
k <- 1 ; rho.type <- 1 ; tr <- 2
my.ctrl <- rlm.control(M=b, calc.full=0)
criterias.s <- criterias.w <- rep(NA,14)
```

Now we fit the full model, and then proceed in a stewpwise manner:

``` r
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
model.shao <- models.s[[ which.min(criterias.s) ]]
model.welsh <- models.w[[ which.min(criterias.w) ]]
model.rfpe
```

    ##  [1]  1  2  3  5  6  7  8  9 10 11 12 13

The models selected with these methods are slighlty different, but noticeably smaller than those found using AIC/BIC or Mallow's Cp:

``` r
dimnames(xx)[[2]][model.shao]
```

    ## [1] "(Intercept)" "rm"          "age"         "dis"         "tax"        
    ## [6] "ptratio"     "black"

``` r
dimnames(xx)[[2]][model.welsh]
```

    ## [1] "(Intercept)" "crim"        "rm"          "age"         "dis"        
    ## [6] "ptratio"     "black"

``` r
dimnames(xx)[[2]][model.rfpe]
```

    ##  [1] "(Intercept)" "crim"        "zn"          "chas"        "nox"        
    ##  [6] "rm"          "age"         "dis"         "rad"         "tax"        
    ## [11] "ptratio"     "black"
