

RFPE <- function(x,y,model,control=rlm.control(), sigma.full)
{
  if(missing(sigma.full))
    sigma.full <- my.rlm(x,y, control=rlm.control())$scale
  # Goodness of fit

  x0 <- x[,model, drop=FALSE]
  fit.mod <- my.rlm(x0,y, control=rlm.control())
  #print(fit.mod)
  beta.mod <- fit.mod$coef
  #print(beta.mod)
  resid <- y - x0 %*% beta.mod
  standres <- resid/sigma.full
  RFPE <- mean(Rho.MM(standres))

  # Add penalty term

  a <- mean(Psi.MM(standres)^2)
  b <- mean(Psiprime.MM(standres))
  q <- (d <- dim(x0))[2]
  n <- d[1]
  RFPE <- RFPE+q*a/(n*b)
  return(RFPE)
}
#RFPE(x,y,model)

Rho.MM <- function(x, cc=4.685061) {
  # Tukey's bisquare rho
  tmp <- x/cc
  tmp2 <- cc^2*tmp^2*(1-tmp^2+tmp^4/3)/2
  tmp2[ abs(tmp) > 1 ] <- cc^2/6
  return(tmp2)
}

Psi.MM <- function(x, cc=4.685061) {
  # Tukey's bisquare psi
  tmp <- x/cc
  tmp2 <- cc*tmp*(1-tmp^2)^2
  tmp2[ abs(tmp) > 1 ] <- 0
  return(tmp2)
}

Psiprime.MM <- function(x, cc=4.685061) {
  # Tukey's bisquare psi prime
  tmp <- x/cc
  tmp2 <-1-6*tmp^2+5*tmp^4
  tmp2[ abs(tmp) > 1 ] <- 0
  return(tmp2)
}


mn <- function(beta, sigma, y, x, cc = 4.685061, rho.type, tr ) {
  # for the data in (x,y) compute the sum of
  # Rho(resid[i] / sigma)
  r <- as.vector(y - x %*% beta) / sigma
  return(sum(Rho(r, cc=cc, type=rho.type)))

}

boot.mean <- function(betas, sigma, y, x, cc = 4.685061, rho.type,
                      tr)
{
  # compute \sum \rho((y-beta'x)/sigma) for
  # each beta in the rows of "betas"
  tmp <- apply(betas, 1, mn, sigma=sigma, y=y, x=x,
               cc=cc, rho.type=rho.type, tr=tr)
  # returns the average of these \sum\rho 's
  return(mean(tmp))
}

criteria <- function(beta, sigma, y, x, betas, cc=4.685061, p,
                     k=1, rho.type=2, tr=2)
{
  # compute the criteria
  n <- length(y)
  # Shao's bootstrap estimator
  # bootstrap mean of \sum \rho( r_i(beta) / sigma)
  # for each beta in the rows of "betas"
  a1 <- boot.mean(betas, sigma, y, x, cc, rho.type, tr)
  # Welsh's criteria = Shao's criteria + penalty + \sum(r_i(beta))
  # for the estimator in "beta"
  a2 <- mn(beta, sigma, y, x, cc, rho.type, tr) + k * log(n) * p + a1
  # a3 <- mn(beta, sigma, y, x, cc, rho.type, tr)

  #print(c(a1, a3))

  a1 <- a1 * sigma^2 / n
  a2 <- a2 * sigma^2 / n

  return(list(s=a1, w=a2))

}


Rho <- function(x, cc=4.685061, type=1, tr=2) {
  # type: 1 = Tukey's bisquare
  # 2 = truncated x^2
  # 3 = x^2 (LS)
  if(type == 3) return(x^2)
  else if(type == 2)
    return(pmin(x^2,tr^2))
  else {
    #tmp <- x/cc
    #tmp2 <- 3 * tmp^2 - 3 * tmp^4 + tmp^6
    #tmp2[abs(tmp) > 1] <- 1
    #return(tmp)
    tmp <- x/cc
    tmp2 <- cc^2*tmp^2*(1-tmp^2+tmp^4/3)/2
    tmp2[ abs(tmp) > 1 ] <- cc^2/6
    return(tmp2)
  }
}


"roboot" <-
  function(x, y, nboot, boot.samp, control=rlm.control())
  {
    M <- control$M
    Nres <- control$Nres
    seed <- control$seed
    fixed <- control$fixed
    tuning.chi <- control$tuning.chi
    tuning.psi <- control$tuning.psi
    calc.full <- control$calc.full
    max.it <- control$max.it
    groups <- control$groups
    n.group <- control$n.group
    k.fast.s <- control$k.fast.s
    n <- nrow(x)
    if(missing(nboot)) nboot <- n
    p <- ncol(x)
    a <- .C("R_rlm_rand",
            x = as.double(x),
            y = as.double(y),
            n = as.integer(n),
            p = as.integer(p),
            boot.samp = as.integer(boot.samp),
            Nres = as.integer(Nres),
            M = as.integer(M),
            nboot = as.integer(nboot),
            ours = double(M * p),
            full = double(M * p),
            coef = double(p),
            s = double(p),
            scale = as.double(0),
            seed = as.integer(seed),
            calc.full = as.integer(calc.full),
            tuning.chi = as.double(tuning.chi),
            tuning.psi = as.double(tuning.psi),
            as.integer(max.it), integer(1),
            as.integer(groups), as.integer(n.group),
            as.integer(k.fast.s))
    a$ours <- scale( matrix(a$ours, nrow = M),
                     center = -a$coef, scale=FALSE)
    a$full <- matrix(a$full, nrow = M)
    #full <- a$full[a$full[, 2] < 9000,  ]
    coef <- a$coef
    rank <- qr(x)$rank
    r1 <- 1:rank
    dn <- colnames(x)
    if (is.matrix(y)) {
      coef[-r1, ] <- NA
      dimnames(coef) <- list(dn, colnames(y))
    }
    else {
      coef[-r1] <- NA
      names(coef) <- dn
    }
    f <- x %*% as.matrix(a$coef)
    r <- as.matrix(y) - f
    z <- list(fitted.value=as.vector(f), residuals=as.vector(r),
              rank=rank, degree.freedom=n-rank, coefficients=coef, s=a$s,
              scale=a$scale, seed=a$seed, cov=var(a$ours,na.rm=TRUE),
              ours=a$ours, full=a$full)
    return(z)
  }



rlm.control <- function(M=2000,seed = 99, Nres = 500,
                        fixed = FALSE, tuning.chi = 1.54764, bb=0.5, tuning.psi = 4.685061,
                        groups = 5, n.group = 400, k.fast.s=1, max.it=50,
                        compute.rd = TRUE, calc.full=0)
{
  return(list(M=M, seed=seed, fixed=fixed, Nres = Nres,
              tuning.chi=tuning.chi, bb=bb, tuning.psi=tuning.psi,
              groups=groups, n.group=n.group,
              k.fast.s=k.fast.s, max.it=max.it,
              compute.rd=compute.rd, calc.full=calc.full))
}


"my.rlm" <-
  function(x, y, control=rlm.control())
  {
    x <- as.matrix(x)
    M <- control$M
    Nres <- control$Nres
    seed <- control$seed
    fixed <- control$fixed
    tuning.chi <- control$tuning.chi
    tuning.psi <- control$tuning.psi
    max.it <- control$max.it
    groups <- control$groups
    n.group <- control$n.group
    k.fast.s <- control$k.fast.s
    n <- nrow(x)
    p <- ncol(x)

    a <- .C('R_S_rlm', as.double(x),
            as.double(y), as.integer(n),
            as.integer(p),
            as.integer(Nres),
            as.integer(max.it),
            s=as.double(0),
            beta.s=as.double(rep(0,p)),
            beta.m=as.double(rep(0,p)),
            as.integer(1),
            as.integer(seed), as.double(tuning.chi),
            as.double(tuning.psi), as.integer(groups),
            as.integer(n.group), as.integer(k.fast.s))

    return(list(coef=a$beta.m, scale=a$s, coef.s = a$beta.s))

  }

