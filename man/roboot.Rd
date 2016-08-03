\name{roboot}
\alias{roboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implements the Fast and Robust Bootstrap
%%  ~~function to do ... ~~
}
\description{Implements the Fast and Robust Bootstrap
for MM-regression estimators as proposed in
Salibian-Barrera and Zamar (2002)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
roboot(x, y, nboot, boot.samp, control = rlm.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{nboot}{
%%     ~~Describe \code{nboot} here~~
}
  \item{boot.samp}{
%%     ~~Describe \code{boot.samp} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Salibian-Barrera, Matias; Zamar, Ruben H. Bootrapping robust estimates of regression. Ann. Statist. 30 (2002), no. 2, 556--582. doi:10.1214/aos/1021379865. http://projecteuclid.org/euclid.aos/1021379865.
%% ~put references to the literature/web site here ~
}
\author{Matias Salibian-Barrera
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{my.rlm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, nboot, boot.samp, control = rlm.control()) 
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
    if (missing(nboot)) 
        nboot <- n
    p <- ncol(x)
    a <- .C("R_rlm_rand", x = as.double(x), y = as.double(y), 
        n = as.integer(n), p = as.integer(p), boot.samp = as.integer(boot.samp), 
        Nres = as.integer(Nres), M = as.integer(M), nboot = as.integer(nboot), 
        ours = double(M * p), full = double(M * p), coef = double(p), 
        s = double(p), scale = as.double(0), seed = as.integer(seed), 
        calc.full = as.integer(calc.full), tuning.chi = as.double(tuning.chi), 
        tuning.psi = as.double(tuning.psi), as.integer(max.it), 
        integer(1), as.integer(groups), as.integer(n.group), 
        as.integer(k.fast.s))
    a$ours <- scale(matrix(a$ours, nrow = M), center = -a$coef, 
        scale = FALSE)
    a$full <- matrix(a$full, nrow = M)
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
    f <- x \%*\% as.matrix(a$coef)
    r <- as.matrix(y) - f
    z <- list(fitted.value = as.vector(f), residuals = as.vector(r), 
        rank = rank, degree.freedom = n - rank, coefficients = coef, 
        s = a$s, scale = a$scale, seed = a$seed, cov = var(a$ours, 
            na.rm = TRUE), ours = a$ours, full = a$full)
    return(z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
