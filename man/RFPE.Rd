\name{RFPE}
\alias{RFPE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust future prediction error
%%  ~~function to do ... ~~
}
\description{This function computes the 
robust future prediction error as
described in Section 5.12 of Maronna et al 2006.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RFPE(x, y, model, control = rlm.control(), sigma.full)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{full design matrix, matrix containing all explanatory variables
%%     ~~Describe \code{x} here~~
}
  \item{y}{vector of responses
%%     ~~Describe \code{y} here~~
}
  \item{model}{vector of integers corresponding to the columns of \code{x} to be
  included
%%     ~~Describe \code{model} here~~
}
  \item{control}{a list of control options as returned by \code{rlm.control}
%%     ~~Describe \code{control} here~~
}
  \item{sigma.full}{estimated residual scale, generally computed using all
  available features (explanatory variables)
%%     ~~Describe \code{sigma.full} here~~
}
}
\details{Robust future prediction error
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a real value
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Maronna, R., Martin, D., Yohai, V.J. (2006). Robust Statistics, theory
and methods. Wiley and Sons.
%% ~put references to the literature/web site here ~
}
\author{Matias Salibian-Barrera
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{rlm.control}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, model, control = rlm.control(), sigma.full)
{
    if (missing(sigma.full))
        sigma.full <- my.rlm(x, y, control = rlm.control())$scale
    x0 <- x[, model, drop = FALSE]
    fit.mod <- my.rlm(x0, y, control = rlm.control())
    beta.mod <- fit.mod$coef
    resid <- y - x0 \%*\% beta.mod
    standres <- resid/sigma.full
    RFPE <- mean(Rho.MM(standres))
    a <- mean(Psi.MM(standres)^2)
    b <- mean(Psiprime.MM(standres))
    q <- (d <- dim(x0))[2]
    n <- d[1]
    RFPE <- RFPE + q * a/(n * b)
    return(RFPE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
