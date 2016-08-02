\name{RFPE}
\alias{RFPE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RFPE(x, y, model, control = rlm.control(), sigma.full)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{sigma.full}{
%%     ~~Describe \code{sigma.full} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
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
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
