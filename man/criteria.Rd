\name{criteria}
\alias{criteria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
criteria(beta, sigma, y, x, betas, cc = 4.685061, p, k = 1, rho.type = 2, tr = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{betas}{
%%     ~~Describe \code{betas} here~~
}
  \item{cc}{
%%     ~~Describe \code{cc} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{rho.type}{
%%     ~~Describe \code{rho.type} here~~
}
  \item{tr}{
%%     ~~Describe \code{tr} here~~
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
function (beta, sigma, y, x, betas, cc = 4.685061, p, k = 1, 
    rho.type = 2, tr = 2) 
{
    n <- length(y)
    a1 <- boot.mean(betas, sigma, y, x, cc, rho.type, tr)
    a2 <- mn(beta, sigma, y, x, cc, rho.type, tr) + k * log(n) * 
        p + a1
    a3 <- mn(beta, sigma, y, x, cc, rho.type, tr)
    a1 <- a1 * sigma^2/n
    a2 <- a2 * sigma^2/n
    return(list(s = a1, w = a2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
