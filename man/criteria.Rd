\name{criteria}
\alias{criteria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes two model selection criteria based
on bootstrapped regression estimates
%%  ~~function to do ... ~~
}
\description{Computes Shao's and Welsh's model selection criteria.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
criteria(beta, sigma, y, x, betas, cc = 4.685061, p, k = 1, rho.type = 2, tr = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{Vector of regression coefficient estimates.
%%     ~~Describe \code{beta} here~~
}
  \item{sigma}{Residual scale estimate.
%%     ~~Describe \code{sigma} here~~
}
  \item{y}{Vector of responses.
%%     ~~Describe \code{y} here~~
}
  \item{x}{\code{n} by \code{p} design matrix.
%%     ~~Describe \code{x} here~~
}
  \item{betas}{Matrix of bootstrapped regression coefficients (one in each row).
%%     ~~Describe \code{betas} here~~
}
  \item{cc}{Tuning constant for the \code{rho} function involved in the
  robust criterion.
%%     ~~Describe \code{cc} here~~
}
  \item{p}{Number of explanatory variables, used in the penalized criterion
  of Welsh.
%%     ~~Describe \code{p} here~~
}
  \item{k}{Penalty coefficient: \code{k * log(n) * p}, where \code{n} is the
  number of observations in the data.
%%     ~~Describe \code{k} here~~
}
  \item{rho.type}{An integer indicating which loss function is to be used:
    \code{1}: Tukey's bisquare, \code{2}: truncated squared loss, \code{3}: squared
    loss
%%     ~~Describe \code{rho.type} here~~
}
  \item{tr}{Truncation constant, only used if \code{rho.type = 2}
%%     ~~Describe \code{tr} here~~
}
}
\details{Computes Shao's and Welsh's model selection criteria.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with two components:
\item{s }{Shao's criterion}
\item{w }{Welsh's criterion}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Matias Salibian-Barrera
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
