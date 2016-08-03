\name{rlm.control}
\alias{rlm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Options for the MM-regression estimator and the fast
and robust bootstrap.
%%  ~~function to do ... ~~
}
\description{This function returns a list of 
tuning options for \code{\link{my.rlm}} and \code{\link{roboot}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rlm.control(M = 2000, seed = 99, Nres = 500, fixed = FALSE, 
tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061, groups = 5, 
n.group = 400, k.fast.s = 1, max.it = 50, compute.rd = TRUE, calc.full = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{Number of bootstrap samples to use.
%%     ~~Describe \code{M} here~~
}
  \item{seed}{Seed used for the pseudo-random numbers generator used to 
compute the MM-regression estimates.
%%     ~~Describe \code{seed} here~~
}
  \item{Nres}{Number of sub-samples used to compute the S-regression estimator. 
%%     ~~Describe \code{Nres} here~~
}
  \item{fixed}{Unused.
%%     ~~Describe \code{fixed} here~~
}
  \item{tuning.chi}{Tuning constant for the M-scale estimator that is
  minimized to compute the S-regression estimator.
%%     ~~Describe \code{tuning.chi} here~~
}
  \item{bb}{Consistency parameter for the M-scale estimator. It should depend
  on the value of \code{tuning.chi} above.
%%     ~~Describe \code{bb} here~~
}
  \item{tuning.psi}{Tuning constant for the MM-regression estimator.
%%     ~~Describe \code{tuning.psi} here~~
}
  \item{groups}{Parameters of the Fast-S algorithm. 
%%     ~~Describe \code{groups} here~~
}
  \item{n.group}{Parameters of the Fast-S algorithm. 
%%     ~~Describe \code{n.group} here~~
}
  \item{k.fast.s}{Parameters of the Fast-S algorithm. 
%%     ~~Describe \code{k.fast.s} here~~
}
  \item{max.it}{Parameters of the Fast-S algorithm. 
%%     ~~Describe \code{max.it} here~~
}
  \item{compute.rd}{Unused.
%%     ~~Describe \code{compute.rd} here~~
}
  \item{calc.full}{logical value indicating whether
  \code{\link{roboot}} should also compute fully bootstrapped
  MM-regression estimators. 
%%     ~~Describe \code{calc.full} here~~
}
}
\details{Returns a list of all tuning constants and their values, which
are used by \code{my.rlm}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a list of all its arguments with their corresponding values.
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

\seealso{\code{\link{my.rlm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (M = 2000, seed = 99, Nres = 500, fixed = FALSE, tuning.chi = 1.54764, 
    bb = 0.5, tuning.psi = 4.685061, groups = 5, n.group = 400, 
    k.fast.s = 1, max.it = 50, compute.rd = TRUE, calc.full = 0) 
{
    return(list(M = M, seed = seed, fixed = fixed, Nres = Nres, 
        tuning.chi = tuning.chi, bb = bb, tuning.psi = tuning.psi, 
        groups = groups, n.group = n.group, k.fast.s = k.fast.s, 
        max.it = max.it, compute.rd = compute.rd, calc.full = calc.full))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
