\name{rlm.control}
\alias{rlm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rlm.control(M = 2000, seed = 99, Nres = 500, fixed = FALSE, tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061, groups = 5, n.group = 400, k.fast.s = 1, max.it = 50, compute.rd = TRUE, calc.full = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{Nres}{
%%     ~~Describe \code{Nres} here~~
}
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{tuning.chi}{
%%     ~~Describe \code{tuning.chi} here~~
}
  \item{bb}{
%%     ~~Describe \code{bb} here~~
}
  \item{tuning.psi}{
%%     ~~Describe \code{tuning.psi} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{n.group}{
%%     ~~Describe \code{n.group} here~~
}
  \item{k.fast.s}{
%%     ~~Describe \code{k.fast.s} here~~
}
  \item{max.it}{
%%     ~~Describe \code{max.it} here~~
}
  \item{compute.rd}{
%%     ~~Describe \code{compute.rd} here~~
}
  \item{calc.full}{
%%     ~~Describe \code{calc.full} here~~
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
