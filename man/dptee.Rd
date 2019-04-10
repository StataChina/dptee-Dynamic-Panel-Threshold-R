\name{dptee}
\alias{dptee}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model estimating
}
\description{
Function to estimate a Dynamic Panels with Threshold Effect and Endogeneity
}
\usage{
dptee(formula, data, con=0, trim=0.1, qn=20, ns=20, jj=20, effects, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{model formula,}
  \item{data}{data to be used,}
  \item{con}{constant in the upper regime. 1 to include, 0 otherwise,}
  \item{trim}{percentage to trim before search, single,}
  \item{qn}{number of quantiles to examine,}
  \item{ns}{iter num in averaging (1st uses analytic formula),}
  \item{jj}{iteration for linearity test p-value,}
  \item{effects}{individual or twoqays,}
  \item{...}{futher arguments.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Seo, Myung & Shin, Yongcheol. (2016). Dynamic panels with threshold effect and endogeneity. Journal of Econometrics. 195. 10.1016/j.jeconom.2016.03.005.
}
\author{
Frederico Uchoa
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{summary.dptee}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
