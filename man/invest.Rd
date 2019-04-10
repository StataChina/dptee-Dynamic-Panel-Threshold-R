\name{invest}
\alias{invest}
\docType{data}
\title{
Data set used in Seo and Shin (2016)
}
\description{
A balanced panel of 565 US firms over the period 1973-1987, which is extracted from an original data set constructed by Hall and Hall (1993) and used in Seo & Shin (2016).
}
\usage{data("invest")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{n}}{firm index}
    \item{\code{t}}{time index}
    \item{\code{y}}{investment/assets}
    \item{\code{c}}{cash-flow/assets}
    \item{\code{Tq}}{Tobin's Q}
    \item{\code{d}}{debt/assets}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
https://drive.google.com/open?id=0B8V9YA1XnGExdGFjay1uQlgyUHM
}
\references{
Seo, Myung & Shin, Yongcheol. (2016). Dynamic panels with threshold effect and endogeneity. Journal of Econometrics. 195.
}
\examples{
data(invest)
## maybe str(invest) ; plot(invest) ...
}
\keyword{datasets}
