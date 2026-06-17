\name{summary.ctimeROC}
\alias{summary.ctimeROC}
\title{
Summary method for \code{ctimeROC} objects
}
\description{
Default summary method for objects fitted with \code{ctimeROC.np()}, \code{ctimeROC.sp()}, or \code{ctimeROC.ps()} functions.
}
\usage{
\method{summary}{ctimeROC}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{ctimeROC} as produced by \code{ctimeROC.np()}, \code{ctimeROC.sp()}, or \code{ctimeROC.ps()}.}
  \item{\dots}{Further arguments passed to or from other methods. Not yet implemented.}
}
\details{
The information printed depends on the method. In all cases, the call to the function, the method, and the sample sizes are printed.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

\seealso{
\code{\link{ctimeROC.np}}, \code{\link{ctimeROC.sp}} or \code{\link{ctimeROC.ps}}.
}
\examples{
  library(CondTimeROC)
  # Simulate the data
  set.seed(123)
  # Sample size
  n <- 600
  # Covariate
  x <- rnorm(n, 1, 1)
  # Biomarker
  y <- rnorm(n, x, 1)
  # Survival time
  u <- runif(n)
  t <- exp(log(-log(1-runif(n))) - y - 0.5*x)/0.1
  # Censoring time
  c <- pmin(rexp(n, rate = 1/(1+4*abs(x))), 20)
  # Censoring indicator
  delta <- as.numeric(t <= c)
  # Observed time
  t <- pmin(t,c)
  # Data
  data <- data.frame(x = x, y = y, t = t, delta = delta)
  
  # Penalised-based approach
  \donttest{
  croc.ps <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ s(y) + s(x) + s(t)", 
    formula.biomarker = list(mean = y ~ s(x), var = ~ s(x)), data = data,
    cutoffs = seq(-3, 5, l = 200), 
    predict.time = 1,
    cov.values = data.frame(x = 0))
  
  summary(croc.ps)
  }
  \dontshow{
  croc.ps <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ y + x + t", 
    formula.biomarker = list(mean = y ~ x, var = ~ x), data = data, 
    predict.time = 1, 
    cov.values = data.frame(x = 0), 
    cutoffs = seq(-3, 5, l = 200))
  
  summary(croc.ps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
