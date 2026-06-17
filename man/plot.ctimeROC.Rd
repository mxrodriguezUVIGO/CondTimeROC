\name{plot.ctimeROC}
\alias{plot.ctimeROC}
\title{
Default ctimeROC plotting
}
\description{
Takes a fitted \code{ctimeROC} object produced by \code{ctimeROC.np()}, \code{ctimeROC.sp()}, or \code{ctimeROC.ps()} and plots the covariate-specific time-dependent ROC curve (ctimeROC) and associated area under the ctimeROC curve. The suitable type of graphic is chosen according to the number and nature of the covariates.
}
\usage{
\method{plot}{ctimeROC}(x, ask = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{ctimeROC} as produced by \code{ctimeROC.np()}, \code{ctimeROC.sp()}, or \code{ctimeROC.ps()}.}
  \item{ask}{A logical value. If TRUE, the default, the user is asked for confirmation, before a new figure is drawn.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
%\details{
%TODO
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ctimeROC.np}}, \code{\link{ctimeROC.sp}} or \code{\link{ctimeROC.ps}}
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
  
  croc.ps$AUC
  
  # Plot the ROC curve manually
  plot(croc.ps$FPR, croc.ps$TPR, xlab = "FPR", 
    ylab = "TPR", type = "l", lwd = 2)
  abline(0, 1, col = "grey")
  
  # Plot the ROC curve using the provided plot method
  plot(croc.ps)
  }
  \dontshow{
    croc.ps <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ y + x + t", 
      formula.biomarker = list(mean = y ~ x, var = ~ x), data = data, 
      predict.time = 1, 
      cov.values = data.frame(x = 0), 
      cutoffs = seq(-3, 5, l = 200))
    
    croc.ps$AUC
    
    # Plot the ROC curve manually
    plot(croc.ps$FPR, croc.ps$TPR, xlab = "FPR", 
      ylab = "TPR", type = "l", lwd = 2)
    abline(0,1,col = "grey")

    # Plot the ROC curve using the provided plot method
    plot(croc.ps)
  }

  # Using the previous results, we now obtain the ROC curve
  # for a different time
  croc.ps.t2 <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ s(y) + s(x) + s(t)", 
    formula.biomarker = list(mean = y ~ s(x), var = ~ s(x)), data = data,
    cutoffs = seq(-3, 5, l = 200), 
    predict.time = 2,
    cov.values = data.frame(x = 0),
    fitted.models = croc.ps$fitted.models)

  croc.ps.t2$AUC
  
  # Plot the ROC curve manually
  plot(croc.ps.t2$FPR, croc.ps.t2$TPR, xlab = "FPR", 
    ylab = "TPR", type = "l", lwd = 2)
  abline(0,1,col = "grey")

  # Plot the ROC curve using the provided plot method
  plot(croc.ps.t2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
