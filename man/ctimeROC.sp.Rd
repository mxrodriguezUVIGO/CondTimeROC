\name{ctimeROC.sp}
\alias{ctimeROC.sp}
\title{
Semiparametric estimation of the cumulative/dynamic time-dependent ROC curve conditional on covariates.
}
\description{
This function implements the semiparametric approach to the estimation of the conditional cumulative/dynamic time-dependent ROC curve proposed by Song and Zhou (2008).
}
\usage{
ctimeROC.sp(formula.coxph, formula.lm, data, cutoffs = NULL, predict.time, cov.values)
}
\arguments{
  \item{formula.coxph}{Formula for the Cox model associated with the proposal by Song and Zhou (2008). It has to be of the form: "Surv(t, delta) ~ y + x", with x being the covariate and y the biomarker.}
  \item{formula.lm}{Formula for the linear model associated with the proposal by Song and Zhou (2008). It has to be of the form: "y ~ x", with x being the covariate and y the biomarker.}
  \item{data}{Data frame containing all needed variables (covariate, biomarker, survival time and censoring indicator)}
  \item{cutoffs}{Vector of cutoffs at which the conditional time-dependent True Positive and the False Positive Rates will be computed. If NULL, the cutoffs are the unique values of the marker in the sample.}
  \item{predict.time}{Time-point at which the conditional time-dependent True Positive and the False Positive Rates will be computed. In the current implementation, only one value is allowed.}
  \item{cov.values}{Data.frame with the covariates values at which the conditional time-dependent True Positive and the False Positive Rates will be computed.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a list with the following components:
\item{cutoffs}{Vector with the cutoffs at which the True Positive Rates and the False Positive Rates have been computed}
\item{TPR}{Matrix with as many columns as covariate values and as many rows as cutoffs with the estimated True Positive Rates.}
\item{FPR}{Matrix with as many columns as covariate values and as many rows as cutoffs with the estimated False Positive Rates.}
\item{AUC}{Vector with the estimated time-dependent AUC at each covariate value.}
}
\references{
X. Song and X-H. Zhou (2008). A semiparametric approach for the covariate specific ROC curve with survival outcome. \emph{Statistica Sinica} \bold{18}, 947--965.
}
\author{
Maria Xose Rodriguez-Alvarez and Luis Filipe Meira-Machado
}
%\note{
%}

\seealso{
\code{\link{ctimeROC.ps}}, \code{\link{ctimeROC.np}}
}
\examples{
	set.seed(123)
	# Simulate the data
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
	
	# Song and Zhou approach
	croc.sz <- ctimeROC.sp(formula.coxph = "Surv(t, delta) ~ y + x", 
		formula.lm = "y ~ x", data = data, predict.time = 1, 
		cov.values = data.frame(x = 0), cutoffs = seq(-3,5, l = 200))
	
	croc.sz$AUC
	
	# Plot the ROC curve
	plot(croc.sz$FPR, croc.sz$TPR, xlab = "FPR", 
		ylab = "TPR", type = "l", lwd = 2)
	abline(0, 1, col = "grey")
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
