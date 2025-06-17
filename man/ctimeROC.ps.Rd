\name{ctimeROC.ps}
\alias{ctimeROC.ps}
\title{
Penalised-based estimation of the cumulative/dynamic time-dependent ROC curve conditional on covariates.
}
\description{
This function implements the penalised-based approach to the estimation of the conditional cumulative/dynamic time-dependent ROC curve proposed by Rodriguez-Alvarez and Inacio (2025).
}
\usage{
ctimeROC.ps(formula.hazard, formula.biomarker, data, cutoffs = NULL, 
	predict.time, cov.values, n.timepoints= NULL, select = FALSE, 
	fitted.models = NULL)
}
\arguments{
  \item{formula.hazard}{Formula for the model of the conditional hazard function. It must be in the form (with the left-hand side being the important part): Surv(t, delta) ~ t + y + s(x1) + x2, where t is the observed/survival time, delta is the censoring indicator, y is the biomarker, and x1 and x2 are the covariates. Nonlinear effects can be incorporated for all variables in the model. For more details, see the \code{mgcv} package, on which this function relies.}
  \item{formula.biomarker}{Named list containing the formulas for the location-scale regression model of the biomarker. The model for the mean function (including the response/biomarker) should be specified under the element mean, and the model for the variance function under var. In the latter case, a right-hand side formula should be provided. See the example below.}
  \item{data}{Data frame containing all needed variables (covariates, biomarker, survival time and censoring indicator)}
  \item{cutoffs}{Vector of cutoffs at which the conditional time-dependent True Positive and the False Positive Rates will be computed. If NULL, the cutoffs are the unique values of the marker in the sample.}
  \item{predict.time}{Time-point at which the conditional time-dependent True Positive and the False Positive Rates will be computed. In the current implementation, only one value is allowed.}
  \item{cov.values}{Data frame with the covariate values at which the conditional time-dependent True Positive and the False Positive Rates will be computed}
  \item{n.timepoints}{For the data augmentation strategy in the piecewise exponential approach, the number of equally spaced breakpoints. If NULL (default), breakpoints are set at the uncensored times. If NULL (default), the break points are set at the uncensored times.}
  \item{select}{Logical. If TRUE, the double penalty approach for model building/selection described in Marra and Wood (2011) is implemented. See \code{mgcv} package for details.}
  \item{fitted.models}{Named list containing the fitted models for the conditional hazard function (hazard) and the biomarker (biomarker). If provided, it can save time by preventing the models from being refitted. See the example below.}
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
\item{fitted.models}{Named list with the fitted models for the conditional hazard function (hazard) and the biomarker (biomarker).}
}
\references{
G. Marra and S.N. Wood (2011). Practical variable selection for generalized additive models. \emph{Computational Statistics and Data Analysis}, \bold{55}, 2372-2387.

M.X. Rodriguez-Alvarez and V. Inacio (2025). Penalised spline estimation of covariate-specific time-dependent ROC curves. ArXiv preprint (https://arxiv.org/abs/2506.13604)
}
\author{
Maria Xose Rodriguez-Alvarez and Vanda Inacio
}
%\note{
%}

\seealso{
\code{\link{ctimeROC.sp}}, \code{\link{ctimeROC.np}}
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
	
	# Penalised-based approach
	\donttest{
	croc.ps <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ s(y) + s(x) + s(t)", 
		formula.biomarker = list(mean = y ~ s(x), var = ~ s(x)), data = data,
		cutoffs = seq(-3, 5, l = 200), 
		predict.time = 1,
		cov.values = data.frame(x = 0))
	
	croc.ps$AUC
	
	# Plot the ROC curve
	plot(croc.ps$FPR, croc.ps$TPR, xlab = "FPR", 
		ylab = "TPR", type = "l", lwd = 2)
	abline(0, 1, col = "grey")
	}
	\dontshow{
	croc.ps <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ y + x + t", 
		formula.biomarker = list(mean = y ~ x, var = ~ x), data = data, 
		predict.time = 1, 
		cov.values = data.frame(x = 0), 
		cutoffs = seq(-3, 5, l = 200))
	
	croc.ps$AUC
	
	# Plot the ROC curve
	plot(croc.ps$FPR, croc.ps$TPR, xlab = "FPR", 
		ylab = "TPR", type = "l", lwd = 2)
	abline(0,1,col = "grey")}

	# Using the previous results, we now obtain the ROC curve
	# for a different time
	croc.ps.t2 <- ctimeROC.ps(formula.hazard = "Surv(t, delta) ~ s(y) + s(x) + s(t)", 
		formula.biomarker = list(mean = y ~ s(x), var = ~ s(x)), data = data,
		cutoffs = seq(-3, 5, l = 200), 
		predict.time = 2,
		cov.values = data.frame(x = 0),
		fitted.models = croc.ps$fitted.models)

	croc.ps.t2$AUC
	
	# Plot the ROC curve
	plot(croc.ps.t2$FPR, croc.ps.t2$TPR, xlab = "FPR", 
		ylab = "TPR", type = "l", lwd = 2)
	abline(0,1,col = "grey")


}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
