\name{ctimeROC.np}
\alias{ctimeROC.np}
\title{
Nonparametric estimation of the cumulative/dynamic time-dependent ROC curve conditional on covariates.
}
\description{
This function implements the nonparametric approaches to the estimation of the conditional cumulative/dynamic time-dependent ROC curve proposed by Rodriguez-Alvarez et al. (2016).
}
\usage{
ctimeROC.np(marker, covariate, time, status, data, cutoffs = NULL, predict.time, 
		   cov.value = NULL, smooth = c("nu", "ns"), 
		   bwsel = c("ALbw", "dpik", "npbw", "NNE"), ipcw = c("Beran", "Akritas"), 
		   lbd.x = NULL, lbd.y = NULL, span = NULL, lambda = NULL, 
		   kernel = c("gaussian", "epanechnikov"), 
		   AUC.calc = c("numaprox", "cexpr"))
}
\arguments{
  \item{marker}{Character string with the name of the biomarker variable.}
  \item{covariate}{Character string with the name of the covariate variable.}
  \item{time}{Character string with the name of the survival time.}
  \item{status}{Character string with the name of the censoring indicator variable.}
  \item{data}{Data frame containing all needed variables (covariate, biomarker, observed/survival time and censoring indicator).}
  \item{cutoffs}{Vector of cutoffs at which the conditional time-dependent True Positive and the False Positive Rates will be computed. If NULL, the cutoffs are the unique values of the marker in the sample.}
  \item{predict.time}{Time-point at which the conditional time-dependent True Positive and the False Positive Rates will be computed.  In the current implementation, only one value is allowed.}
  \item{cov.value}{Covariate value at which the conditional time-dependent True Positive and the False Positive Rates will be computed. In the current implementation, only one value is allowed.}
  \item{smooth}{Character string with the nonparametric estimator: "nu": unsmoothed estimator - "ns": smoothed estimator. By default "nu".}
  \item{bwsel}{Character string with the bandwidth selector. "ALbw": Altman and Leger (1995). "dpik": Sheater and Jones (1991). "npbw": Li et al. (2013). "NNE": Nearest neighbor.}
  \item{ipcw}{Character string with the approach for computing inverse probability of censoring weights. "Beran": the Beran's estimator (1981). "Akritas": Akritas' estimator (1994).}
  \item{lbd.x}{Numerical value with the bandwidth for the covariate. If NULL it is computed according to the chosen bandwidth selector.}
  \item{lbd.y}{Numerical value with the bandwidth for the marker. If NULL it is computed according to the chosen bandwidth selector.}
  \item{span}{Numerical value with the span for the NNE. Also needed when using the icpw based on Akritas' approach.}
  \item{lambda}{Numerical value with the smoothing parameter for NNE. Also needed when using the ibpw based on Akritas' approach. Need either lambda or span for NNE.}
  \item{kernel}{Character string with the kernel function. Either "gaussian" or "epanechnikov".}
  \item{AUC.calc}{Character string with the approach for the estimation of the time-dependent AUC. Either numerical integration methods ("numaprox") or the closed-form expression presented in the paper by Rodriguez-Alvarez et al. (2016).}
}
%\details{
%}
\value{
Returns a list with the following components:
\item{cutoffs}{Vector with the cutoffs at which the True Positive Rates and the False Positive Rates have been computed}
\item{TPR}{Vector with the estimated True Positive Rates.}
\item{FPR}{Vector with the estimated False Positive Rates.}
\item{AUC}{Estimated time-dependent AUC.}
\item{lbd}{Used bandwidths.}
}
\references{
M. Akritas (1994). Nearest neighbor estimation of a bivariate distribution under random censoring. \emph{Annals of Statistics} \bold{22}, 1299--1327.

N. Altman and C. Leger (1995). Bandwidth selection for kernel distribution function estimation. \emph{Journal of Statistical Planning and Inference} \bold{46}, 195--214.

R. Beran (1982). Nonparametric Regression with Randomly Censored Survival Data. \emph{Technical report: University of California, Berkeley}.

Q. Li and J. Lin and J. S. Racine (2013). Optimal bandwidth selection for nonparametric conditional distribution and quantile functions. \emph{Journal of Business and Economic Statistics} \bold{31}, 57--65.

M.X. Rodriguez-Alvarez, L. Meira-Machado, E. Abu-Assi and S. Raposeiras-Roubin (2016). Nonparametric estimation of time-dependent ROC curves conditional on a continuous covariate. \emph{Statistics in Medicine} \bold{35}, 1090--1102.

S. J. Sheater and M. C. Jones (1991). A reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society, Series B} \bold{53}, 683--690.
}
\author{
Maria Xose Rodriguez-Alvarez and Luis Filipe Meira-Machado
}
%\note{
%}

\seealso{
\code{\link{ctimeROC.ps}}, \code{\link{ctimeROC.sp}}
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
	
	# Non-parametric approach
	\donttest{
	croc.np <- ctimeROC.np(marker = "y", covariate = "x", time = "t", 
		status = "delta", data = data, cutoffs = seq(-3, 5, l = 200), 
		predict.time = 1, cov.value = 0, smooth = "ns", 
		bwsel = "npbw", AUC.calc = "cexpr")
	
	croc.np$AUC
	
	# Plot the ROC curve
	plot(croc.np$FPR, croc.np$TPR, xlab = "FPR", 
		ylab = "TPR", type = "l", lwd = 2)
	abline(0,1,col = "grey")}

	\dontshow{
	croc.np <- ctimeROC.np(marker = "y", covariate = "x", time = "t", 
		status = "delta", data = data, cutoffs = c(-3,5), 
		predict.time = 1, cov.value = 0, 
		bwsel = "dpik", AUC.calc = "cexpr")
	
	croc.np$AUC
	
	# Plot the ROC curve
	plot(croc.np$FPR, croc.np$TPR, xlab = "FPR", 
		ylab = "TPR", type = "l", lwd = 2)
	abline(0,1,col = "grey")}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
