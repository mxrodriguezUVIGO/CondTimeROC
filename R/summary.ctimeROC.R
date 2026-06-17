summary.ctimeROC <-
function(object, ...) {
	res <- list()
	res$call <- object$call
	method <- switch(class(object)[1], "ctimeROC.np" = "Covariate-specific time-dependent ROC curve - Nonparametric", 
								  "ctimeROC.sp" = "Covariate-specific time-dependent ROC curve - Semiparametric", 
								  "ctimeROC.ps" = "Covariate-specific time-dependent ROC curve - Penalised-spline")

	res$method <- method

	if(class(object)[1] == "ctimeROC.np") {
		res$np.results <- list()

		# Estimator
		tmp <- ifelse(object$smooth == "nu", "Unsmoothed", "Smoothed")
		res$np.results$estimator <- paste0("\nKernel Estimator: ", tmp)

		# Bandwidth selector
		if(object$bwsel == "ALbw") {
			tmp <- "ALbw (Altman and Leger, 1995)"
		} else if (object$bwsel == "dpik") {
			tmp <- "dpik (Sheater and Jones, 1991)"
		} else if (object$bwsel == "npbw") {
			tmp <- "npbw (Li et al., 2013)"
		} else if (object$bwsel == "NNE") {
			tmp <- "Nearest neighbor"
		}
		res$np.results$bandwidth.sel <- paste0("\nBandwidth selector: ", tmp)

		# Inverse probability of censoring weights
		if(object$ipcw == "Beran") {
			tmp <- "Beran's estimator (1981)"
		} else if (object$ipcw == "Akritas") {
			tmp <- "Akritas' estimator (1994)"
		}
		res$np.results$ipwc <- paste0("\nInverse probability of censoring weights: ", tmp)

		# Kernel funcion
		if(object$kernel == "gaussian") {
			tmp <- "Gaussian"
		} else if (object$kernel == "epanechnikov") {
			tmp <- "Epanechnikov"
		}
		res$np.results$kernel <- paste0("\nKernel function: ", tmp)

		# Bandwidths
		m <- matrix(ncol = 2, nrow = 1, dimnames = list(c("Bandwidth:"), c("Covariate", "Biomarker")))
		m[1,] <- c(sprintf("%.6f", object$lbd[1]), sprintf("%.6f", object$lbd[2]))

		res$np.results$lbd <- m
	}

	if(class(object)[1] == "ctimeROC.sp") {
		res$hazard 	  <- summary(object$fitted.models$hazard)
		res$biomarker <- summary(object$fitted.models$biomarker)
	}
	if(class(object)[1] == "ctimeROC.ps") {
		res$hazard 	  <- summary(object$fitted.models$hazard)
		res$biomarker_mean <- summary(object$fitted.models$biomarker$fit_mean)
		res$biomarker_variance <- summary(object$fitted.models$biomarker$fit_var)
	}

	# Summary observations
	m <- matrix(ncol = 1, nrow = 2, dimnames = list(c("Number of observations", "Number of missing data"), c("")))
	m[1,] <- c(sprintf("%.0f", nrow(object$data)))
	m[2,] <- c(sprintf("%.0f", sum(object$missing.ind)))
	res$sz <- m

	print.summary.ctimeROC(res, ...)
	invisible(res)   	   	
}
