print.ctimeROC <-
function(x, ...) {
	method <- switch(class(x)[1], "ctimeROC.np" = "Covariate-specific time-dependent ROC curve - Nonparametric", 
								  "ctimeROC.sp" = "Covariate-specific time-dependent ROC curve - Semiparametric", 
								  "ctimeROC.ps" = "Covariate-specific time-dependent ROC curve - Penalised-spline")

	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", method))
	cat("\n")

	invisible(x)   
}
