print.summary.ctimeROC <-
function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat(paste0("\nApproach: ", x$method))

	if(!is.null(x$hazard)) {
		cat("\n----------------------------------------------------------\n")
		cat("Model for the conditional hazard function")
		cat("\n----------------------------------------------------------\n")
		print(x$hazard)
	}
	if(!is.null(x$biomarker)) {
		cat("\n----------------------------------------------------------\n")
		cat("Model for the biomarker")
		cat("\n----------------------------------------------------------\n")
		print(x$biomarker)
	}
	if(!is.null(x$biomarker_mean)) {
		cat("\n----------------------------------------------------------\n")
		cat("Model for the biomarker (mean)")
		cat("\n----------------------------------------------------------\n")
		print(x$biomarker_mean)
	}

	if(!is.null(x$biomarker_var)) {
		cat("\n----------------------------------------------------------\n")
		cat("Model for the biomarker (variance)")
		cat("\n----------------------------------------------------------\n")
		print(x$biomarker_var)
	}

	if(!is.null(x$np.results)) {
		cat("\n")
		cat(x$np.results$estimator)
		cat(x$np.results$bandwidth.sel)
		cat(x$np.results$ipwc)
		cat(x$np.results$kernel)
		cat("\n\n")
		print(x$np.results$lbd, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	if(!is.null(x$sz)) {
		cat("\n----------------------------------------------------------\n")
		cat("Sample sizes")
		cat("\n----------------------------------------------------------\n")
		print(x$sz, quote = FALSE, right = TRUE, na.print = "", print.gap = 5)
	}
	invisible(x)
}
