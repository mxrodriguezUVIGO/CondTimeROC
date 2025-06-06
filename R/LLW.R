LLW <-
function(covariate, x, kernel = "gaussian", lbd) {	
	len <- length(covariate)
	listg  <- .C(LLWeightsKernel, as.double(covariate), as.integer(nobs), as.double(x), as.double(lbd), as.character(kernel), w = double(nobs))
	return(listg$weight)
}
