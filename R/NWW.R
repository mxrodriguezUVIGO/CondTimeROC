NWW <-
function(covariate, x, kernel = "gaussian", lbd) {
	len <- length(covariate)	
	listg <- .C(NWWeightsKernel, as.double(covariate), as.integer(len), as.double(x), as.double(lbd), as.character(kernel), weight = double(len))
	return(listg$weight)
}
