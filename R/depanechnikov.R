depanechnikov <-
function(x) {	
	res <- (3/4)*(1-x^2)*as.numeric(abs(x)<=1)
	res
}
