BeranKernel <-
function(time, status, covariate, predict.time, cov.value, kernel, lbd) {
	len <- length(time)
	delta <- rep(1, len)
	return(.C(SurvBeranKernel, as.double(time), as.integer(status), as.double(covariate), as.integer(delta), as.integer(len), as.double(predict.time), as.double(cov.value), as.double(lbd), as.character(kernel), p = as.double(1))$p)
}
