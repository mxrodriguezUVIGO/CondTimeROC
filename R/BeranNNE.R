BeranNNE <-
function(time, status, covariate, predict.time, cov.value, lambda = NULL, span = NULL, window = "symmetric") {   #P(T>t|X=x)
	len <- length(time)
	delta <- rep(1, len)
	if (is.null(lambda) & is.null(span)) {
		cat("method = NNE requires either lambda or span! /n")
		stop(0)
	}
	if (!is.null(lambda) & !is.null(span)) {
		cat("method = NNE requires either lambda or span! /n")
		stop(0)
	}
	return(.C(SurvBeranNNE, as.double(time), as.integer(status), as.double(covariate), as.integer(delta), as.integer(len), as.double(predict.time), as.double(cov.value), as.double(span), as.character("symmetric"), p = as.double(1))$p)
}
