AkritasNNE <-
function (times, status, covariate, predict.time, cov.value, lambda = NULL, span = NULL, window = "symmetric") {	
	n <- length(covariate)
	
	bad <- is.na(times) | is.na(status) | is.na(covariate)
	
	times <- times[!bad]
	status <- status[!bad]
	covariate <- covariate[!bad]
	
	if (sum(bad) > 0) 
		cat(paste("\n", sum(bad), "records with missing values dropped. \n"))
   
	ord <- order(times)
	times <- times[ord]
	status <- status[ord]
	covariate <- covariate[ord]
	if (is.null(lambda) & is.null(span)) {
		cat("method = NNE requires either lambda or span! \n")
		stop(0)
	}
	covariate.unique <- cov.value
	S.t.covariate <- rep(0, length(covariate.unique))
	t.evaluate <- unique(times[status == 1])
	t.evaluate <- t.evaluate[order(t.evaluate)]
	t.evaluate <- t.evaluate[t.evaluate <= predict.time]
	for (i in 1:length(covariate.unique)) {
		wt <- NNE.weights(covariate, covariate.unique[i], span = span, lambda = lambda, window = window)$w		
		s0 <- 1
		for (k in 1:length(t.evaluate)) {
			n <- sum(wt*(times >= t.evaluate[k]))
			d <- sum(wt*((times == t.evaluate[k])*(status == 1)))
			if (!is.na(n) & n > 0)
				s0 <- s0 * (1 - d/n)
		}
		S.t.covariate[i] <- s0
	}
	res <- list(S = S.t.covariate, lambda = lambda, span = span)
	res
}
