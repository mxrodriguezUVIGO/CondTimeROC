NNE.weights <-
function(x, xx, span = NULL, lambda = NULL, window = "symmetric") {
	n <- length(x)
	if (!is.null(span)) {
		if (window == "symmetric") {
			ddd <- (x - xx)
			ddd <- ddd[order(ddd)]
			index0 <- sum(ddd < 0) + 1
			index1 <- index0 + trunc(n * span/2)
			if (index1 > n) 
				index1 <- n
			lambda <- ddd[index1]
			wt <- as.integer(((x - xx) <= lambda) & ((x - xx) >= 0))
			index0 <- sum(ddd <= 0)
			index2 <- index0 - trunc(n * span/2)
			if (index2 < 1) 
				index2 <- 1
			lambda <- abs(ddd[index2])
			set.index <- ((x - xx) >= -lambda) & ((x - xx) <= 0)
			wt[set.index] <- 1
		}
		if (window == "asymmetric") {
			ddd <- (x - xx)
			ddd <- ddd[order(ddd)]
			index0 <- sum(ddd < 0) + 1
			index <- index0 + trunc(n * span)
			if (index > n) 
				index <- n
			lambda <- ddd[index]
			wt <- as.integer(((x - xx) <= lambda) & ((x - xx) >= 0))
		}
	} else {
		wt <- exp(-(x - xx)^2/lambda^2)
	}
	res <- list(w = wt, span = span, lambda = lambda)
	res
}
