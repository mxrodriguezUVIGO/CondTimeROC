simpson <- 
function(fun.int, set.int) {
	l.set.int <- length(set.int)
	integral <- (set.int[l.set.int] - set.int[1])/(l.set.int - 1)/3*(fun.int[1] + fun.int[l.set.int] + 4*sum(fun.int[seq(2,l.set.int - 1, by = 2)]) + 2*sum(fun.int[seq(3, l.set.int - 2, by = 2)]))
}
obtain.ROCCurve <-
function(set.p, FPR, TPR) {	
	ROC <- rep(NA, length(set.p))
	for(i in 1:length(set.p)) {
		n <- length(FPR)
		LowIndex <- max(c(1:n)[(FPR <= set.p[i])])
		HighIndex <- min(c(1:n)[(FPR >= set.p[i])])
		delta.x <- FPR[HighIndex] - FPR[LowIndex]
		delta.y <- TPR[HighIndex] - TPR[LowIndex]
		delta.2.eval <- set.p[i] - FPR[LowIndex]
		target.y <- TPR[LowIndex]
		if(delta.x > 0) target.y <- target.y + (delta.2.eval/delta.x)*delta.y
		ROC[i] <- target.y
	}
	AUC <- simpson(ROC, set.p)
	res <- list(ROC = ROC, AUC = AUC)
	res
}
