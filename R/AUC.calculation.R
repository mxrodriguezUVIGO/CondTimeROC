AUC.calculation <-
function(smooth = c("nu", "ns"), lbd.y, wt, wx, ipcw.FPR, ipcw.TPR, status, marker, kernel = c("gaussian", "epanechnikov")) {
	smooth = match.arg(smooth)
	kernel = match.arg(kernel)
	n <- length(marker)
	FPR <- vector(length = n)
	for(i in 1:n) {
		if(smooth == "ns") {  	
			wy <- 1 - eval(parse(text = switch(kernel, "gaussian"="pnorm", "epanechnikov"="pepanechnikov")))((marker[i] - marker[])/lbd.y)
		} else {
			wy <- as.numeric((marker > marker[i])) + 0.5*as.numeric((marker == marker[i]))
		}
		num.fpr <- (1-wt)*(1-wy)*wx*(1/ipcw.FPR)
		num.fpr <- replace(num.fpr, !is.finite(num.fpr), 0)
	
		den.fpr <- (1-wt)*wx*(1/ipcw.FPR)
		den.fpr <- replace(den.fpr, !is.finite(den.fpr), 0)
	
		FPR[i] <- 1 - sum(num.fpr, na.rm = TRUE)/sum(den.fpr, na.rm = TRUE)
	}
	num.auc <- status*wt*wx*(1-FPR)/ipcw.TPR
	num.auc <- replace(num.auc, !is.finite(num.auc), 0)
	
	den.auc <- status*wt*wx/ipcw.TPR
	den.auc <- replace(den.auc, !is.finite(den.auc), 0)
	
	AUC <- sum(num.auc)/sum(den.auc)
	AUC	
}
