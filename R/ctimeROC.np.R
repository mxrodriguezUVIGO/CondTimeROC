ctimeROC.np <-
function(marker, covariate, time, status, data, cutoffs = NULL, predict.time, cov.value = NULL, smooth = c("nu","ns"), bwsel = c("ALbw","dpik","npbw", "NNE"), ipcw = c("Beran", "Akritas"), lbd.x = NULL, lbd.y = NULL, span = NULL, lambda = NULL, kernel = c("gaussian", "epanechnikov"), AUC.calc = c("numaprox","cexpr")) {
	
	smooth = match.arg(smooth)
	bwsel = match.arg(bwsel)
	kernel = match.arg(kernel)
	ipcw = match.arg(ipcw)
	AUC.calc = match.arg(AUC.calc)
	
	if ((bwsel == "NNE" || ipcw == "Akritas") & is.null(lambda) & is.null(span)) {
		cat("NNE and/or Akritas requires either lambda or span! \n")
		stop(0)
	}
	# Delete NAs
	data <- data[ ,c(marker, covariate, time, status)]
	data <- na.omit(data)
	
	x <- data[,covariate]
	y <- data[,marker]
	ttilde <- data[,time]
	status <- data[,status]	
	
	n <- length(y)
	
	xx <- cov.value
	pvt <- predict.time
	
	tc <- if(is.null(cutoffs)) {
		unique(sort(y))
	} else {
		cutoffs
	}
	nc <- length(tc)
	FPR <- TPR <- rep(0, nc)	
	
	if(is.null(lbd.x)) {
		lbd.x <- ifelse(bwsel == "ALbw", ALbw(type_kernel = switch(kernel, "gaussian"="n", "epanechnikov"="e"), vec_data = x[]), 
			ifelse(bwsel == "dpik", dpik(x[], scalest = "iqr", level = 1L, kernel = switch(kernel, "gaussian"="normal", "epanechnikov"="epanech")), -1))		
	}
	if(is.null(lbd.y) & smooth  == "ns") {
		lbd.y <- ifelse(bwsel == "ALbw", ALbw(type_kernel = switch(kernel, "gaussian"="n", "epanechnikov"="e"), vec_data = y[]), 
			ifelse(bwsel == "dpik", dpik(y[], scalest = "iqr", level = 1L, kernel = switch(kernel, "gaussian"="normal", "epanechnikov" = "epanech")), -1))
	}
	
	if(bwsel == "npbw" && smooth == "nu" && lbd.x == -1) {
		bw <- npcdistbw(y~x, data.frame(x=x, y=y), cxkertype = kernel)
		lbd.x <- bw$xbw
	} else if(bwsel == "npbw" && smooth == "ns" && lbd.x == -1 && lbd.y == -1) {				
		bw <- npcdistbw(y~x, data.frame(x=x, y=y), cxkertype = kernel)
		lbd.x <- bw$xbw
		lbd.y <- bw$ybw				
	}
	
	if((!identical(bwsel, "NNE") & smooth == "nu" & lbd.x == -1) || (!identical(bwsel, "NNE") & smooth == "ns" && (lbd.x == -1 || lbd.y == -1))) {
		stop("Incorrect bandwidth selector")
	}
	
	ipcw.FPR <- ipcw.TPR <- vector(length = n)
	
	if(!is.null(xx)) {
		if(identical(bwsel, "NNE")) {  # NNE weights
			wx <- NNE.weights(x, xx, span = span, lambda = lambda, window = "symmetric")$w
		} else { # Kernel weights
			wx <- lbd.x^{-1}*eval(parse(text = switch(kernel, "gaussian"="dnorm", "epanechnikov"="depanechnikov")))((xx-x)/lbd.x)
		}		
		wx <- wx/sum(wx)
	} else {
		wx <- rep(1, n)
	}			
	wt.fpr <- wt.tpr <- ttilde <= predict.time
	
	# IPCW
	for (k in 1:n){
		if(ipcw == "Beran") {
			if(identical(bwsel, "NNE")) {
				ipcw.TPR[k] <- BeranNNE(ttilde[], 1-status[], x[], ttilde[k], x[k], span = span, lambda = lambda)
				ipcw.FPR[k] <- BeranNNE(ttilde[], 1-status[], x[], predict.time, x[k], span = span, lambda = lambda)
			} else { 
				ipcw.TPR[k] <- BeranKernel(ttilde[], 1-status[], x[], ttilde[k], x[k], kernel = kernel, lbd = lbd.x)
				ipcw.FPR[k] <- BeranKernel(ttilde[], 1-status[], x[], predict.time, x[k], kernel = kernel, lbd = lbd.x)
			}
		} else {
			ipcw.TPR[k] <- AkritasNNE(ttilde[], 1-status[], x[], ttilde[k], x[k], span = span, lambda = lambda)$S
			ipcw.FPR[k] <- AkritasNNE(ttilde[], 1-status[], x[], predict.time, x[k], span = span, lambda = lambda)$S
		}
	}
	for (i in 1:nc){		
		if(smooth == "ns") {  	
			wy <- eval(parse(text = switch(kernel, "gaussian"="pnorm", "epanechnikov"="pepanechnikov")))((tc[i] - y[])/lbd.y)
		} else {
			wy <- y <= tc[i]
		}
		num.tpr <- wt.tpr*(1-wy)*wx*status*(1/ipcw.TPR)
		num.tpr <- replace(num.tpr, !is.finite(num.tpr), 0)
	
		den.tpr <- wt.tpr*wx*status*(1/ipcw.TPR)
		den.tpr <- replace(den.tpr, !is.finite(den.tpr), 0)
			
		TPR[i] <- sum(num.tpr, na.rm = TRUE)/sum(den.tpr, na.rm = TRUE)
				
		num.fpr <- (1-wt.fpr)*(1-wy)*wx*(1/ipcw.FPR)
		num.fpr <- replace(num.fpr, !is.finite(num.fpr), 0)
	
		den.fpr <- (1-wt.fpr)*wx*(1/ipcw.FPR)
		den.fpr <- replace(den.fpr, !is.finite(den.fpr), 0)
	
		FPR[i] <- sum(num.fpr, na.rm = TRUE)/sum(den.fpr, na.rm = TRUE)				
	}
	TPR <- replace(TPR, !is.finite(TPR), 0)
	FPR <- replace(FPR, !is.finite(FPR), 0)
	
	TPR <- c(1, TPR, 0)
	FPR <- c(1, FPR, 0)
	
	# AUC
	if(AUC.calc == "cexpr") {
		if(smooth == "nu" & is.null(cutoffs)) {
			n = length(FPR)
			dx <- FPR[-n] - FPR[-1]
			mid.y <- (TPR[-n] + TPR[-1])/2
			area <- sum(dx * mid.y)
		} else {
			area <- AUC.calculation(smooth = smooth, lbd.y = lbd.y, wt = wt.tpr, wx = wx, ipcw.FPR = ipcw.FPR, ipcw.TPR = ipcw.TPR, status = status, marker = y, kernel = kernel)
		}
	} else {
		area <- obtain.ROCCurve(seq(0,1, l = 101), rev(FPR), rev(TPR))$AUC
	}
	
	res <- list(cutoffs = c(-Inf, tc, Inf), TPR = TPR, FPR = FPR, AUC = area, lbd = c(ifelse(is.null(lbd.x), -1, lbd.x), ifelse(is.null(lbd.y), -1, lbd.y)))
	res
}
