plot.ctimeROC <-
function(x, ask = TRUE, ...) {
	change.ROC.format <- function(p, ROC) {
		temp <- reshape(ROC, varying = paste("p", round(p, 3), sep = ""), sep = "",
		v.names = "ROC", timevar = "p", times = p, idvar = "comb", direction = "long")
		temp[order(temp$comb),]
	}

	dots <- list(...)	
	cutoffs <- x$cutoffs
	n.cutoffs <- length(cutoffs)
	n.cov <- length(x$AUC)

	# Main
	main.legend <- paste0("Time point: ", x$predict.time)

	if(is.numeric(x$cov.values)) {
		x$cov.values <- data.frame(x$cov.values)
		x$FPR <- matrix(x$FPR, ncol = 1)
		x$TPR <- matrix(x$TPR, ncol = 1)
	}

	# Compute ROC curves
	# FPF
	p <- seq(0, 1, l = 101)
	n.p <- length(p)
	ROC <- matrix(NA, ncol = n.p, nrow = n.cov)
	for(i in 1:n.cov) {
		ROC[i,] <- obtain.ROCCurve(p, rev(x$FPR[,i]), rev(x$TPR[,i]))$ROC
	}

	colnames(ROC) <- paste("p", round(p, 3), sep = "")
	ROC <- cbind(x$cov.values, ROC)	   
	set.accuracy <- "AUC"
	set.accuracy.legend <- "AUC"
	ind.accuracy <- is.element(set.accuracy, set.accuracy[is.element(set.accuracy, names(x))])

	if(any(ind.accuracy)) { 
		accuracy <- set.accuracy[ind.accuracy]
		accuracy.legend <- set.accuracy.legend[ind.accuracy] 		  
	} else {
		accuracy <- NULL
	}   
	if(!is.null(accuracy)) {
		for (i in 1:length(accuracy)){
			aux <- names(ROC)		   
			ROC <- cbind(ROC, x[[accuracy[i]]])
			names(ROC) <- c(aux, set.accuracy.legend[i]) 
		}
	}
				
	names.cov <- names(ROC[, 1:(ncol(ROC) - n.p - sum(ind.accuracy)), drop = FALSE])
	ind.cat <- unlist(lapply(ROC[ ,names.cov, drop = FALSE], is.factor))
	names(ind.cat) <- names.cov	 
	names.cont <- names.cov[!ind.cat]		
	n.cont <- length(names.cont)
	n.cat <- sum(ind.cat)
	names.cat <- if(n.cat > 0) names.cov[ind.cat]
	
	delete.obs <- duplicated(x$cov.values)
	ROC <- ROC[!delete.obs,,drop = FALSE]		 
	if (n.cont > 1) {
		ROC[, names.cont] <- apply(round(ROC[ , names.cont, drop = FALSE], 3), 2, factor)
	}
	if (n.cat > 0) {
		exp.cat <- unique(ROC[, names.cat, drop = FALSE])
		exp.cat.matrix <- as.matrix(exp.cat)
		dim.exp.cat <- nrow(exp.cat)
		levels.cat <- if(n.cat > 0) lapply(ROC[, names.cat, drop = FALSE], levels)			  
		n.levels <- as.numeric(unlist(lapply(levels.cat, length)))		  
		if(n.cont == 0) {
			ROC.long <- change.ROC.format(p, ROC)
			print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cat, collapse = "+"))),
			data = ROC.long,
			ylim = c(0, 1),
			xlim = c(0, 1),
			main = main.legend,
			#ylim = c(-0.1, 1.05),
			xlab = "FPF",
			ylab = "TPF",
			strip = strip.custom(strip.names = TRUE, strip.levels = TRUE, sep = " = ",
			par.strip.text = list(cex = if(!is.null(dots$cex.par.strip.text)) dots$cex.par.strip.text else 0.75)),
			panel = function(x, y, subscripts) {
				panel.abline(0, 1, lty = 2, col = "grey")
				panel.xyplot(x, y, type = "l")
				for (i in 1:length(set.accuracy)) {
					if(ind.accuracy[i]) { 
						acc.val <- round(unique(ROC.long[subscripts, set.accuracy[i]]),2)						
						ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
						labels = paste(accuracy.legend[i],"=",acc.val), adj = c(1,0.5),
						cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
					}
				}
			}))
		} else {				
			cat.cont <- vector("list", dim.exp.cat)					 
			for(i in 1:dim.exp.cat) {
				cat.cont[[i]] <- vector("list", n.cont)
				for(j in 1:n.cont) {
					ind <- t(apply(ROC[ , names.cat, drop = F], 1, function(x) x == exp.cat[i, ]))				  
					if(dim(ind)[1] == 1) ind <- t(ind)			  
					cat.cont[[i]][[j]] <- unique(ROC[apply(ind, 1, all), names.cont[j]])
				}
			}		 
			n.comb <- c(0, as.numeric(cumsum(unlist(lapply(cat.cont, function(x)cumprod(lapply(x, length))[n.cont]))))*n.p)			 
			ROC.long <- change.ROC.format(p, ROC)
			for (i in 1:dim.exp.cat) {
				if(i > 1 && ask) {
					readline("Press return for next page....")				  
				}
				print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cont, collapse = "+"))),
				data = ROC.long,
				ylim = c(0, 1),
				xlim = c(0, 1),
				main = main.legend,
				xlab = "FPF",
				ylab = "TPF",
				#ylim = c(-0.1,1.05),
				subset = (1 + n.comb[i]):n.comb[i + 1],
				strip = strip.custom(style = 3, strip.names = TRUE, strip.levels = TRUE, sep = " = ",
				par.strip.text = list(cex = if(!is.null(dots$par.strip.text)) dots$par.strip.text else 0.75)),
				panel = function(x, y, subscripts) {
					panel.abline(0, 1, lty = 2, col = "grey")
					panel.xyplot(x, y, type = "l")
					for (j in 1:length(set.accuracy)) {
						if(ind.accuracy[j]) { 
							acc.val <- round(unique(ROC.long[subscripts, set.accuracy[j]]),2)
							ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (j-1),
							labels = paste(accuracy.legend[j],"=",acc.val), adj = c(1,0.5),
							cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 1)
						}
					}
				},
				main = paste(names.cat, "=", exp.cat.matrix[i, ])))
			}
		}
	} else {		
		ROC.long <- change.ROC.format(p, ROC)
		print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cont, collapse = "+"))),
		data = ROC.long,
		ylim = c(0, 1),
		xlim = c(0, 1),
		main = main.legend,
		xlab = "FPF",
		ylab = "TPF",
		strip = strip.custom(style = 3, strip.names = TRUE, strip.levels = TRUE, sep = " = ",
		par.strip.text = list(cex = if(!is.null(dots$par.strip.text)) dots$par.strip.text else 0.75)),
		panel = function(x, y, subscripts) {
			panel.abline(0, 1, lty = 2, col = "grey")
			panel.xyplot(x, y, type = "l")
			for (i in 1:length(set.accuracy))
				if(ind.accuracy[i]) { 
					acc.val <- round(unique(ROC.long[subscripts, set.accuracy[i]]),2)						
					ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
					labels = paste(accuracy.legend[i],"=",acc.val), adj = c(1, 0.5),
					cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
				}
		 }
		))
	}
}
