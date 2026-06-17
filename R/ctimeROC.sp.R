ctimeROC.sp <-
function(formula.coxph, formula.lm, data, cutoffs = NULL, predict.time, cov.values) {
	# Obtain the marker, covariates, survival time and status
	if(is.character(formula.coxph)) {
		formula.coxph <- as.formula(formula.coxph)
	}
	if(is.character(formula.lm)) {
		formula.lm <- as.formula(formula.lm)
	}
	tf_lm <- terms.formula(formula.lm)
	if (attr(tf_lm, "response") > 0) {
		marker <- as.character(attr(tf_lm, "variables")[2])
	} else {
		stop("The formula for the model of the biomarker should include the response variable (left hand side)")
	}
	vars_lm <- get_vars_formula(formula.lm)
	vars_coxph <- get_vars_formula(formula.coxph) # May include time
	
	time_var <- all.vars(formula.coxph)[1]
	delta_var <- all.vars(formula.coxph)[2]

	all_vars <- unique(c(marker, time_var, delta_var, vars_lm, vars_coxph))
	all_cov  <- all_vars[!(all_vars %in% c(marker, time_var, delta_var))]

	# Old code
	#marker 	 <- all.vars(formula.lm)[1] #as.character(attr(terms.formula(formula.lm), "variables")[2])
	#terms.lm <- all.vars(formula.lm) 	#attr(terms.formula(formula.lm), "term.labels")
	#terms.coxph <- all.vars(formula.coxph)
	##terms.coxph <- attr(terms.formula(formula.coxph), "term.labels")	
	##terms.surv <- attr(terms.formula(formula.coxph), "variables")
	##cov.names <- unique(c(terms.lm, terms.coxph, as.character(terms.surv[[2]][[2]]), as.character(terms.surv[[2]][[3]])))
	#cov.names <- unique(c(terms.lm, terms.coxph))
	#cov.names <- cov.names[!(cov.names %in% marker)]

	if (inherits(data, what = 'data.frame')) {
        data <- as.data.frame(data)
    } else {
        stop("The object specified in argument 'data' is not a data frame")
    }
	# Check that all variables are in the data
	if(sum(is.na(match(all_vars, names(data))))) {
        stop("Not all needed variables are supplied in data")
    }

    # Check that all covariates are in cov.values
	if(sum(is.na(match(all_cov, names(cov.values))))) {
        stop("Not all needed variables are supplied in cov.values")
    }

    data.new <- data[,all_vars]
    omit <- apply(data.new, 1, anyNA)

	# Delete NAs
    data.new <- data.new[!omit,,drop = FALSE]

	# Delete NAs
	# data <- data[,all_vars]
	# data <- na.omit(data)

	# The cutoffs to compute the ROC curve	
	if(is.null(cutoffs)) {
		cutoffs <- sort(unique(data.new[,marker]))
	}
	
	# Now is a data.frame (new version)
	# cov.values <- matrix(cov.values, ncol = 1, nrow = length(cov.values))
	# colnames(cov.values) <- attr(terms.formula(formula.lm), "term.labels")
	# cov.values <- data.frame(cov.values)

	# TPR and FPR (we assume that only one predict.time is considered)
	TPR <- FPR <- matrix(ncol = nrow(cov.values), nrow = length(cutoffs))
	 
	# Fit the coxph model	
	fit.coxph <- coxph(formula.coxph, data = data.new)
	
	# Fit linear regression model
	fit.lm <- lm(formula.lm, data = data.new)
	residuals <- fit.lm$model[,marker] - fit.lm$fitted
	
	# Obtain the log of the cumulative baseline hazard function
	#sfit <- survfit(fit.coxph, newdata = data.new)
	sfit <- survfit(fit.coxph)
	h <- log(-log(sfit$surv))
	h_time <- approxfun(sfit$time, h, yleft = min(h), yright = max(h))(predict.time) # log(Delta_0(predict.time))
	
	# Predict on the new dataframe
	p.lm <- predict(fit.lm, newdata = cov.values)
	
	# Obtain residuals + \alpha^T X
	# as many columns as covariate values
	# as many rows as sample size
	p.roc <- outer(residuals, p.lm, '+')
	
	# Conditional survival function
	# as many columns as covariate values	
	# as many rows as sample size
	p.coxph <- do.call("cbind", sapply(1:ncol(p.roc), function(i, fit, data, data.roc) {
		df <- cov.values[i,,drop = FALSE]
		df <- do.call("rbind", replicate(length(p.roc[,i]), df, simplify = FALSE))
		df[,marker] <- p.roc[,i]
		p.aux <- predict(fit, newdata = df)
		surv.aux <- list(exp(-exp(h_time + p.aux)))
	}, fit = fit.coxph, data = cov.values, data.roc = p.roc))
	
	den.fpr <- apply(p.coxph, 2, sum)
	den.tpr <- length(residuals) - den.fpr
	
	num.fpr <- t(sapply(1:length(cutoffs), function(i){
		apply((p.roc >= cutoffs[i])*p.coxph, 2, sum)	
	}))
	num.tpr <- t(sapply(1:length(cutoffs), function(i){
		apply((p.roc >= cutoffs[i])*(1-p.coxph), 2, sum)	
	}))
	
	if(ncol(num.fpr) != nrow(cov.values)) num.fpr = t(num.fpr)
	if(ncol(num.tpr) != nrow(cov.values)) num.tpr = t(num.tpr)
	
	FPR <- rbind(rep(1,nrow(cov.values)),t(t(num.fpr)/den.fpr),rep(0,nrow(cov.values)))
	TPR <- rbind(rep(1,nrow(cov.values)),t(t(num.tpr)/den.tpr),rep(0,nrow(cov.values)))

	# AUC
	area <- vector(length = nrow(cov.values))
	names(area) <- rownames(cov.values)
	for(i in 1:length(area)) { 
		area[i] <- obtain.ROCCurve(seq(0,1, l = 101), rev(FPR[,i]), rev(TPR[,i]))$AUC
	}
	
	res <- list()
	res$call 			<- match.call()
	res$data 			<- data
	res$missing.ind 	<- omit
    res$marker 			<- marker
    res$covariates		<- all_cov
    res$time			<- time_var
    res$status			<- delta_var
	res$cutoffs 		<- c(-Inf, cutoffs, Inf)
	res$predict.time 	<- predict.time
	res$cov.values		<- cov.values
	res$TPR 			<- TPR
	res$FPR 			<- FPR
	res$AUC 			<- area
	res$fitted.models 	<- list(hazard = fit.coxph, biomarker = fit.lm)
	class(res) <- c("ctimeROC.sp", "ctimeROC")
	res
}
