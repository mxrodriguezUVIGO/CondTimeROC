#####################################################################################################
# Obtain the variables in a (possibly mgcv) formula
#####################################################################################################
get_vars_formula <-
function(formula) {
	gp <- mgcv::interpret.gam(formula)
	vars.formula <- gp$pred.names
	vars.formula
}
#####################################################################################################
# Location and scale model: mgcv
#####################################################################################################
loc_scale_model_mgcv <- function(formulas, data, select = TRUE) {
	fit <- gam(as.formula(formulas$mean), data = data, method = 'REML', select = select)
	reg_fun_fitted <- fit$fitted.values
	
	# Extract response variable from formulas$mean
	response <- as.character(attr(terms.formula(as.formula(formulas$mean)), "variables")[2])
	if(!is.null(formulas$var)) {
		data$working_responses <- log((data[,response] - reg_fun_fitted)^2)
		#formulas$var <- paste0("working_responses", as.character(formulas$var))
		formulas$var <- update.formula(formulas$var, working_responses ~ .)
		fit_wr <- gam(as.formula(formulas$var), data = data, method = 'REML', select = select)
		theta_hat <- sum(((data[,response] - reg_fun_fitted)^2)*exp(fit_wr$fitted.values))/sum((exp(fit_wr$fitted.values))^2)
		var_fun_fitted <- theta_hat*exp(fit_wr$fitted.values)

		std_residuals <- (data[,response] - reg_fun_fitted)/sqrt(var_fun_fitted)
	} else {
		fit_wr <- NULL
		theta_hat <- NULL
		std_residuals <- data[,response] - reg_fun_fitted
	}

	res <- list(fit_mean = fit, fit_var = fit_wr, theta = theta_hat, std_residuals = std_residuals)
	res
}
predict_loc_scale_model_mgcv <- function(object, newdata) {
	reg_fun_pred <- predict(object$fit_mean, newdata)
	if(!is.null(object$fit_var)) {
		sd_fun_pred <- sqrt(object$theta*exp(predict(object$fit_var, newdata)))
	} else {
		sd_fun_pred <- rep(1, nrow(newdata))
	}
	res <- list()
	res$p_mean <- reg_fun_pred
	res$p_sd <- sd_fun_pred
	res
}
#####################################################################################################
# Location and scale model: mgcv
#####################################################################################################
ctimeROC.ps <- function(formula.hazard, formula.biomarker, data,
	cutoffs = NULL, 
	predict.time, 
	cov.values,
	n.timepoints = NULL,
	select = FALSE,	
	fitted.models = NULL) {
	
	if(!is.list(formula.biomarker)) {
		aux <- list()
		aux$mean <- formula.biomarker
		aux$var <- NULL
		formula.biomarker <- aux
	}
	if(is.character(formula.hazard)) {
		formula.hazard <- as.formula(formula.hazard)
	}
	if(is.character(formula.biomarker$mean)) {
		formula.biomarker$mean <- as.formula(formula.biomarker$mean)
	}
	if(is.character(formula.biomarker$var)) {
		formula.biomarker$var <- as.formula(formula.biomarker$var)
	}

	# Obtain marker, time, delta and covariates for processing
		tf_mean <- terms.formula(formula.biomarker$mean)
		if (attr(tf_mean, "response") > 0) {
			marker <- as.character(attr(tf_mean, "variables")[2])
		} else {
			stop("The formula for the mean function of the biomarker should include the response variable (left hand side)")
		}
		vars_mean <- get_vars_formula(formula.biomarker$mean)
		vars_variance <- get_vars_formula(formula.biomarker$var)
		vars_hazard <- get_vars_formula(formula.hazard) # May include time
		
		time_var <- all.vars(formula.hazard)[1]
		delta_var <- all.vars(formula.hazard)[2]

		all_vars <- unique(c(marker, time_var, delta_var, vars_mean, vars_variance, vars_hazard))
		all_cov  <- all_vars[!(all_vars %in% c(marker, time_var, delta_var))]

	# Check that all variables are in the data
	if(sum(is.na(match(all_vars, names(data))))) {
        stop("Not all needed variables are supplied in data")
    }

    # Check that all covariates are in cov.values
	if(sum(is.na(match(all_cov, names(cov.values))))) {
        stop("Not all needed variables are supplied in cov.values")
    }

	# Delete NAs
	data <- data[,all_vars]
	data <- na.omit(data)	
	
	# Surv for pamm
	vars_hazard_wo_time <- vars_hazard[!(vars_hazard %in% time_var)]  	
	tf <- terms.formula(formula.hazard)
	formula.surv <- as.formula(paste(as.character(attr(tf, "variables")[2]), "~", paste(vars_hazard_wo_time, collapse = "+")))

	# The cutoffs to compute the ROC curve	
	if(is.null(cutoffs)) {
		cutoffs <- sort(unique(data[,marker]))
	}

	# TPR and FPR (we assume that only a time is considered)
	TPR <- FPR <- matrix(ncol = nrow(cov.values), nrow = length(cutoffs))
	
	if(is.null(fitted.models)) {
		# Hazard model
		if(sum(data[[delta_var]] == 0)) {
			max.time <- max(unique(c(data[[time_var]][data[[delta_var]] == 1], max(data[[time_var]][data[[delta_var]] == 0]))))
		} else {
			max.time <- max(unique(c(data[[time_var]][data[[delta_var]] == 1])))
		}		

		# Obtain the augmented data using package pammtools 
		if(is.null(n.timepoints)) {
			ped_data <- data %>% as_ped(formula.surv, max_time = max.time)
		} else {
			ped_data <- data %>% as_ped(formula.surv, cut = seq(0, max.time, length = n.timepoints))
		}
		ped_data[,time_var] <- ped_data$tend
		
		if(length(mgcv::interpret.gam(formula.hazard)$smooth.spec) == 0) {
			discrete <- FALSE
		} else {
			discrete <- TRUE
		}
		formula.hazard.fake <- as.formula(paste0("ped_status ~ ", paste0(attr(terms.formula(formula.hazard), "term.labels"), collapse = "+")))
		fit.hazard <- bam(formula = formula.hazard.fake,
				data = ped_data, family = poisson(link = "log"), offset = ped_data$offset, discrete = discrete, select = select)

		# Location/scale model
		fit.biomarker <- loc_scale_model_mgcv(formula.biomarker, data = data)
	} else {
		fit.hazard <- fitted.models$hazard
		fit.biomarker   <- fitted.models$biomarker
	}

	# Location/scale model
	# Obtain the residuals	
		residuals <- fit.biomarker$std_residuals
	
	# Predict on the new dataframe
		p.gam <- predict_loc_scale_model_mgcv(fit.biomarker, newdata = cov.values)
	
	# Obtain \sigma(X)*residuals + \mu(X)
	# as many columns as covariate values
	# as many rows as sample size
	p.roc <- outer(residuals, p.gam$p_sd, '*')
	p.roc <- sweep(p.roc, 2, p.gam$p_mean, FUN = "+")

	p.coxph <- do.call("cbind", sapply(1:ncol(p.roc), function(i, fit, data, data.roc, predict.time) {
		# Covariate value at which to obtain predictions
		df <- data[i,,drop = FALSE]
		# Replicate the former for each individual in the dataset
		df <- do.call("rbind", replicate(length(data.roc[,i]), df, simplify = FALSE))
		# Add the "marker" value (here residuals + \alpha^T X)
		df[,marker] <- data.roc[,i]
		# Replicate df as many times as times used to obtain the integral
		df <- do.call("rbind", replicate(101, df, simplify = FALSE))
		# Add times used to obtain the integral
		df[,time_var] <- rep(seq(0, predict.time, l = 101), each = length(data.roc[,i]))
		p.cox <- matrix(predict(fit, newdata = df, type = "response"), ncol = 101, byrow = FALSE)
		list(exp(-apply(p.cox, 1, simpson, set.int = seq(0, predict.time, l = 101))))
	}, fit = fit.hazard, data = cov.values, data.roc = p.roc, predict.time = predict.time))

	den.fpr <- apply(p.coxph, 2, sum)
	den.tpr <- length(residuals) - den.fpr
	
	num.fpr <- matrix(c(sapply(1:length(cutoffs), function(i){
		apply((p.roc >= cutoffs[i])*p.coxph, 2, sum)
	})), ncol = nrow(cov.values), byrow = TRUE)	

	num.tpr <- matrix(c(sapply(1:length(cutoffs), function(i){
		apply((p.roc >= cutoffs[i])*(1-p.coxph), 2, sum)
	})), ncol = nrow(cov.values), byrow = TRUE)	


	FPR <- rbind(1, t(t(num.fpr)/den.fpr), 0)
	TPR <- rbind(1, t(t(num.tpr)/den.tpr), 0)

	area <- sapply(1:nrow(cov.values), function(x, fpr.seq, FPR, TPR) {
			obtain.ROCCurve(fpr.seq, rev(FPR[,x]), rev(TPR[,x]))$AUC
		}, fpr.seq = seq(0,1, l = 101), FPR = FPR, TPR = TPR)

	res <- list(cutoffs = c(-Inf, cutoffs, Inf), 
		TPR = TPR, 
		FPR = FPR, 
		AUC = area, 
		fitted.models = list(hazard = fit.hazard, biomarker = fit.biomarker))
}