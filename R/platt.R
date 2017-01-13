# Copyright (C) 2014 Klambauer Guenter 


#' @title Platt scaling
#' @description Identifies coefficients of the sigmoid for Platt scaling.
#' 
#' @param predictions A numeric vector containing the prediction values 
#' from any machine learning method.
#' @param labels The labels of the data. Should be a vector containing only
#' zeros and ones, or a logical vector
#' 
#' @return list A list containing the parameters of the sigmoid A and B, and
#' the Platt-scaled prediction values.
#' 
#' @useDynLib platt
#' @author Guenter Klambauer
#' @export

plattScaling <- function(predictions, labels){
	if (!is.numeric(predictions)) stop("\"predictions\" must be numeric! ")
	
	if (is.numeric(labels)){
		if (all(labels %in% c(0,1))){
			#ok
 		} else {
			stop("Labels must be logical or numeric containing only 0 and 1.")
		}
	} else if (is.logical(labels)){
		# ok
	} else {
		stop("Labels must be logical or numeric containing only 0 and 1.")
	}
	
	if (length(predictions) != length(labels)){
		stop("Predictions and labels must be vectors of the same length.")
	}
	
	if (var(predictions)==0) stop("Variance of predictions is zero.")
	
	if (any(is.na(predictions)) | any(is.na(labels))){
		message("Detected NAs in predictions or labels. Removing.")
		predictions <- predictions[which(!is.na(predictions) & !is.na(labels))]
		labels <- labels[which(!is.na(predictions) & !is.na(labels))]
	}
	
	mm <- median(predictions)
	ss <- sd(predictions)
	predictions <- (predictions - mm)/ss
	
	pS=.Call("plattScaling", predictions, as.logical(labels), as.numeric(sum(!as.logical(labels))), as.numeric(sum(as.logical(labels))))
	
	A <- pS$A
	B <- pS$B
	success <- TRUE
		
	if (A > (-1e-7)) {
		warning("Curve-fitting was not successful")
		A <- (-1)
		B <- 0
		success <- FALSE
	}
	
	#newPred <- predictProb(A,B,predictions)
	newPred <- 1/(1+exp(A*predictions+B))
	
	if (sd(newPred)==0 | sd(predictions)==0){
		warning("Standard deviation of predictions is zero.")
		A <- (-1)
		B <- 0
		success <- FALSE
	} else {
		T <- try(K <- cor(newPred,predictions,method="spearman",use="pairwise.complete.obs"))
		if (inherits(T,"try-error")){
			A <- (-1)
			B <- 0
			success <- FALSE
		} else {
			if (K < 0.99){
				warning("Platt Scaling changed the ranking of the values.")
				success <- FALSE
			}
		} 
		
	}
	
	plattScalingResult <- list(pred=newPred,A=A,B=B,norm=c(mm,ss),success=success)
	class(plattScalingResult) <- "plattScalingResult"
	
	return(plattScalingResult)	
	
}

