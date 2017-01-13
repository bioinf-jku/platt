# Copyright (C) 2014 Klambauer Guenter 


#' @title Platt scaling of prediction values using precalculated coefficients. 
#' 
#' @description This function performs the mapping of the prediction values
#' to the sigmoid using the coefficients determined by the function "plattScaling".
#' 
#' @param pS An object returned from the function "platt". Class must be
#' "plattScalingResult".
#' @param predictions A vector of prediction values typically from a test set.
#' A numeric vector ist expected
#' @return numeric A numeric vector containing the Platt-scaled values.
#' 
#' @author Guenter Klambauer
#' @export

predictProb <- function(pS, predictions) {
	if (class(pS)!="plattScalingResult") stop("Input must be \"plattScalingResult\"")
	A <- pS$A; B <- pS$B; mm <- pS$norm[1]; ss <- pS$norm[2]
	if (!is.numeric(A)) stop("Parameter \"A\" must be numeric! ")
	if (!is.numeric(B)) stop("Parameter \"A\" must be numeric! ")
	if (!is.numeric(predictions)) stop("\"predictions\" must be numeric! ")
	if (length(A)!=1) stop("A single numeric value for \"A\" is expected.")
	if (length(B)!=1) stop("A single numeric value for \"B\" is expected.")
	predictions <- (predictions - mm)/ss
	vec=1/(1+exp(A*predictions+B))
	return(vec)
}