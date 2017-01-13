# Copyright (C) 2014 Klambauer Guenter 


#' A probabilistic method to make ensemble predictions.
#' 
#' @param df A data frame containing the probabilisitic predictions of 
#' different methods. Columns are assumed to be methods.
#' @param class1Prob A numeric value that gives the prior probability of
#' observing class 1. Does not change the order of the output.
#' 
#' @return numeric A vector that contains the ensemble predictions 
#' 
#' @author Guenter Klambauer
#' @export

ensemble <- function(df,class1Prob=0.5){
	g <- class1Prob
	rr <- apply(df,1,function(x){
				pp <- prod(x,na.rm=TRUE)
				n <- length(which(!is.na(x)))
				return(  pp/( pp + prod(1-x,na.rm=TRUE)*(g/(1-g))^(n-1))  )
			})
	return(rr)
}


