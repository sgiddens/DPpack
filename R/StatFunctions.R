#' Differentially Private Mean
#'
#' This function computes the differentially private mean(s) of a given dataset
#' at user-specified levels of epsilon and delta. If the given dataset is a
#' matrix or data frame, differentially private means are computed over columns
#' and collectively satisfy differential privacy at the specified level.
#'
#' @param x Numeric vector, matrix, or data frame. Means taken over columns
#'   (when applicable).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bounds Numeric vector of lower bounds on each column of x. The
#'   length of lower.bounds must match the number of columns of x (length 1 if x
#'   is a vector). If not given, it is computed from the data to be the min
#'   value of each column of x (or the min of x if x is a vector). Note this
#'   computation may result in additional privacy loss.
#' @param upper.bounds Numeric vector of upper bounds on each column of x. The
#'   length of upper.bounds must match the number of columns of x (length 1 if x
#'   is a vector). If not given, it is computed from the data to be the max
#'   value of each column of x (or the max of x if x is a vector). Note this
#'   computation may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta) to the statistics. For example, if this function is run
#'   on a two-column matrix and alloc.proportions = c(.75, .25), then 75% of the
#'   privacy budget eps (and delt) is allocated to the statistical computation
#'   for column 1, and the remaining 25% is allocated to the statistical
#'   computation for column 2. This ensures (eps, delt)-level privacy across all
#'   computations. By default, it distributes eps and delt evenly among the
#'   calculations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return List of bounded and/or unbounded sanitized means.
#' @examples
#' meanDP(c(1,4,-2,8,-6),1,lower.bounds=-10,upper.bounds=10)
#' meanDP(c(1,4,-2,8,-6),1,which.sensitivity='unbounded',
#'   lower.bounds=-10,upper.bounds=10,mechanism='gaussian',
#'   delt=0.5,type.DP='aDP')
#' meanDP(matrix(c(1,4,-2,8,-6,0),ncol=2),1,which.sensitivity='bounded',
#'   lower.bounds=c(-10,-10),upper.bounds=c(10,10),alloc.proportions=c(1,2))
#'
#' @export
meanDP <- function (x, eps, which.sensitivity='bounded', lower.bounds=NULL,
                    upper.bounds=NULL, mechanism='laplace', delt=NULL,
                    type.DP='pDP', alloc.proportions=NULL){
  #### INPUT CHECKING ####
  {if (is.null(dim(x))){
    if (is.null(upper.bounds)){
      warning(paste("Upper bound missing and will be calculated from the data.",
              "This may represent additional privacy loss."));
      upper.bounds <- max(x);
    } else{
      if (length(upper.bounds)!=1) stop("Length of upper.bounds must be 1.");
    }
    if (is.null(lower.bounds)){
      warning(paste("Lower bound missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      lower.bounds <- min(x);
    } else{
      if (length(lower.bounds)!=1) stop("Length of lower.bounds must be 1.");
    }
    x[x<lower.bounds] <- lower.bounds;
    x[x>upper.bounds] <- upper.bounds;

  } else{
    if (is.null(upper.bounds)){
      warning(paste("Upper bounds missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      upper.bounds <- apply(x,2,max);
    } else{
      if (length(upper.bounds)!=ncol(x)) stop("Length of upper.bounds must be
                                              equal to the number of columns of x.");
    }
    if (is.null(lower.bounds)){
      warning(paste("Lower bounds missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      lower.bounds <- apply(x,2,min);
    } else{
      if (length(lower.bounds)!=ncol(x)) stop("Length of lower.bounds must be
                                              equal to the number of columns of x.");
    }
    for (i in 1:length(upper.bounds)){
      x[x[,i]<lower.bounds[i]] <- lower.bounds[i];
      x[x[,i]>upper.bounds[i]] <- upper.bounds[i];
    }
  }
  if (mechanism=='gaussian'){
    if (is.null(delt)){
      print("Must specify delta for Gaussian mechanism.");
    }
  } else if (mechanism!='laplace'){
    stop("Mechanism must be one of {'laplace', 'gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- meanDataAccess(x, lower.bounds=lower.bounds,
                               upper.bounds=upper.bounds);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  if (mechanism=='laplace'){
    sanitized.means <- laplaceMechanism(tv,eps,bs,us,which.sensitivity,
                                        alloc.proportions);
  } else if (mechanism=='gaussian'){
    sanitized.means <- gaussianMechanism(tv,eps,delt,bs,us,which.sensitivity,
                                         type.DP,alloc.proportions);
  }
  ##########

  ########## Postprocessing layer
  class(sanitized.means)<-"Sanitized Mean";
  return(sanitized.means);
  ##########
}

#' Differentially Private Variance
#'
#' This function computes the differentially private variance(s) of a given
#' dataset at user-specified levels of epsilon and delta. If the given dataset
#' is a matrix or data frame, differentially private variances are computed over
#' columns and collectively satisfy differential privacy at the specified level.
#'
#' @param x Numeric vector, matrix, or data frame. Variances taken over columns
#'   (when applicable).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bounds Numeric vector of lower bounds on each column of x. The
#'   length of lower.bounds must match the number of columns of x (length 1 if x
#'   is a vector). If not given, it is computed from the data to be the min
#'   value of each column of x (or the min of x if x is a vector). Note this
#'   computation may result in additional privacy loss.
#' @param upper.bounds Numeric vector of upper bounds on each column of x. The
#'   length of upper.bounds must match the number of columns of x (length 1 if x
#'   is a vector). If not given, it is computed from the data to be the max
#'   value of each column of x (or the max of x if x is a vector). Note this
#'   computation may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta) to the statistics. For example, if this function is run
#'   on a two-column matrix and alloc.proportions = c(.75, .25), then 75% of the
#'   privacy budget eps (and delt) is allocated to the statistical computation
#'   for column 1, and the remaining 25% is allocated to the statistical
#'   computation for column 2. This ensures (eps, delt)-level privacy across all
#'   computations. By default, it distributes eps and delt evenly among the
#'   calculations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return List of bounded and/or unbounded sanitized variances.
#' @examples
#' varDP(c(1,4,-2,8,-6),1,lower.bounds=-10,upper.bounds=10)
#' varDP(c(1,4,-2,8,-6),1,which.sensitivity='unbounded',
#'   lower.bounds=-10,upper.bounds=10,mechanism='gaussian',
#'   delt=0.5,type.DP='aDP')
#' varDP(matrix(c(1,4,-2,8,-6,0),ncol=2),1,which.sensitivity='bounded',
#'   lower.bounds=c(-10,-10),upper.bounds=c(10,10),alloc.proportions=c(1,2))
#'
#' @export
varDP <- function (x, eps, which.sensitivity='bounded', lower.bounds=NULL,
                   upper.bounds=NULL, mechanism='laplace', delt=NULL,
                   type.DP='pDP', alloc.proportions=NULL){
  #### INPUT CHECKING ####
  {if (is.null(dim(x))){
    if (is.null(upper.bounds)){
      warning(paste("Upper bound missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      upper.bounds <- max(x);
    } else{
      if (length(upper.bounds)!=1) stop("Length of upper.bounds must be 1.");
    }
    if (is.null(lower.bounds)){
      warning(paste("Lower bound missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      lower.bounds <- min(x);
    } else{
      if (length(lower.bounds)!=1) stop("Length of lower.bounds must be 1.");
    }
    x[x<lower.bounds] <- lower.bounds;
    x[x>upper.bounds] <- upper.bounds;

  } else{
    if (is.null(upper.bounds)){
      warning(paste("Upper bounds missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      upper.bounds <- apply(x,2,max);
    } else{
      if (length(upper.bounds)!=ncol(x)) stop("Length of upper.bounds must be
                                              equal to the number of columns of x.");
    }
    if (is.null(lower.bounds)){
      warning(paste("Lower bounds missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      lower.bounds <- apply(x,2,min);
    } else{
      if (length(lower.bounds)!=ncol(x)) stop("Length of lower.bounds must be
                                              equal to the number of columns of x.");
    }
    for (i in 1:length(upper.bounds)){
      x[x[,i]<lower.bounds[i]] <- lower.bounds[i];
      x[x[,i]>upper.bounds[i]] <- upper.bounds[i];
    }
  }
    if (mechanism=='gaussian'){
      if (is.null(delt)){
        print("Must specify delta for Gaussian mechanism.");
      }
    } else if (mechanism!='laplace'){
      stop("Mechanism must be one of {'laplace', 'gaussian'}.");
    }
  }
  ##########

  ########## Data access layer
  results <- varDataAccess(x, lower.bounds=lower.bounds,
                            upper.bounds=upper.bounds);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  if (mechanism=='laplace'){
    while (TRUE){ # Make sure variance is > 0 after noise
      sanitized.vars <- laplaceMechanism(tv,eps,bs,us,which.sensitivity,
                                         alloc.proportions);
      done = TRUE;
      if (!is.null(sanitized.vars$Bounded) && any(sanitized.vars$Bounded<=0)) done=FALSE;
      if (!is.null(sanitized.vars$Unbounded) && any(sanitized.vars$Unbounded<=0)) done=FALSE;
      if (done) break;
    }
  }  else if (mechanism=='gaussian'){
    while (TRUE){
      sanitized.vars <- gaussianMechanism(tv,eps,delt,bs,us,which.sensitivity,
                                          type.DP,alloc.proportions);
      done = TRUE;
      if (!is.null(sanitized.vars$Bounded) && any(sanitized.vars$Bounded<=0)) done=FALSE;
      if (!is.null(sanitized.vars$Unbounded) && any(sanitized.vars$Unbounded<=0)) done=FALSE;
      if (done) break;
    }
  }
  ##########

  ########## Postprocessing layer
  class(sanitized.vars)<-"Sanitized Variance";
  return(sanitized.vars);
  ##########
}

#' Differentially Private Standard Deviation
#'
#' This function computes the differentially private standard deviation(s) of a
#' given dataset at user-specified levels of epsilon and delta. If the given
#' dataset is a matrix or data frame, differentially private standard deviations
#' are computed over columns and collectively satisfy differential privacy at
#' the specified level.
#'
#' @param x Numeric vector, matrix, or data frame. Standard deviations taken
#'   over columns (when applicable).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bounds Numeric vector of lower bounds on each column of x. The
#'   length of lower.bounds must match the number of columns of x (length 1 if x
#'   is a vector). If not given, it is computed from the data to be the min
#'   value of each column of x (or the min of x if x is a vector). Note this
#'   computation may result in additional privacy loss.
#' @param upper.bounds Numeric vector of upper bounds on each column of x. The
#'   length of upper.bounds must match the number of columns of x (length 1 if x
#'   is a vector). If not given, it is computed from the data to be the max
#'   value of each column of x (or the max of x if x is a vector). Note this
#'   computation may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta) to the statistics. For example, if this function is run
#'   on a two-column matrix and alloc.proportions = c(.75, .25), then 75% of the
#'   privacy budget eps (and delt) is allocated to the statistical computation
#'   for column 1, and the remaining 25% is allocated to the statistical
#'   computation for column 2. This ensures (eps, delt)-level privacy across all
#'   computations. By default, it distributes eps and delt evenly among the
#'   calculations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return List of bounded and/or unbounded sanitized standard deviations.
#' @examples
#' sdDP(c(1,4,-2,8,-6),1,lower.bounds=-10,upper.bounds=10)
#' sdDP(c(1,4,-2,8,-6),1,which.sensitivity='unbounded',
#'   lower.bounds=-10,upper.bounds=10,mechanism='gaussian',
#'   delt=0.5,type.DP='aDP')
#' sdDP(matrix(c(1,4,-2,8,-6,0),ncol=2),1,which.sensitivity='bounded',
#'   lower.bounds=c(-10,-10),upper.bounds=c(10,10),alloc.proportions=c(1,2))
#'
#' @export
sdDP <- function (x, eps, which.sensitivity='bounded', lower.bounds=NULL,
                  upper.bounds=NULL, mechanism='laplace', delt=NULL,
                  type.DP='pDP', alloc.proportions=NULL){
  ########## Input checking

  ##########

  ########## Data Access/privacy layer
  sanitized.variances <- varDP(x,eps,which.sensitivity,lower.bounds,upper.bounds,
                               mechanism,delt,type.DP,alloc.proportions);
  ##########

  ########## Postprocessing layer
  sanitized.sds = list();
  if (which.sensitivity=='bounded' || which.sensitivity=='both') {
    sanitized.sds[["Bounded"]] <- sqrt(sanitized.variances$Bounded);
  }
  if (which.sensitivity=='unbounded' || which.sensitivity=='both') {
    sanitized.sds[["Unbounded"]] <- sqrt(sanitized.variances$Unbounded);
  }
  class(sanitized.sds) <- "Sanitized Std Dev";
  return(sanitized.sds);
  ##########
}

#' Differentially Private Covariance
#'
#' This function computes the differentially private covariance of a pair of
#' vectors at user-specified levels of epsilon and delta.
#'
#' @param x1,x2 Numeric vectors.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bound1,lower.bound2 Real numbers giving the lower bounds of x1
#'   and x2, respectively. If not given, they are computed from the data to be
#'   the min value of x1 and x2, respectively. Note this computation may result
#'   in additional privacy loss.
#' @param upper.bound1,upper.bound2 Real numbers giving the upper bounds of x1
#'   and x2, respectively. If not given, they are computed from the data to be
#'   the max value of x1 and x2, respectively. Note this computation may result
#'   in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @return List of bounded and/or unbounded sanitized covariances.
#' @examples
#' covDP(c(1,4,-2,8,-6),c(1,3,2,2,4),1,which.sensitivity='bounded',
#'   lower.bound1=-10,upper.bound1=10,lower.bound2=0,upper.bound2=5,
#'   mechanism='laplace')
#' covDP(c(1,4,-2,8,-6),c(1,3,2,2,4),1,which.sensitivity='unbounded',
#'   lower.bound1=-10,upper.bound110,lower.bound2=0,upper.bound2=5,
#'   mechanism='gaussian',delt=0.5,type.DP='aDP')
#'
#' @export
covDP <- function (x1, x2, eps, which.sensitivity='bounded',
                  lower.bound1=NULL, upper.bound1=NULL,
                  lower.bound2=NULL, upper.bound2=NULL,
                  mechanism='laplace', delt=NULL, type.DP='pDP'){
  #### INPUT CHECKING ####
  {
  if (is.null(upper.bound1)){
    warning("Upper bound on x1 missing and will be calculated from the data.
            This may represent additional privacy loss.");
    upper.bound1 <- max(x1);
  } else if (length(upper.bound1)!=1) stop("Length of upper.bound1 must be 1.");
  if (is.null(lower.bound1)){
    warning("Lower bound on x1 missing and will be calculated from the data.
            This may represent additional privacy loss.");
    lower.bound1 <- min(x1);
  } else if (length(lower.bound1)!=1) stop("Length of lower.bound1 must be 1.");
  if (is.null(upper.bound2)){
    warning("Upper bound on x2 missing and will be calculated from the data.
          This may represent additional privacy loss.");
    upper.bound2 <- max(x2);
  } else if (length(upper.bound2)!=1) stop("Length of upper.bound2 must be 1.");
  if (is.null(lower.bound2)){
    warning("Lower bound on x2 missing and will be calculated from the data.
          This may represent additional privacy loss.");
    lower.bound2 <- min(x2);
  } else if (length(lower.bound2)!=1) stop("Length of lower.bound2 must be 1.");
  x1[x1<lower.bound1] <- lower.bound1;
  x1[x1>upper.bound1] <- upper.bound1;
  x2[x2<lower.bound2] <- lower.bound2;
  x2[x2>upper.bound2] <- upper.bound2;
  if (mechanism=='gaussian'){
    if (is.null(delt)){
      print("Must specify delta for Gaussian mechanism.");
    }
  } else if (mechanism!='laplace'){
    stop("Mechanism must be one of {'laplace', 'gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- covDataAccess(x1, x2, lower.bound1, upper.bound1,
                           lower.bound2, upper.bound2);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  if (mechanism=='laplace'){
    sanitized.cov <- laplaceMechanism(tv,eps,bs,us,which.sensitivity);
  } else if (mechanism=='gaussian'){
    sanitized.cov <- gaussianMechanism(tv,eps,delt,bs,us,which.sensitivity,
                                       type.DP);
  }

  ##########

  ########## Postprocessing layer
  class(sanitized.cov)<-"Sanitized Covariance";
  return(sanitized.cov);
  ##########
}

#' Differentially Private Histogram
#'
#' This function computes a differentially private histogram from a vector at
#' user-specified levels of epsilon and delta. A histogram object is returned
#' with sanitized values for the counts for easy plotting (see examples).
#'
#' @param x Numeric vector from which the histogram will be formed.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param breaks Identical to the argument with the same name from
#'   \code{\link[graphics]{hist}}.
#' @param normalize Logical value. If FALSE (default), returned histogram counts
#'   correspond to frequencies. If TRUE, returned histogram counts correspond to
#'   densities (i.e. area of histogram is one).
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bound Real number giving the lower bound of x. If not given, it
#'   is computed from the data to be the min value of x. Note this computation
#'   may result in additional privacy loss.
#' @param upper.bound Real number giving the upper bound of x. If not given, it
#'   is computed from the data to be the max value of x. Note this computation
#'   may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param allow.negative Logical value. If FALSE (default), any negative values
#'   in the sanitized histogram due to the added noise will be set to 0. If
#'   TRUE, the negative values (if any) will be returned.
#' @return List of bounded and/or unbounded sanitized histograms.
#' @examples
#' result <- histogramDP(c(1,1,-2,8,-6),1,which.sensitivity='bounded',
#'   lower.bound=-10,upper.bound=10,mechanism='laplace')
#' plot(result$Bounded)
#' result <- histogramDP(c(1,1,-2,8,-6),1,normalize=TRUE,
#'   which.sensitivity='unbounded',lower.bound=-10,upper.bound=10,
#'   mechanism='gaussian',delt=0.5,type.DP='aDP',allow.negative=FALSE)
#' plot(result$Unbounded)
#'
#' @export
histogramDP <- function(x, eps, breaks="Sturges", normalize=FALSE,
                        which.sensitivity='bounded',
                        lower.bound=NULL, upper.bound=NULL,
                        mechanism='laplace', delt=NULL, type.DP='pDP',
                        allow.negative=FALSE){
  #### INPUT CHECKING ####
  {
    if (is.null(upper.bound)){
      warning("Upper bound missing and will be calculated from the data.
            This may represent additional privacy loss.");
      upper.bound <- max(x);
    } else if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.");
    if (is.null(lower.bound)){
      warning("Lower bound missing and will be calculated from the data.
            This may represent additional privacy loss.");
      lower.bound <- min(x);
    } else if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.");
    x[x<lower.bound] <- lower.bound;
    x[x>upper.bound] <- upper.bound;
    if (mechanism=='gaussian'){
      if (is.null(delt)){
        print("Must specify delta for Gaussian mechanism.");
      }
    } else if (mechanism!='laplace'){
      stop("Mechanism must be one of {'laplace', 'gaussian'}.");
    }
  }
  ##########

  ########## Data access layer
  results <- histogramDataAccess(x, breaks);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  counts <- tv$counts;
  # Might need to verify this is right later (also see tableDP)
  # This means that each param[i] in the mechanism becomes bs/eps rather
  #       than bs[i]/(alloc.proportions[i]*eps)
  bs <- rep(bs, length(counts))/length(counts);
  us <- rep(us, length(counts))/length(counts);
  if (mechanism=='laplace'){
    sanitized.counts <- laplaceMechanism(counts,eps,bs,us,which.sensitivity);
  } else if (mechanism=='gaussian'){
    sanitized.counts <- gaussianMechanism(counts,eps,delt,bs,us,which.sensitivity,
                                          type.DP);
  }
  ##########

  ########## Postprocessing layer
  sanitized.hist <- list();
  if (!is.null(sanitized.counts$Bounded)){
    sc <- sanitized.counts$Bounded;
    if (!allow.negative) sc[sc<0] <- 0;
    if (normalize) {
      sc <- sc/sum(sc);
    } else{
      sc <- round(sc);
    }
    tv$counts <- sc;
    sanitized.hist[["Bounded"]] <- tv;
  }
  if (!is.null(sanitized.counts$Unbounded)){
    sc <- sanitized.counts$Unbounded;
    if (!allow.negative) sc[sc<0] <- 0;
    if (normalize) {
      sc <- sc/sum(sc);
    } else{
      sc <- round(sc);
    }
    tv$counts <- sc;
    sanitized.hist[["Unbounded"]] <- tv;
  }

  class(sanitized.hist)<-"Sanitized Histogram";
  return(sanitized.hist);
  ##########
}

#' Differentially Private Contingency Table
#'
#' This function computes a differentially private contingency table from two
#' vectors of data at user-specified levels of epsilon and delta.
#'
#' @param x,y Vectors of data from which to create the contingency table.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param allow.negative Logical value. If FALSE (default), any negative values
#'   in the sanitized table due to the added noise will be set to 0. If TRUE,
#'   the negative values (if any) will be returned.
#' @return List of bounded and/or unbounded sanitized contingency tables.
#' @examples
#' x <- MASS::Cars93$Type;
#' y <- MASS::Cars93$Origin;
#' tableDP(x,y,1,which.sensitivity='bounded',mechanism='laplace')
#' tableDP(x,y,1,which.sensitivity='unbounded',mechanism='gaussian',delt=0.5,
#'   type.DP='aDP',allow.negative=FALSE)
#'
#' @export
tableDP <- function(x, y, eps, which.sensitivity='bounded', mechanism='laplace',
                    delt=NULL, type.DP='pDP', allow.negative=FALSE){
  #### INPUT CHECKING ####
  {
    if (mechanism=='gaussian'){
      if (is.null(delt)){
        print("Must specify delta for Gaussian mechanism.");
      }
    } else if (mechanism!='laplace'){
      stop("Mechanism must be one of {'laplace', 'gaussian'}.");
    }
  }
  ##########

  ########## Data access layer
  results <- tableDataAccess(x, y);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;

  # Flatten tv, but keep values to recreate later
  rnames <- row.names(tv);
  cnames <- colnames(tv);
  dims <- dim(tv);
  dim(tv) <- NULL;
  ##########

  ########## Privacy layer
  # Might need to verify this is right later (also see histogramDP)
  # This means that each param[i] in the mechanism becomes bs/eps rather
  #       than bs[i]/(alloc.proportions[i]*eps)
  bs <- rep(bs, length(tv))/length(tv);
  us <- rep(us, length(tv))/length(tv);
  if (mechanism=='laplace'){
    sanitized.tables <- laplaceMechanism(tv,eps,bs,us,which.sensitivity);
  } else if (mechanism=='gaussian'){
    sanitized.tables <- gaussianMechanism(tv,eps,delt,bs,us,which.sensitivity,
                                          type.DP);
  }
  ##########

  ########## Postprocessing layer
  # Unflatten and round tables
  sanitized.table <- list();
  if (!is.null(sanitized.tables$Bounded)){
    bounded.table <- sanitized.tables$Bounded;
    dim(bounded.table) <- dims;
    row.names(bounded.table) <- rnames;
    colnames(bounded.table) <- cnames;
    bounded.table <- round(bounded.table);
    if (!allow.negative) bounded.table[bounded.table<0] <- 0;
    sanitized.table[["Bounded"]] <- bounded.table;
  }
  if (!is.null(sanitized.tables$Unbounded)){
    unbounded.table <- sanitized.tables$Unbounded;
    dim(unbounded.table) <- dims;
    row.names(unbounded.table) <- rnames;
    colnames(unbounded.table) <- cnames;
    unbounded.table <- round(unbounded.table);
    if (!allow.negative) unbounded.table[unbounded.table<0] <- 0;
    sanitized.table[["Unbounded"]] <- unbounded.table;
  }

  class(sanitized.table)<-"Sanitized Contingency Table";
  return(sanitized.table);
  ##########
}

#' Differentially Private Pooled Variance
#'
#' This function computes the differentially private pooled variance from two or
#' more vectors of data at user-specified levels of epsilon and delta.
#'
#' @param ... Two or more vectors from which to compute the pooled variance.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bound Real number giving the lower bound of the input data. If
#'   not given, it is computed from the data to be the min value over all given
#'   data. Note this computation may result in additional privacy loss.
#' @param upper.bound Real number giving the upper bound of the input data. If
#'   not given, it is computed from the data to be the max value over all given
#'   data. Note this computation may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param approx.n.max Logical indicating whether to approximate n.max, which is
#'   defined to be the length of the largest input vector. Approximation is best
#'   if n.max is very large.
#' @return List of bounded and/or unbounded sanitized pooled variances.
#' @examples
#' pooledVarDP(c(1,4,-2,8,-6),c(1,2),c(-5,-7),eps=1,which.sensitivity='bounded',
#'   lower.bound=-10,upper.bound=10,mechanism='laplace')
#' pooledVarDP(c(1,4,-2,8,-6),c(1,2),c(-5,-7),eps=1,
#'   which.sensitivity='unbounded',lower.bound=-10,upper.bound=10,
#'   mechanism='gaussian',delt=0.5,type.DP='aDP',approx.n.max=TRUE)
#'
#' @export
pooledVarDP <- function(..., eps=1, which.sensitivity='bounded',
                        lower.bound=NULL, upper.bound=NULL,
                        mechanism='laplace', delt=NULL, type.DP='pDP',
                        approx.n.max=FALSE){
  samples <- list(...);
  #### INPUT CHECKING ####
  {
  J = length(samples);
  if (is.null(lower.bound)){
    warning(paste("Lower bound missing and will be calculated from the data.",
                  "This may represent additional privacy loss."));
    lower.bounds <- numeric(J);
    for (j in 1:J){
      lower.bounds[j] <- min(samples[[j]])
    }
    lower.bound <- min(lower.bounds);
  } else if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.");
  if (is.null(upper.bound)){
    warning(paste("Upper bound missing and will be calculated from the data.",
                  "This may represent additional privacy loss."));
    upper.bounds <- numeric(J);
    for (j in 1:J){
      upper.bounds[j] <- max(samples[[j]])
    }
    upper.bound <- max(upper.bounds);
  } else if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.");

  for (j in 1:J){
    samples[[j]][samples[[j]]<lower.bound] <- lower.bound;
    samples[[j]][samples[[j]]>upper.bound] <- upper.bound;
    # if (any(samples[[j]]<lower.bound || any(samples[[j]]>upper.bound))){
    #   stop("Each element in samples must be contained in (lower.bounds, upper.bounds).")
    # }
  }

  if (mechanism=='gaussian'){
    if (is.null(delt)){
      print("Must specify delta for Gaussian mechanism.");
    }
  } else if (mechanism!='laplace'){
    stop("Mechanism must be one of {'laplace', 'gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- pooledVarDataAccess(samples, lower.bound=lower.bound,
                                 upper.bound=upper.bound,
                                 approx.n.max=approx.n.max);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  if (mechanism=='laplace'){
    while (TRUE){ # Make sure variance is > 0 after noise
      sanitized.vars <- laplaceMechanism(tv,eps,bs,us,which.sensitivity);
      done = TRUE;
      if (!is.null(sanitized.vars$Bounded) && any(sanitized.vars$Bounded<=0)) done=FALSE;
      if (!is.null(sanitized.vars$Unbounded) && any(sanitized.vars$Unbounded<=0)) done=FALSE;
      if (done) break;
    }
  } else if (mechanism=='gaussian'){
    while (TRUE){ # Make sure variance is > 0 after noise
      sanitized.vars <- gaussianMechanism(tv,eps,delt,bs,us,which.sensitivity,
                                          type.DP);
      done = TRUE;
      if (!is.null(sanitized.vars$Bounded) && any(sanitized.vars$Bounded<=0)) done=FALSE;
      if (!is.null(sanitized.vars$Unbounded) && any(sanitized.vars$Unbounded<=0)) done=FALSE;
      if (done) break;
    }
  }
  ##########

  ########## Postprocessing layer
  class(sanitized.vars)<-"Sanitized Pooled Variance";
  return(sanitized.vars);
  ##########
}

#' Differentially Private Pooled Covariance
#'
#' This function computes the differentially private pooled covariance from two
#' or more two-column matrices of data at user-specified levels of epsilon and
#' delta.
#'
#' @param ... Two or more matrices, each with two columns from which to compute
#'   the pooled covariance.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bound1,lower.bound2 Real numbers giving the lower bounds over
#'   the first and second columns of all input data, respectively. If not given,
#'   they are computed from the data to be the min value over all first and
#'   second columns of the given data, respectively. Note this computation may
#'   result in additional privacy loss.
#' @param upper.bound1,upper.bound2 Real numbers giving the upper bounds over
#'   the first and second columns of all input data, respectively. If not given,
#'   they are computed from the data to be the max value over all first and
#'   second columns of the given data, respectively. Note this computation may
#'   result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'laplace',
#'   'gaussian'}. See \code{\link{laplaceMechanism}} and
#'   \code{\link{gaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter
#'   (necessary if using Gaussian mechanism).
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). See \code{\link{gaussianMechanism}} for a more
#'   detailed description of this parameter. Only used if mechanism is
#'   'gaussian'.
#' @param approx.n.max Logical indicating whether to approximate n.max, which is
#'   defined to be the length of the largest input vector. Approximation is best
#'   if n.max is very large.
#' @return List of bounded and/or unbounded sanitized pooled covariances.
#' @examples
#' x1 <- matrix(c(1,4,-2,8,-6,-3),ncol=2)
#' x2 <- matrix(c(1,2,-5,7),ncol=2)
#' pooledVarDP(x1,x2,eps=1,which.sensitivity='bounded',
#'   lower.bound1=-10,upper.bound1=10,lower.bound2=-10,upper.bound2=10,
#'   mechanism='laplace')
#' pooledVarDP(x1,x2,eps=1,which.sensitivity='unbounded',
#'   lower.bound1=-10,upper.bound1=10,lower.bound2=-10,upper.bound2=10,
#'   mechanism='gaussian',delt=0.5,type.DP='aDP',approx.n.max=TRUE)
#'
#' @export
pooledCovDP <- function(..., eps=1, which.sensitivity='bounded',
                        lower.bound1=NULL, upper.bound1=NULL,
                        lower.bound2=NULL, upper.bound2=NULL,
                        mechanism='laplace', delt=NULL, type.DP='pDP',
                        approx.n.max=FALSE){
  samples <- list(...);
  #### INPUT CHECKING ####
  {
  J = length(samples);
  if (is.null(lower.bound1)){
    warning("Lower bound for first column missing and will be calculated from the data.
            This may represent additional privacy loss.");
    lower.bounds1 <- numeric(J);
    for (j in 1:J){
      lower.bounds1[j] <- min(samples[[j]][,1]);
    }
    lower.bound1 <- min(lower.bounds1);
  } else if (length(lower.bound1)!=1) stop("Length of lower.bound1 must be 1.");
  if (is.null(lower.bound2)){
    warning("Lower bound for second column missing and will be calculated from the data.
            This may represent additional privacy loss.");
    lower.bounds2 <- numeric(J);
    for (j in 1:J){
      lower.bounds2[j] <- min(samples[[j]][,2]);
    }
    lower.bound2 <- min(lower.bounds2);
  } else if (length(lower.bound2)!=1) stop("Length of lower.bound2 must be 1.");
  if (is.null(upper.bound1)){
    warning("Upper bound for first column missing and will be calculated from the data.
            This may represent additional privacy loss.");
    upper.bounds1 <- numeric(J);
    for (j in 1:J){
      upper.bounds1[j] <- max(samples[[j]][,1]);
    }
    upper.bound1 <- max(upper.bounds1);
  } else if (length(upper.bound1)!=1) stop("Length of upper.bound1 must be 1.");
  if (is.null(upper.bound2)){
    warning("Upper bound for second column missing and will be calculated from the data.
            This may represent additional privacy loss.");
    upper.bounds2 <- numeric(J);
    for (j in 1:J){
      upper.bounds2[j] <- max(samples[[j]][,2]);
    }
    upper.bound2 <- max(upper.bounds2);
  } else if (length(upper.bound2)!=1) stop("Length of upper.bound2 must be 1.");

  for (j in 1:J){
    samples[[j]][samples[[j]][,1]<lower.bound1,1] <- lower.bound1;
    samples[[j]][samples[[j]][,1]>upper.bound1,1] <- upper.bound1;
    samples[[j]][samples[[j]][,2]<lower.bound2,2] <- lower.bound2;
    samples[[j]][samples[[j]][,2]>upper.bound2,2] <- upper.bound2;
    # if (any(samples[[j]][,1]<lower.bound1 || any(samples[[j]][,1]>upper.bound1))){
    #   stop("Each element in first column of samples must be contained in
    #        (lower.bound1, upper.bound1).")
    # }
  }
  # for (j in 1:J){
  #   if (any(samples[[j]][,2]<lower.bound2 || any(samples[[j]][,2]>upper.bound2))){
  #     stop("Each element in second column of samples must be contained in
  #          (lower.bound2, upper.bound2).")
  #   }
  # }

  if (mechanism=='gaussian'){
    if (is.null(delt)){
      print("Must specify delta for Gaussian mechanism.");
    }
  } else if (mechanism!='laplace'){
    stop("Mechanism must be one of {'laplace', 'gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- pooledCovDataAccess(samples, lower.bound1=lower.bound1,
                                 upper.bound1=upper.bound1,
                                 lower.bound2=lower.bound2,
                                 upper.bound2=upper.bound2,
                                 approx.n.max=approx.n.max);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  if (mechanism=='laplace'){
    sanitized.cov <- laplaceMechanism(tv,eps,bs,us,which.sensitivity);
  } else if (mechanism=='gaussian'){
    sanitized.cov <- gaussianMechanism(tv,eps,delt,bs,us,which.sensitivity,
                                       type.DP);
  }
  ##########

  ########## Postprocessing layer
  class(sanitized.cov)<-"Sanitized Pooled Covariance";
  return(sanitized.cov);
  ##########
}

#' Differentially Private Quantile
#'
#' This function computes the differentially private quantile of an input
#' vector at user-specified levels of epsilon and delta.
#'
#' @param x Numeric vector of which the quantile will be taken.
#' @param quant Real number between 0 and 1 indicating which quantile to return.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bound Real number giving the lower bound of x. If not given, it
#'   is computed from the data to be the min value of x. Note this computation
#'   may result in additional privacy loss.
#' @param upper.bound Real number giving the upper bound of x. If not given, it
#'   is computed from the data to be the max value of x. Note this computation
#'   may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'exponential'}.
#'   See \code{\link{exponentialMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter (not
#'   currently used).
#' @return List of bounded and/or unbounded sanitized quantiles.
#' @examples
#' quantileDP(c(1,1,-2,8,-6),.25,1,which.sensitivity='bounded',
#'   lower.bound=-10,upper.bound=10,mechanism='exponential')
#' quantileDP(c(1,1,-2,8,-6),.75,1,which.sensitivity='unbounded',
#'   lower.bound=-10,upper.bound=10,mechanism='exponential')
#'
#' @export
quantileDP <- function (x, quant, eps, which.sensitivity='bounded',
                        lower.bound=NULL, upper.bound=NULL,
                        mechanism='exponential', delt=NULL){
  # NOTE: See
  # https://github.com/IBM/differential-privacy-library/blob/main/diffprivlib/tools/quantiles.py
  #
  # NOTE: The data access and privacy layers are somewhat mixed.

  #### INPUT CHECKING ####
  {
  if (is.null(upper.bound)){
    warning(paste("Upper bound missing and will be calculated from the data.",
                  "This may represent additional privacy loss."));
    upper.bound <- max(x);
  } else{
    if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.");
  }
  if (is.null(lower.bound)){
    warning(paste("Lower bound missing and will be calculated from the data.",
                  "This may represent additional privacy loss."));
    lower.bound <- min(x);
  } else{
    if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.");
  }
  x[x<lower.bound] <- lower.bound;
  x[x>upper.bound] <- upper.bound;
  if (quant<0 || quant>1) stop("quant must be between 0 and 1.")
  if (mechanism!='exponential'){
    stop("Mechanism must be one of {'exponential'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- quantileDataAccess(x, quant, lower.bound=lower.bound,
                            upper.bound=upper.bound);
  utility <- results$Utility;
  sorted <- results$Sorted;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;

  ##########

  ########## Privacy layer
  if (mechanism=='exponential'){
    sanitized.indices <- exponentialMechanism(utility, eps, bs, us,
                                              which.sensitivity,
                                              measure=diff(sorted));
    sanitized.quantile <- list();
    if (!is.null(sanitized.indices$Bounded)){
      bounded.idx <- sanitized.indices$Bounded;
      bounded.sanitized <- runif(1)*(sorted[bounded.idx+1]-sorted[bounded.idx]) +
        sorted[bounded.idx];
      sanitized.quantile[["Bounded"]] <- bounded.sanitized;
    }
    if (!is.null(sanitized.indices$Unbounded)){
      unbounded.idx <- sanitized.indices$Unbounded;
      unbounded.sanitized <- runif(1)*(sorted[unbounded.idx+1]-sorted[unbounded.idx]) +
        sorted[unbounded.idx];
      sanitized.quantile[["Unbounded"]] <- unbounded.sanitized;
    }
  }
  ##########

  ########## Postprocessing layer
  class(sanitized.quantile)<-"Sanitized Quantile";
  return(sanitized.quantile);
  ##########
}

#' Differentially Private Median
#'
#' This function computes the differentially private median of an input
#' vector at user-specified levels of epsilon and delta.
#'
#' @param x Numeric vector of which the median will be taken.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param lower.bound Real number giving the lower bound of x. If not given, it
#'   is computed from the data to be the min value of x. Note this computation
#'   may result in additional privacy loss.
#' @param upper.bound Real number giving the upper bound of x. If not given, it
#'   is computed from the data to be the max value of x. Note this computation
#'   may result in additional privacy loss.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'exponential'}.
#'   See \code{\link{exponentialMechanism}} for a description of the supported
#'   mechanisms.
#' @param delt Positive real number defining the delta privacy parameter (not
#'   currently used).
#' @return List of bounded and/or unbounded sanitized medians.
#' @examples
#' medianDP(c(1,1,-2,8,-6),1,which.sensitivity='bounded',
#'   lower.bound=-10,upper.bound=10,mechanism='exponential')
#' medianDP(c(1,1,-2,8,-6),1,which.sensitivity='unbounded',
#'   lower.bound=-10,upper.bound=10,mechanism='exponential')
#'
#' @export
medianDP <- function (x, eps, which.sensitivity='bounded',
                      lower.bound=NULL, upper.bound=NULL,
                      mechanism='exponential', delt=NULL){
  sanitized.median <- quantileDP(x,.5,eps,which.sensitivity,lower.bound,
                                 upper.bound,mechanism,delt);
  class(sanitized.median)<-"Sanitized Median";
  return(sanitized.median);
  ##########
}





