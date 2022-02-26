#' Differentially Private Mean
#'
#' This function computes the differentially private mean(s) of a given dataset
#' at user-specified privacy levels of epsilon and delta. If the given dataset
#' is a matrix or data frame, differentially private means are computed over
#' columns and collectively satisfy differential privacy at the specified level.
#'
#' @param x Numeric vector, matrix, or data frame. Means taken over columns
#'   (when applicable).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bounds Numeric vector of global or public lower bounds on each
#'   column of x. The length of lower.bounds must match the number of columns of
#'   x (length 1 if x is a vector).
#' @param upper.bounds Numeric vector of global or public upper bounds on each
#'   column of x. The length of upper.bounds must match the number of columns of
#'   x (length 1 if x is a vector).
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta) to the statistics. For example, if this function is run
#'   on a two-column matrix and alloc.proportions = c(.75, .25), then 75% of the
#'   privacy budget eps (and delta) is allocated to the statistical computation
#'   for column 1, and the remaining 25% is allocated to the statistical
#'   computation for column 2. This ensures (eps, delta)-level privacy across
#'   all computations. By default, it distributes eps and delta evenly among the
#'   calculations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return Sanitized mean(s) based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' D <- rnorm(500, mean=3, sd=2)
#' lb <- -3 # 3 std devs below mean
#' ub <- 9 # 3 std devs above mean
#' meanDP(D, 1, lb, ub)
#' meanDP(D,.5, lb, ub, which.sensitivity='unbounded', mechanism='Gaussian',
#'   delta=0.01)
#' D.2col <- matrix(D,ncol=2)
#' meanDP(D.2col, 1, lower.bounds=c(lb,lb), upper.bounds=c(ub,ub),
#'   which.sensitivity='bounded', type.DP='pDP', alloc.proportions=c(1,2))
#'
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#' @export
meanDP <- function (x, eps, lower.bounds, upper.bounds,
                    which.sensitivity='bounded', mechanism='Laplace',
                    delta=0, type.DP='aDP', alloc.proportions=NULL){
  #### INPUT CHECKING ####
  {if (is.null(dim(x))){
    if (length(upper.bounds)!=1) stop("Length of upper.bounds must be 1.");
    if (length(lower.bounds)!=1) stop("Length of lower.bounds must be 1.");

    x[x<lower.bounds] <- lower.bounds;
    x[x>upper.bounds] <- upper.bounds;
  } else{
    if (length(upper.bounds)!=ncol(x)) {
      stop("Length of upper.bounds must be equal to the number of columns of x.");
    }
    if (length(lower.bounds)!=ncol(x)) {
      stop("Length of lower.bounds must be equal to the number of columns of x.");
    }
    for (i in 1:length(upper.bounds)){
      x[x[,i]<lower.bounds[i]] <- lower.bounds[i];
      x[x[,i]>upper.bounds[i]] <- upper.bounds[i];
    }
  }
  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
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
  if (mechanism=='Laplace'){
    sanitized.means <- LaplaceMechanism(tv,eps,bs,us,which.sensitivity,
                                        alloc.proportions);
  } else if (mechanism=='Gaussian'){
    sanitized.means <- GaussianMechanism(tv,eps,delta,bs,us,which.sensitivity,
                                         type.DP,alloc.proportions);
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both') class(sanitized.means)<-"Sanitized Mean";
  return(sanitized.means);
  ##########
}

#' Differentially Private Variance
#'
#' This function computes the differentially private variance(s) of a given
#' dataset at user-specified privacy levels of epsilon and delta. If the given
#' dataset is a matrix or data frame, differentially private variances are
#' computed over columns and collectively satisfy differential privacy at the
#' specified level.
#'
#' @param x Numeric vector, matrix, or data frame. Variances taken over columns
#'   (when applicable).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bounds Numeric vector of global or public lower bounds on each
#'   column of x. The length of lower.bounds must match the number of columns of
#'   x (length 1 if x is a vector).
#' @param upper.bounds Numeric vector of global or public upper bounds on each
#'   column of x. The length of upper.bounds must match the number of columns of
#'   x (length 1 if x is a vector).
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta) to the statistics. For example, if this function is run
#'   on a two-column matrix and alloc.proportions = c(.75, .25), then 75% of the
#'   privacy budget eps (and delta) is allocated to the statistical computation
#'   for column 1, and the remaining 25% is allocated to the statistical
#'   computation for column 2. This ensures (eps, delta)-level privacy across
#'   all computations. By default, it distributes eps and delta evenly among the
#'   calculations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return Sanitized variance(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' D <- rnorm(500, mean=3, sd=2)
#' lb <- -3 # 3 std devs below mean
#' ub <- 9 # 3 std devs above mean
#' varDP(D, 1, lb, ub)
#' varDP(D,.5, lb, ub, which.sensitivity='unbounded', mechanism='Gaussian',
#'   delta=0.01)
#' D.2col <- matrix(D,ncol=2)
#' varDP(D.2col, 1, lower.bounds=c(lb,lb), upper.bounds=c(ub,ub),
#'   which.sensitivity='bounded', type.DP='pDP', alloc.proportions=c(1,2))
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
varDP <- function (x, eps, lower.bounds, upper.bounds,
                   which.sensitivity='bounded', mechanism='Laplace', delta=0,
                   type.DP='aDP', alloc.proportions=NULL){
  #### INPUT CHECKING ####
  {if (is.null(dim(x))){
    if (length(upper.bounds)!=1) stop("Length of upper.bounds must be 1.");
    if (length(lower.bounds)!=1) stop("Length of lower.bounds must be 1.");
    x[x<lower.bounds] <- lower.bounds;
    x[x>upper.bounds] <- upper.bounds;
  } else{
    if (length(upper.bounds)!=ncol(x)) {
      stop("Length of upper.bounds must be equal to the number of columns of x.");
    }
    if (length(lower.bounds)!=ncol(x)) {
      stop("Length of lower.bounds must be equal to the number of columns of x.");
    }
    for (i in 1:length(upper.bounds)){
      x[x[,i]<lower.bounds[i]] <- lower.bounds[i];
      x[x[,i]>upper.bounds[i]] <- upper.bounds[i];
    }
  }
  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
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
  if (mechanism=='Laplace'){
    while (TRUE){ # Make sure variance is > 0 after noise
      sanitized.vars <- LaplaceMechanism(tv,eps,bs,us,which.sensitivity,
                                         alloc.proportions)
      done <- TRUE
      if (which.sensitivity=='both'){
        if (any(sanitized.vars$Bounded<=0)) done <- FALSE
        if (any(sanitized.vars$Unbounded<=0)) done <- FALSE
      } else if (any(sanitized.vars<=0)) done <- FALSE
      if (done) break;
    }
  }  else if (mechanism=='Gaussian'){
    while (TRUE){
      sanitized.vars <- GaussianMechanism(tv,eps,delta,bs,us,which.sensitivity,
                                          type.DP,alloc.proportions);
      done <- TRUE
      if (which.sensitivity=='both'){
        if (any(sanitized.vars$Bounded<=0)) done <- FALSE
        if (any(sanitized.vars$Unbounded<=0)) done <- FALSE
      } else if (any(sanitized.vars<=0)) done <- FALSE
      if (done) break;
    }
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both') class(sanitized.vars)<-"Sanitized Variance";
  return(sanitized.vars);
  ##########
}

#' Differentially Private Standard Deviation
#'
#' This function computes the differentially private standard deviation(s) of a
#' given dataset at user-specified privacy levels of epsilon and delta. If the
#' given dataset is a matrix or data frame, differentially private standard
#' deviations are computed over columns and collectively satisfy differential
#' privacy at the specified level.
#'
#' @param x Numeric vector, matrix, or data frame. Standard deviations taken
#'   over columns (when applicable).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bounds Numeric vector of global or public lower bounds on each
#'   column of x. The length of lower.bounds must match the number of columns of
#'   x (length 1 if x is a vector).
#' @param upper.bounds Numeric vector of global or public upper bounds on each
#'   column of x. The length of upper.bounds must match the number of columns of
#'   x (length 1 if x is a vector).
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta) to the statistics. For example, if this function is run
#'   on a two-column matrix and alloc.proportions = c(.75, .25), then 75% of the
#'   privacy budget eps (and delta) is allocated to the statistical computation
#'   for column 1, and the remaining 25% is allocated to the statistical
#'   computation for column 2. This ensures (eps, delta)-level privacy across
#'   all computations. By default, it distributes eps and delta evenly among the
#'   calculations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return Sanitized standard deviation(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' D <- rnorm(500, mean=3, sd=2)
#' lb <- -3 # 3 std devs below mean
#' ub <- 9 # 3 std devs above mean
#' sdDP(D, 1, lb, ub)
#' sdDP(D,.5, lb, ub, which.sensitivity='unbounded', mechanism='Gaussian',
#'   delta=0.01)
#' D.2col <- matrix(D,ncol=2)
#' sdDP(D.2col, 1, lower.bounds=c(lb,lb), upper.bounds=c(ub,ub),
#'   which.sensitivity='bounded', type.DP='pDP', alloc.proportions=c(1,2))
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
sdDP <- function (x, eps, lower.bounds, upper.bounds,
                  which.sensitivity='bounded', mechanism='Laplace', delta=0,
                  type.DP='aDP', alloc.proportions=NULL){
  ########## Input checking

  ##########

  ########## Data Access/privacy layer
  sanitized.variances <- varDP(x,eps,lower.bounds,upper.bounds,which.sensitivity,
                               mechanism,delta,type.DP,alloc.proportions);
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    sanitized.sds <- list()
    sanitized.sds[["Bounded"]] <- sqrt(sanitized.variances$Bounded)
    sanitized.sds[["Unbounded"]] <- sqrt(sanitized.variances$Unbounded)
    class(sanitized.sds) <- "Sanitized Standard Deviation";
  } else sanitized.sds <- sqrt(sanitized.variances)

  return(sanitized.sds);
  ##########
}

#' Differentially Private Covariance
#'
#' This function computes the differentially private covariance of a pair of
#' vectors at user-specified privacy levels of epsilon and delta.
#'
#' @param x1,x2 Numeric vectors.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound1,lower.bound2 Real numbers giving the global or public
#'   lower bounds of x1 and x2, respectively.
#' @param upper.bound1,upper.bound2 Real numbers giving the global or public
#'   upper bounds of x1 and x2, respectively.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @return Sanitized covariance(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' D1 <- sort(rnorm(500, mean=3, sd=2))
#' D2 <- sort(rnorm(500, mean=-1,sd=0.5))
#' lb1 <- -3 # 3 std devs below mean
#' lb2 <- -2.5 # 3 std devs below mean
#' ub1 <- 9 # 3 std devs above mean
#' ub2 <- .5 # 3 std devs above mean
#' covDP(D1, D2, 1, lb1, ub1, lb2, ub2)
#' covDP(D1, D2, .5, lb1, ub1, lb2, ub2, which.sensitivity='unbounded',
#'   mechanism='Gaussian', delta=0.01)
#'
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
covDP <- function (x1, x2, eps, lower.bound1, upper.bound1, lower.bound2,
                  upper.bound2, which.sensitivity='bounded',
                  mechanism='Laplace', delta=0, type.DP='aDP'){
  #### INPUT CHECKING ####
  {if (length(upper.bound1)!=1) stop("Length of upper.bound1 must be 1.");
  if (length(lower.bound1)!=1) stop("Length of lower.bound1 must be 1.");
  if (length(upper.bound2)!=1) stop("Length of upper.bound2 must be 1.");
  if (length(lower.bound2)!=1) stop("Length of lower.bound2 must be 1.");
  x1[x1<lower.bound1] <- lower.bound1;
  x1[x1>upper.bound1] <- upper.bound1;
  x2[x2<lower.bound2] <- lower.bound2;
  x2[x2>upper.bound2] <- upper.bound2;
  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
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
  if (mechanism=='Laplace'){
    sanitized.cov <- LaplaceMechanism(tv,eps,bs,us,which.sensitivity);
  } else if (mechanism=='Gaussian'){
    sanitized.cov <- GaussianMechanism(tv,eps,delta,bs,us,which.sensitivity,
                                       type.DP);
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both') class(sanitized.cov)<-"Sanitized Covariance";
  return(sanitized.cov)
  ##########
}

#' Differentially Private Histogram
#'
#' This function computes a differentially private histogram from a vector at
#' user-specified privacy levels of epsilon and delta. A histogram object is
#' returned with sanitized values for the counts for easy plotting (see
#' examples).
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
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param allow.negative Logical value. If FALSE (default), any negative values
#'   in the sanitized histogram due to the added noise will be set to 0. If
#'   TRUE, the negative values (if any) will be returned.
#' @return Sanitized histogram(s) based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' x <- rnorm(500)
#' hist(x) # Non-private histogram
#' result <- histogramDP(x, 1)
#' plot(result) # Private histogram
#'
#' hist(x, freq=FALSE) # Normalized non-private histogram
#' result <- histogramDP(x, .5, normalize=TRUE, which.sensitivity='unbounded',
#'   mechanism='Gaussian',delta=0.01, allow.negative=TRUE)
#' plot(result) # Normalized private histogram (note negative values allowed)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#' @export
histogramDP <- function(x, eps, breaks="Sturges", normalize=FALSE,
                        which.sensitivity='bounded', mechanism='Laplace',
                        delta=0, type.DP='aDP', allow.negative=FALSE){
  #### INPUT CHECKING ####
  {if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- histogramDataAccess(x, breaks, mechanism)
  tv <- results$True.Values
  bs <- results$Bounded.Sensitivities
  us <- results$Unbounded.Sensitivities
  ##########

  ########## Privacy layer
  counts <- tv$counts
  if (mechanism=='Laplace'){
    sanitized.counts <- LaplaceMechanism(counts,eps,bs,us,which.sensitivity);
  } else if (mechanism=='Gaussian'){
    sanitized.counts <- GaussianMechanism(counts,eps,delta,bs,us,which.sensitivity,
                                          type.DP);
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    sanitized.hist <- list()

    sc <- sanitized.counts$Bounded
    if (!allow.negative) sc[sc<0] <- 0
    if (normalize) {
      sc <- sc/sum(sc)
    } else{
      sc <- round(sc)
    }
    tv$counts <- sc
    sanitized.hist[["Bounded"]] <- tv

    sc <- sanitized.counts$Unbounded;
    if (!allow.negative) sc[sc<0] <- 0;
    if (normalize) {
      sc <- sc/sum(sc);
    } else{
      sc <- round(sc);
    }
    tv$counts <- sc;
    sanitized.hist[["Unbounded"]] <- tv;

    class(sanitized.hist)<-"Sanitized Histogram";
  } else{
    sc <- sanitized.counts
    if (!allow.negative) sc[sc<0] <- 0
    if (normalize) {
      sc <- sc/sum(sc)
    } else{
      sc <- round(sc)
    }
    tv$counts <- sc
    sanitized.hist <- tv
  }

  return(sanitized.hist);
  ##########
}

#' Differentially Private Contingency Table
#'
#' This function computes a differentially private contingency table from given
#' vectors of data at user-specified privacy levels of epsilon and delta.
#'
#' @param ... Vectors of data from which to create the contingency table.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param allow.negative Logical value. If FALSE (default), any negative values
#'   in the sanitized table due to the added noise will be set to 0. If TRUE,
#'   the negative values (if any) will be returned.
#' @return Sanitized contingency table(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' x <- MASS::Cars93$Type;
#' y <- MASS::Cars93$Origin;
#' z <- MASS::Cars93$AirBags;
#' tableDP(x,y,eps=1,which.sensitivity='bounded',mechanism='Laplace',
#'   type.DP='pDP')
#' tableDP(x,y,z,eps=.5,which.sensitivity='unbounded',mechanism='Gaussian',
#'   delta=0.01)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#' @export
tableDP <- function(..., eps=1, which.sensitivity='bounded',
                    mechanism='Laplace', delta=0, type.DP='aDP',
                    allow.negative=FALSE){
  #### INPUT CHECKING ####
  {if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- tableDataAccess(..., mechanism=mechanism);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;

  # Flatten tv, but keep values to recreate later
  named.table <- tv - tv;
  dims <- dim(tv);
  dim(tv) <- NULL;
  ##########

  ########## Privacy layer
  # This means that each param[i] in the mechanism becomes bs/eps rather
  #       than bs[i]/(alloc.proportions[i]*eps)
  # bs <- rep(bs, length(tv))/length(tv)
  # us <- rep(us, length(tv))/length(tv)
  if (mechanism=='Laplace'){
    sanitized.tables <- LaplaceMechanism(tv,eps,bs,us,which.sensitivity);
  } else if (mechanism=='Gaussian'){
    sanitized.tables <- GaussianMechanism(tv,eps,delta,bs,us,which.sensitivity,
                                          type.DP);
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    sanitized.table <- list()

    # Unflatten and round tables
    bounded.table <- sanitized.tables$Bounded;
    dim(bounded.table) <- dims;
    bounded.table <- named.table + bounded.table;
    bounded.table <- round(bounded.table);
    if (!allow.negative) bounded.table[bounded.table<0] <- 0;
    sanitized.table[["Bounded"]] <- bounded.table;

    unbounded.table <- sanitized.tables$Unbounded;
    dim(unbounded.table) <- dims;
    unbounded.table <- named.table + unbounded.table;
    unbounded.table <- round(unbounded.table);
    if (!allow.negative) unbounded.table[unbounded.table<0] <- 0;
    sanitized.table[["Unbounded"]] <- unbounded.table;

    class(sanitized.table)<-"Sanitized Contingency Table";
  } else{
    dim(sanitized.tables) <- dims;
    sanitized.tables <- named.table + sanitized.tables;
    sanitized.tables <- round(sanitized.tables);
    if (!allow.negative) sanitized.tables[sanitized.tables<0] <- 0;
    sanitized.table <- sanitized.tables;
  }

  return(sanitized.table);
  ##########
}

#' Differentially Private Pooled Variance
#'
#' This function computes the differentially private pooled variance from two or
#' more vectors of data at user-specified privacy levels of epsilon and delta.
#'
#' @param ... Two or more vectors from which to compute the pooled variance.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Real number giving the global or public lower bound of the
#'   input data.
#' @param upper.bound Real number giving the global or public upper bound of the
#'   input data.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param approx.n.max Logical indicating whether to approximate n.max (defined
#'   to be the length of the largest input vector) in the computation of the
#'   global sensitivity based on the upper and lower bounds of the data
#'   \insertCite{Liu2019b}{DPpack}. Approximation is best if n.max is very
#'   large.
#' @return Sanitized pooled variance(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' # Build datasets
#' D1 <- rnorm(500, mean=3, sd=2)
#' D2 <- rnorm(200, mean=3, sd=2)
#' D3 <- rnorm(100, mean=3, sd=2)
#' lower.bound <- -3 # 3 standard deviations below mean
#' upper.bound <- 9 # 3 standard deviations above mean
#'
#' # Get private pooled variance without approximating n.max
#' private.pooled.var <- pooledVarDP(D1, D2, D3, eps=1, lower.bound=lower.bound,
#'                                   upper.bound = upper.bound)
#' private.pooled.var
#'
#' # If n.max is sensitive, we can also use
#' private.pooled.var <- pooledVarDP(D1, D2, D3, eps=1, lower.bound=lower.bound,
#'                                   upper.bound = upper.bound,
#'                                   approx.n.max = FALSE)
#' private.pooled.var
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
pooledVarDP <- function(..., eps=1, lower.bound, upper.bound,
                        which.sensitivity='bounded', mechanism='Laplace',
                        delta=0, type.DP='aDP', approx.n.max=FALSE){
  samples <- list(...);
  #### INPUT CHECKING ####
  {J = length(samples);
  if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.");
  if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.");

  for (j in 1:J){
    samples[[j]][samples[[j]]<lower.bound] <- lower.bound;
    samples[[j]][samples[[j]]>upper.bound] <- upper.bound;
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
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
  if (mechanism=='Laplace'){
    while (TRUE){ # Make sure variance is > 0 after noise
      sanitized.vars <- LaplaceMechanism(tv,eps,bs,us,which.sensitivity);
      done <- TRUE
      if (which.sensitivity=='both'){
        if (any(sanitized.vars$Bounded<=0)) done <- FALSE
        if (any(sanitized.vars$Unbounded<=0)) done <- FALSE
      } else if (any(sanitized.vars<=0)) done <- FALSE
      if (done) break;
    }
  } else if (mechanism=='Gaussian'){
    while (TRUE){ # Make sure variance is > 0 after noise
      sanitized.vars <- GaussianMechanism(tv,eps,delta,bs,us,which.sensitivity,
                                          type.DP);
      done <- TRUE
      if (which.sensitivity=='both'){
        if (any(sanitized.vars$Bounded<=0)) done <- FALSE
        if (any(sanitized.vars$Unbounded<=0)) done <- FALSE
      } else if (any(sanitized.vars<=0)) done <- FALSE
      if (done) break;
    }
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both') class(sanitized.vars)<-"Sanitized Pooled Variance";
  return(sanitized.vars);
  ##########
}

#' Differentially Private Pooled Covariance
#'
#' This function computes the differentially private pooled covariance from two
#' or more two-column matrices of data at user-specified privacy levels of
#' epsilon and delta.
#'
#' @param ... Two or more matrices, each with two columns from which to compute
#'   the pooled covariance.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound1,lower.bound2 Real numbers giving the global or public
#'   lower bounds over the first and second columns of all input data,
#'   respectively.
#' @param upper.bound1,upper.bound2 Real numbers giving the global or public
#'   upper bounds over the first and second columns of all input data,
#'   respectively.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param approx.n.max Logical indicating whether to approximate n.max (defined
#'   to be the length of the largest input vector) in the computation of the
#'   global sensitivity based on the upper and lower bounds of the data
#'   \insertCite{Liu2019b}{DPpack}. Approximation is best if n.max is very
#'   large.
#' @return Sanitized pooled covariance(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' # Build datasets
#' D1 <- sort(rnorm(500, mean=3, sd=2))
#' D2 <- sort(rnorm(500, mean=-1, sd=0.5))
#' D3 <- sort(rnorm(200, mean=3, sd=2))
#' D4 <- sort(rnorm(200, mean=-1, sd=0.5))
#' M1 <- matrix(c(D1, D2), ncol=2)
#' M2 <- matrix(c(D3, D4), ncol=2)
#'
#' lb1 <- -3 # 3 std devs below mean
#' lb2 <- -2.5 # 3 std devs below mean
#' ub1 <- 9 # 3 std devs above mean
#' ub2 <- .5 # 3 std devs above mean
#' # Pooled covariance satisfying pure 1-differential privacy
#' private.pooled.cov <- pooledCovDP(M1, M2, eps = 1, lower.bound1 = lb1,
#'                                   lower.bound2 = lb2, upper.bound1 = ub1,
#'                                   upper.bound2 = ub2)
#' private.pooled.cov
#'
#' # Pooled covariance satisfying approximate (0.9, 0.01)-differential privacy
#' # and approximating n.max in the sensitivity calculation
#' private.pooled.cov <- pooledCovDP(M1, M2, eps = 0.9, lower.bound1 = lb1,
#'                                   lower.bound2 = lb2, upper.bound1 = ub1,
#'                                   upper.bound2 = ub2, mechanism = 'Gaussian',
#'                                   delta = 0.01, approx.n.max = TRUE)
#' private.pooled.cov
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{DPtextbook}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
pooledCovDP <- function(..., eps=1, lower.bound1, upper.bound1, lower.bound2,
                        upper.bound2, which.sensitivity='bounded',
                        mechanism='Laplace', delta=0, type.DP='aDP',
                        approx.n.max=FALSE){
  samples <- list(...);
  #### INPUT CHECKING ####
  {J = length(samples);
  if (length(lower.bound1)!=1) stop("Length of lower.bound1 must be 1.");
  if (length(lower.bound2)!=1) stop("Length of lower.bound2 must be 1.");
  if (length(upper.bound1)!=1) stop("Length of upper.bound1 must be 1.");
  if (length(upper.bound2)!=1) stop("Length of upper.bound2 must be 1.");

  for (j in 1:J){
    samples[[j]][samples[[j]][,1]<lower.bound1,1] <- lower.bound1;
    samples[[j]][samples[[j]][,1]>upper.bound1,1] <- upper.bound1;
    samples[[j]][samples[[j]][,2]<lower.bound2,2] <- lower.bound2;
    samples[[j]][samples[[j]][,2]>upper.bound2,2] <- upper.bound2;
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
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
  if (mechanism=='Laplace'){
    sanitized.cov <- LaplaceMechanism(tv,eps,bs,us,which.sensitivity);
  } else if (mechanism=='Gaussian'){
    sanitized.cov <- GaussianMechanism(tv,eps,delta,bs,us,which.sensitivity,
                                       type.DP);
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both') class(sanitized.cov)<-"Sanitized Pooled Covariance";
  return(sanitized.cov);
  ##########
}

#' Differentially Private Quantile
#'
#' This function computes the differentially private quantile of an input vector
#' at user-specified privacy levels of epsilon and delta.
#'
#' @param x Numeric vector of which the quantile will be taken.
#' @param quant Real number between 0 and 1 indicating which quantile to return.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Real number giving the global or public lower bound of x.
#' @param upper.bound Real number giving the global or public upper bound of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'exponential'}.
#'   See \code{\link{ExponentialMechanism}} for a description of the supported
#'   mechanisms.
#' @param uniform.sampling Boolean indicating whether to sample uniformly
#'   between sorted dataset values when returning the private quantile. If TRUE,
#'   it is possible for this function to return any number in the range
#'   [lower.bound, upper.bound]. If FALSE, only a value present in the dataset
#'   or the lower bound can be returned.
#' @return Sanitized quantile(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' D <- rnorm(500)
#' lower.bound <- -3 # 3 standard deviations below mean
#' upper.bound <- 3 # 3 standard deviations above mean
#'
#' quant <- 0.25
#' eps <- 1
#' # Get 25th quantile satisfying pure 1-differential privacy
#' private.quantile <- quantileDP(D, quant, eps, lower.bound, upper.bound)
#' private.quantile
#'
#' # Get 75th quantile requiring released value to be in dataset
#' quant <- 0.75
#' private.quantile <- quantileDP(D, quant, eps, lower.bound, upper.bound,
#'                                uniform.sampling = FALSE)
#' private.quantile
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Smith2011a}{DPpack}
#'
#' @export
quantileDP <- function (x, quant, eps, lower.bound, upper.bound,
                        which.sensitivity='bounded', mechanism='exponential',
                        uniform.sampling=TRUE){
  # NOTE: See
  # https://github.com/IBM/differential-privacy-library/blob/main/diffprivlib/tools/quantiles.py
  #
  # NOTE: The data access and privacy layers are somewhat mixed.

  #### INPUT CHECKING ####
  {if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.");
  if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.");
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
    sanitized.indices <- ExponentialMechanism(utility, eps, bs, us,
                                              which.sensitivity,
                                              measure=diff(sorted));

    if (which.sensitivity=='both'){
      sanitized.quantile <- list();

      bounded.idx <- sanitized.indices$Bounded
      if (uniform.sampling){
        bounded.sanitized <- runif(1)*(sorted[bounded.idx+1]-sorted[bounded.idx]) +
          sorted[bounded.idx]
      } else bounded.sanitized <- sorted[bounded.idx]
      sanitized.quantile[["Bounded"]] <- bounded.sanitized

      unbounded.idx <- sanitized.indices$Unbounded
      if (uniform.sampling){
        unbounded.sanitized <- runif(1)*(sorted[unbounded.idx+1]-sorted[unbounded.idx]) +
          sorted[unbounded.idx]
      } else unbounded.sanitized <- sorted[unbounded.idx]
      sanitized.quantile[["Unbounded"]] <- unbounded.sanitized;

      class(sanitized.quantile)<-"Sanitized Quantile";
    } else{
      if (uniform.sampling){
        sanitized.quantile <- runif(1)*(sorted[sanitized.indices+1]-sorted[sanitized.indices]) +
          sorted[sanitized.indices]
      } else sanitized.quantile <- sorted[sanitized.indices]
    }
  }
  ##########

  ########## Postprocessing layer
  return(sanitized.quantile);
  ##########
}

#' Differentially Private Median
#'
#' This function computes the differentially private median of an input vector
#' at user-specified privacy levels of epsilon and delta.
#'
#' @param x Numeric vector of which the median will be taken.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Real number giving the global or public lower bound of x.
#' @param upper.bound Real number giving the global or public upper bound of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively and in composition. Care must be taken not to violate
#'   differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'exponential'}.
#'   See \code{\link{ExponentialMechanism}} for a description of the supported
#'   mechanisms.
#' @param uniform.sampling Boolean indicating whether to sample uniformly
#'   between sorted dataset values when returning the private median. If TRUE,
#'   it is possible for this function to return any number in the range
#'   [lower.bound, upper.bound]. If FALSE, only a value present in the dataset
#'   or the lower.bound can be returned.
#' @return Sanitized median(s) based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' D <- rnorm(500)
#' lower.bound <- -3 # 3 standard deviations below mean
#' upper.bound <- 3 # 3 standard deviations above mean
#'
#' eps <- 1
#' # Get median satisfying pure 1-differential privacy
#' private.median <- medianDP(D, eps, lower.bound, upper.bound)
#' private.median
#'
#' # Require released value to be in dataset
#' private.median <- medianDP(c(1,0,3,3,2), eps, 0, 4, uniform.sampling = FALSE)
#' private.median
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Smith2011a}{DPpack}
#'
#' @export
medianDP <- function (x, eps, lower.bound, upper.bound,
                      which.sensitivity='bounded', mechanism='exponential',
                      uniform.sampling=TRUE){
  sanitized.median <- quantileDP(x,.5,eps,lower.bound,upper.bound,
                                 which.sensitivity,mechanism,uniform.sampling)
  if (which.sensitivity=='both') class(sanitized.median)<-"Sanitized Median"
  return(sanitized.median);
  ##########
}





