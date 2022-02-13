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
#' meanDP(c(1,4,-2,8,-6),1,lower.bounds=-10,upper.bounds=10,type.DP='pDP')
#' meanDP(c(1,4,-2,8,-6),.5,lower.bounds=-10,upper.bounds=10,
#'   which.sensitivity='unbounded',mechanism='Gaussian',
#'   delta=0.5)
#' meanDP(matrix(c(1,4,-2,8,-6,0),ncol=2),1,lower.bounds=c(-10,-10),
#'   upper.bounds=c(10,10),which.sensitivity='bounded',type.DP='pDP',
#'   alloc.proportions=c(1,2))
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
#' varDP(c(1,4,-2,8,-6),1,lower.bounds=-10,upper.bounds=10,type.DP='pDP')
#' varDP(c(1,4,-2,8,-6),.5,lower.bounds=-10,upper.bounds=10,
#'   which.sensitivity='unbounded',mechanism='Gaussian',delta=0.5)
#' varDP(matrix(c(1,4,-2,8,-6,0),ncol=2),1,lower.bounds=c(-10,-10),
#'   upper.bounds=c(10,10),which.sensitivity='bounded',type.DP='pDP',
#'   alloc.proportions=c(1,2))
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
#' sdDP(c(1,4,-2,8,-6),1,lower.bounds=-10,upper.bounds=10,type.DP='pDP')
#' sdDP(c(1,4,-2,8,-6),.5,lower.bounds=-10,upper.bounds=10,
#'   which.sensitivity='unbounded',mechanism='Gaussian',delta=0.5)
#' sdDP(matrix(c(1,4,-2,8,-6,0),ncol=2),1,lower.bounds=c(-10,-10),
#'   upper.bounds=c(10,10),which.sensitivity='bounded',type.DP='pDP',
#'   alloc.proportions=c(1,2))
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
#' covDP(c(1,4,-2,8,-6),c(1,3,2,2,4),1,lower.bound1=-10,upper.bound1=10,
#'   lower.bound2=0,upper.bound2=5,which.sensitivity='bounded',
#'   mechanism='Laplace',type.DP='pDP')
#' covDP(c(1,4,-2,8,-6),c(1,3,2,2,4),.5,lower.bound1=-10,upper.bound110,
#'   lower.bound2=0,upper.bound2=5,which.sensitivity='unbounded',
#'   mechanism='Gaussian',delta=0.5)
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
#' result <- histogramDP(c(1,1,-2,8,-6),1,which.sensitivity='bounded',
#'   mechanism='Laplace',type.DP='pDP')
#' plot(result)
#' result <- histogramDP(c(1,1,-2,8,-6),.5,normalize=TRUE,
#'   which.sensitivity='unbounded',mechanism='Gaussian',delta=0.5,
#'   allow.negative=FALSE)
#' plot(result)
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
  results <- histogramDataAccess(x, breaks);
  tv <- results$True.Values;
  bs <- results$Bounded.Sensitivities;
  us <- results$Unbounded.Sensitivities;
  ##########

  ########## Privacy layer
  counts <- tv$counts;
  # This means that each param[i] in the mechanism becomes bs/eps rather
  #       than bs[i]/(alloc.proportions[i]*eps)
  # bs <- rep(bs, length(counts))/length(counts)
  # us <- rep(us, length(counts))/length(counts)
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
#'   delta=0.5,allow.negative=FALSE)
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
tableDP <- function(..., eps, which.sensitivity='bounded', mechanism='Laplace',
                    delta=0, type.DP='aDP', allow.negative=FALSE){
  #### INPUT CHECKING ####
  {if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace';
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.");
  }
  }
  ##########

  ########## Data access layer
  results <- tableDataAccess(...);
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
#' pooledVarDP(c(1,4,-2,8,-6),c(1,2),c(-5,-7),eps=1,lower.bound=-10,
#'   upper.bound=10,which.sensitivity='bounded',mechanism='Laplace',
#'   type.DP='pDP')
#' pooledVarDP(c(1,4,-2,8,-6),c(1,2),c(-5,-7),eps=.5,
#'   lower.bound=-10,upper.bound=10,which.sensitivity='unbounded',
#'   mechanism='Gaussian',delta=0.5,approx.n.max=TRUE)
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
#' x1 <- matrix(c(1,4,-2,8,-6,-3),ncol=2)
#' x2 <- matrix(c(1,2,-5,7),ncol=2)
#' pooledVarDP(x1,x2,eps=1,lower.bound1=-10,upper.bound1=10,lower.bound2=-10,
#'   upper.bound2=10,which.sensitivity='bounded',mechanism='Laplace',
#'   type.DP='pDP')
#' pooledVarDP(x1,x2,eps=.5,lower.bound1=-10,upper.bound1=10,lower.bound2=-10,
#'   upper.bound2=10,which.sensitivity='unbounded',mechanism='Gaussian',
#'   delta=0.5,approx.n.max=TRUE)
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
#' @return Sanitized quantile(s) based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' quantileDP(c(1,1,-2,8,-6),.25,1,lower.bound=-10,upper.bound=10,
#'   which.sensitivity='bounded',mechanism='exponential')
#' quantileDP(c(1,1,-2,8,-6),.75,1,lower.bound=-10,upper.bound=10,
#'   which.sensitivity='unbounded',mechanism='exponential')
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Smith2011a}{DPpack}
#'
#' @export
quantileDP <- function (x, quant, eps, lower.bound, upper.bound,
                        which.sensitivity='bounded', mechanism='exponential'){
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

      bounded.idx <- sanitized.indices$Bounded;
      bounded.sanitized <- runif(1)*(sorted[bounded.idx+1]-sorted[bounded.idx]) +
        sorted[bounded.idx];
      sanitized.quantile[["Bounded"]] <- bounded.sanitized;

      unbounded.idx <- sanitized.indices$Unbounded;
      unbounded.sanitized <- runif(1)*(sorted[unbounded.idx+1]-sorted[unbounded.idx]) +
        sorted[unbounded.idx];
      sanitized.quantile[["Unbounded"]] <- unbounded.sanitized;

      class(sanitized.quantile)<-"Sanitized Quantile";
    } else{
      sanitized.quantile <- runif(1)*(sorted[sanitized.indices+1]-sorted[sanitized.indices]) +
        sorted[sanitized.indices]
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
#' @return Sanitized median(s) based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' medianDP(c(1,1,-2,8,-6),1,lower.bound=-10,upper.bound=10,
#'   which.sensitivity='bounded',mechanism='exponential')
#' medianDP(c(1,1,-2,8,-6),1,lower.bound=-10,upper.bound=10,
#'   which.sensitivity='unbounded',mechanism='exponential')
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Smith2011a}{DPpack}
#'
#' @export
medianDP <- function (x, eps, lower.bound, upper.bound,
                      which.sensitivity='bounded', mechanism='exponential'){
  sanitized.median <- quantileDP(x,.5,eps,lower.bound,upper.bound,
                                 which.sensitivity,mechanism)
  if (which.sensitivity=='both') class(sanitized.median)<-"Sanitized Median"
  return(sanitized.median);
  ##########
}





