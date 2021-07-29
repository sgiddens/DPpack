#' Compute differentially private mean
#'
#' TODO: Add description.
#'
#' @param x The dataset over which the mean will be taken. Means taken over columns.
#' @param eps The epsilon privacy budget.
#' @param which.sensitivity Can be one of {'bounded', 'unbounded', 'both'}. If
#'      'bounded' (default), returns result plus noise based on bounded
#'      sensitivities. If 'unbounded', returns result plus noise based on
#'      unbounded sensitivities.
#' @param lower.bounds The list of lower bounds (minimums) of the full dataset
#'      from which x was taken. If NULL, it is computed to be the min value of
#'      each column of x.
#' @param upper.bounds The list of upper bound (maximum) of the full dataset
#'      from which x was taken. If NULL, it is computed to be the min value of
#'      each column of x.
#' @param mechanism A string indicating which mechanism to use for differential
#'      privacy. Currently the following mechanisms are supported
#'      {'laplace', 'gaussian'}.
#' @param delt The delta privacy parameter (necessary if using Gaussian mechanism).
#' @param type.DP The type of differential privacy desired for Gaussian mechanism.
#'      Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
#' @param alloc.proportions This gives the proportional allocation of epsilon (and
#'       delta) to the statistics. By default, it distributes eps and delt
#'       evenly among the calculations. Input does not need to be normalized.
#' @return List of bounded and/or unbounded sanitized means.
#' @examples
#' meanDP(c(1,4,-2,8,-6),1,'bounded',-10,10,'laplace')
#' meanDP(c(1,4,-2,8,-6),1,'unbounded',-10,10,'gaussian',0.5,'aDP')
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

#' @export
varDP <- function (x, eps, which.sensitivity='bounded', lower.bounds=NULL,
                   upper.bounds=NULL, mechanism='laplace', delt=NULL,
                   type.DP='pDP', alloc.proportions=NULL){
  # Computes differentially private variance.
  #
  # x: the dataset over which the variance will be taken.
  #     It could be a subset of the full dataset. Variances taken over columns.
  # eps: the epsilon privacy budget.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bounds: the list of lower bound (minimum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of each
  #     column of x.
  # upper.bounds: the list of upper bound (maximum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of each
  #     column of x.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  # alloc.proportions: This gives the proportional allocation of epsilon (and
  #       delta) to the statistics. By default, it distributes eps and delt
  #       evenly among the calculations. Input does not need to be normalized.
  #
  # Returns: List of bounded and/or unbounded sanitized variances.

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

#' @export
sdDP <- function (x, eps, which.sensitivity='bounded', lower.bounds=NULL,
                  upper.bounds=NULL, mechanism='laplace', delt=NULL,
                  type.DP='pDP', alloc.proportions=NULL){
  # Computes differentially private sample standard deviation.
  #
  # x: the dataset over which the standard deviation will be taken.
  #     It could be a subset of the full dataset. Standard deviations taken over
  #     columns.
  # eps: the epsilon privacy budget
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bounds: the list of lower bound (minimum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of each
  #     column of x.
  # upper.bounds: the list of upper bound (maximum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of each
  #     column of x.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  # alloc.proportions: This gives the proportional allocation of epsilon (and
  #       delta) to the statistics. By default, it distributes eps and delt
  #       evenly among the calculations. Input does not need to be normalized.
  #
  # Returns: List of bounded and/or unbounded sanitized standard deviations.
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

#' @export
covDP <- function (x1, x2, eps, which.sensitivity='bounded',
                  lower.bound1=NULL, upper.bound1=NULL,
                  lower.bound2=NULL, upper.bound2=NULL,
                  mechanism='laplace', delt=NULL, type.DP='pDP'){
  # Computes differentially private sample covariance between x1 and x2.
  #
  # x1, x2: the datasets between which the covariance will be taken.
  #     Each could be a subset of the full dataset. Must be single dimensional.
  # eps: the epsilon privacy budget
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bound1, lower_bound2: the lower bound (minimum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of x1, x2
  # upper.bound1, upper_bound2: the upper bound (maximum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the max value of x1, x2
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  #
  # Returns: List of bounded and/or unbounded sanitized covariance.

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

#' @export
histogramDP <- function(x, eps, breaks="Sturges", normalize=FALSE,
                        which.sensitivity='bounded',
                        lower.bound=NULL, upper.bound=NULL,
                        mechanism='laplace', delt=NULL, type.DP='pDP',
                        allow.negative=FALSE){
  # Computes a differentially private histogram of x.
  #
  # x: the dataset over which to take this histogram.
  # eps: the epsilon privacy budget.
  # breaks: same argument as with hist function.
  # normalize: Whether or not to normalize the histogram.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bound: the lower bound (minimum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of x.
  # upper.bound: the upper bound (maximum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the max value of x.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  # allow.negative: Whether to allow negative values in the sanitized counts.
  #     These values may occur as a result of noise and may be desired for
  #     some applications to avoid biasing results. If set to FALSE,
  #     forces negative counts to be 0.
  #
  # Returns: List of bounded and/or unbounded sanitized histograms.

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

#' @export
tableDP <- function(x, y, eps, which.sensitivity='bounded', mechanism='laplace',
                    delt=NULL, type.DP='pDP', allow.negative=FALSE){
  # Computes a differentially private contingency table of x and y.
  #
  # x: the first set of data over which to create the contingency table.
  # y: the second set of data over which to create the contingency table.
  # eps: the epsilon privacy budget.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  # allow.negative: Whether to allow negative values in the sanitized counts.
  #     These values may occur as a result of noise and may be desired for
  #     some applications to avoid biasing results. If set to FALSE,
  #     forces negative counts to be 0.
  #
  # Returns: List of bounded and/or unbounded sanitized contingency tables.

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

#' @export
pooledVarDP <- function(..., eps=1, which.sensitivity='bounded',
                        lower.bound=NULL, upper.bound=NULL,
                        mechanism='laplace', delt=NULL, type.DP='pDP',
                        approx.n.max=FALSE){
  # Computes a differentially private pooled variance of a collection of samples.
  #
  # ...: the collection of samples from which to compute the pooled variance.
  # eps: the epsilon privacy budget.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bound: the lower bound (minimum) of the collections.
  #     If NULL, it is computed to be the min value across all of them.
  # upper.bound: the upper bound (maximum) of the collections.
  #     If NULL, it is computed to be the max value across all of them.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  # approx.n.max: Boolean indicating whether to approximate n.max. Approximation
  #     is best if n.max is very large.
  #
  # Returns: List of bounded and/or unbounded sanitized pooled variances.
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

#' @export
pooledCovDP <- function(..., eps=1, which.sensitivity='bounded',
                        lower.bound1=NULL, upper.bound1=NULL,
                        lower.bound2=NULL, upper.bound2=NULL,
                        mechanism='laplace', delt=NULL, type.DP='pDP',
                        approx.n.max=FALSE){
  # Computes a differentially private pooled covariance of a collection of samples.
  #
  # ...: the collection of samples from which to compute the pooled covariance.
  #       It is expected that each passed in input is of shape (nj,2), where
  #       nj>=2.
  # eps: the epsilon privacy budget.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bound1: the lower bounds (minimums) of each of the x1 collections.
  #     If NULL, it is computed to be the min values of each of them.
  # upper.bound1: the upper bounds (maximums) of each of the x1 collections.
  #     If NULL, it is computed to be the max values of each of them.
  # lower.bound2: the lower bounds (minimums) of each of the x2 collections.
  #     If NULL, it is computed to be the min values of each of them.
  # upper.bound2: the upper bounds (maximums) of each of the x2 collections.
  #     If NULL, it is computed to be the max values of each of them.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'laplace', 'gaussian'}.
  # delt: The delta privacy parameter (for pDP and aDP with Gaussian mechanism).
  # type.DP: The type of differential privacy desired for Gaussian mechanism.
  #     Can be either probabilistic DP ('pDP') or approximate DP ('aDP').
  # approx.n.max: Boolean indicating whether to approximate n.max. Approximation
  #     is best if n.max is very large.
  #
  # Returns: List of bounded and/or unbounded sanitized pooled covariances.
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

#' @export
quantileDP <- function (x, quant, eps, which.sensitivity='bounded',
                        lower.bound=NULL, upper.bound=NULL,
                        mechanism='exponential', delt=NULL){
  # Computes differentially private quantile.
  # NOTE: See
  # https://github.com/IBM/differential-privacy-library/blob/main/diffprivlib/tools/quantiles.py
  #
  # NOTE: The data access and privacy layers are somewhat mixed.
  #
  # x: the dataset of which the quantile will be taken.
  #     It could be a subset of the full dataset.
  # quant: Real number between 0 and 1 indicating which quantile to return.
  # eps: the epsilon privacy budget.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bound: the lower bound (minimum) of the full dataset from
  #     which x was taken.
  # upper.bound: the upper bound (maximum) of the full dataset from
  #     which x was taken.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'exponential'}.
  # delt: The delta privacy parameter (not currently used).
  #
  # Returns: List of bounded and/or unbounded sanitized quantiles.

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

#' @export
medianDP <- function (x, eps, which.sensitivity='bounded',
                      lower.bounds=NULL, upper.bounds=NULL,
                      mechanism='exponential', delt=NULL){
  # Computes differentially private median
  #
  # x: the dataset of which the median will be taken.
  #     It could be a subset of the full dataset.
  # eps: the epsilon privacy budget.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus noise based on bounded
  #       sensitivities. If 'unbounded', returns result plus noise based on
  #       unbounded sensitivities.
  # lower.bounds: the list of lower bound (minimum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of each
  #     column of x.
  # upper.bounds: the list of upper bound (maximum) of the full dataset from
  #     which x was taken. If NULL, it is computed to be the min value of each
  #     column of x.
  # mechanism: a string indicating which mechanism to use for differential
  #     privacy. Currently the following mechanisms are supported
  #     {'exponential'}.
  # delt: The delta privacy parameter (not currently used).
  #
  # Returns: List of bounded and unbounded sanitized medians.
  sanitized.median <- quantileDP(x,.5,eps,which.sensitivity,lower.bounds,
                                 upper.bounds,mechanism,delt);
  class(sanitized.median)<-"Sanitized Median";
  return(sanitized.median);
  ##########
}

#' @export
# logisticRegDP <- function(X, y, ){
# testing
# }






