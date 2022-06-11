#' Differentially Private Mean
#'
#' This function computes the differentially private mean of a given dataset at
#' user-specified privacy levels of epsilon and delta.
#'
#' @param x Dataset whose mean is desired.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Scalar representing the global or public lower bound on
#'   values of x.
#' @param upper.bound Scalar representing the global or public upper bound on
#'   values of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @return Sanitized mean based on the bounded and/or unbounded definitions of
#'   differential privacy.
#' @examples
#' D <- stats::rnorm(500, mean=3, sd=2)
#' lb <- -3 # 3 std devs below mean
#' ub <- 9 # 3 std devs above mean
#' meanDP(D,  1, lb, ub)
#' meanDP(D, .5, lb, ub, which.sensitivity='unbounded', mechanism='Gaussian',
#'   delta=0.01)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#' @export
meanDP <- function (x, eps, lower.bound, upper.bound,
                    which.sensitivity='bounded', mechanism='Laplace', delta=0,
                    type.DP='aDP'){
  #### INPUT CHECKING ####
  {
  if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.")
  if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.")
  x[x<lower.bound] <- lower.bound
  x[x>upper.bound] <- upper.bound

  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }
  if (which.sensitivity=='both'){
    warning("Sensitivity based on bounded and unbounded differential privacy is identical for this statistic. Only one value will be returned.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- meanDataAccess(x, lower.bound=lower.bound, upper.bound=upper.bound)
  tv <- results$True.Values
  sens <- results$Sensitivity
  ##########

  ########## Privacy layer
  if (mechanism=='Laplace'){
    sanitized.mean <- LaplaceMechanism(tv, eps, sens)
  } else if (mechanism=='Gaussian'){
    sanitized.mean <- GaussianMechanism(tv, eps, delta, sens, type.DP)
  }
  ##########

  ########## Postprocessing layer
  return(sanitized.mean)
  ##########
}

#' Differentially Private Variance
#'
#' This function computes the differentially private variance of a given dataset
#' at user-specified privacy levels of epsilon and delta.
#'
#' @param x Numeric vector whose variance is desired.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Scalar representing the global or public lower bound on
#'   values of x.
#' @param upper.bound Scalar representing the global or public upper bound on
#'   values of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @return Sanitized variance based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' D <- stats::rnorm(500, mean=3, sd=2)
#' lb <- -3 # 3 std devs below mean
#' ub <- 9 # 3 std devs above mean
#' varDP(D, 1, lb, ub)
#' varDP(D,.5, lb, ub, which.sensitivity='unbounded', mechanism='Gaussian',
#'   delta=0.01)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
varDP <- function (x, eps, lower.bound, upper.bound,
                   which.sensitivity='bounded', mechanism='Laplace', delta=0,
                   type.DP='aDP'){
  #### INPUT CHECKING ####
  {
  if (!is.vector(x)) stop("x must be a numeric vector.")

  if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.")
  if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.")
  x[x<lower.bound] <- lower.bound
  x[x>upper.bound] <- upper.bound

  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }
  if (which.sensitivity=='both'){
    warning("Sensitivity based on bounded and unbounded differential privacy is identical for this statistic. Only one value will be returned.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- varDataAccess(x, lower.bound=lower.bound, upper.bound=upper.bound)
  tv <- results$True.Values
  sens <- results$Sensitivity
  ##########

  ########## Privacy layer
  if (mechanism=='Laplace'){
    sanitized.var <- -1
    while (sanitized.var<=0){ # Make sure variance is > 0 after noise
      sanitized.var <- LaplaceMechanism(tv, eps, sens)
    }
  }  else if (mechanism=='Gaussian'){
    sanitized.var <- -1
    while (sanitized.var<=0){
      sanitized.var <- GaussianMechanism(tv, eps, delta, sens, type.DP)
    }
  }
  ##########

  ########## Postprocessing layer
  return(sanitized.var)
  ##########
}

#' Differentially Private Standard Deviation
#'
#' This function computes the differentially private standard deviation of a
#' given dataset at user-specified privacy levels of epsilon and delta.
#'
#' @param x Numeric vector whose variance is desired.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Scalar representing the global or public lower bound on
#'   values of x.
#' @param upper.bound Scalar representing the global or public upper bound on
#'   values of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @return Sanitized standard deviation based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' D <- stats::rnorm(500, mean=3, sd=2)
#' lb <- -3 # 3 std devs below mean
#' ub <- 9 # 3 std devs above mean
#' sdDP(D, 1, lb, ub)
#' sdDP(D,.5, lb, ub, which.sensitivity='unbounded', mechanism='Gaussian',
#'   delta=0.01)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
sdDP <- function (x, eps, lower.bound, upper.bound,
                  which.sensitivity='bounded', mechanism='Laplace', delta=0,
                  type.DP='aDP'){
  ########## Input checking

  ##########

  ########## Data Access/privacy layer
  sanitized.variance <- varDP(x, eps, lower.bound, upper.bound,
                              which.sensitivity ,mechanism, delta, type.DP)
  ##########

  ########## Postprocessing layer
  sanitized.sd <- sqrt(sanitized.variance)

  return(sanitized.sd)
  ##########
}

#' Differentially Private Covariance
#'
#' This function computes the differentially private covariance of a pair of
#' vectors at user-specified privacy levels of epsilon and delta.
#'
#' @param x1,x2 Numeric vectors whose covariance is desired.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound1,lower.bound2 Real numbers giving the global or public
#'   lower bounds of x1 and x2, respectively.
#' @param upper.bound1,upper.bound2 Real numbers giving the global or public
#'   upper bounds of x1 and x2, respectively.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @return Sanitized covariance based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' D1 <- sort(stats::rnorm(500, mean=3, sd=2))
#' D2 <- sort(stats::rnorm(500, mean=-1,sd=0.5))
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
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
covDP <- function (x1, x2, eps, lower.bound1, upper.bound1, lower.bound2,
                  upper.bound2, which.sensitivity='bounded',
                  mechanism='Laplace', delta=0, type.DP='aDP'){
  #### INPUT CHECKING ####
  {
  if (!is.vector(x1)) stop("x1 must be a numeric vector.")
  if (!is.vector(x2)) stop("x2 must be a numeric vector.")

  if (length(upper.bound1)!=1) stop("Length of upper.bound1 must be 1.")
  if (length(lower.bound1)!=1) stop("Length of lower.bound1 must be 1.")
  if (length(upper.bound2)!=1) stop("Length of upper.bound2 must be 1.")
  if (length(lower.bound2)!=1) stop("Length of lower.bound2 must be 1.")
  x1[x1<lower.bound1] <- lower.bound1
  x1[x1>upper.bound1] <- upper.bound1
  x2[x2<lower.bound2] <- lower.bound2
  x2[x2>upper.bound2] <- upper.bound2

  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }
  if (which.sensitivity=='both'){
    warning("Sensitivity based on bounded and unbounded differential privacy is identical for this statistic. Only one value will be returned.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- covDataAccess(x1, x2, lower.bound1, upper.bound1,
                           lower.bound2, upper.bound2)
  tv <- results$True.Values
  sens <- results$Sensitivity
  ##########

  ########## Privacy layer
  if (mechanism=='Laplace'){
    sanitized.cov <- LaplaceMechanism(tv, eps, sens)
  } else if (mechanism=='Gaussian'){
    sanitized.cov <- GaussianMechanism(tv, eps, delta, sens, type.DP)
  }
  ##########

  ########## Postprocessing layer
  return(sanitized.cov)
  ##########
}

#' Differentially Private Histogram
#'
#' This function computes a differentially private histogram from a vector at
#' user-specified privacy levels of epsilon and delta. A histogram object is
#' returned with sanitized values for the counts for easy plotting.
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
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param allow.negative Logical value. If FALSE (default), any negative values
#'   in the sanitized histogram due to the added noise will be set to 0. If
#'   TRUE, the negative values (if any) will be returned.
#' @return Sanitized histogram based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' x <- stats::rnorm(500)
#' graphics::hist(x) # Non-private histogram
#' result <- histogramDP(x, 1)
#' plot(result) # Private histogram
#'
#' graphics::hist(x, freq=FALSE) # Normalized non-private histogram
#' result <- histogramDP(x, .5, normalize=TRUE, which.sensitivity='unbounded',
#'   mechanism='Gaussian',delta=0.01, allow.negative=TRUE)
#' plot(result) # Normalized private histogram (note negative values allowed)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#' @export
histogramDP <- function(x, eps, breaks="Sturges", normalize=FALSE,
                        which.sensitivity='bounded', mechanism='Laplace',
                        delta=0, type.DP='aDP', allow.negative=FALSE){
  #### INPUT CHECKING ####
  {
  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
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
    if (which.sensitivity=='bounded'){
      sanitized.counts <- LaplaceMechanism(counts, eps, bs)
    } else if (which.sensitivity=='unbounded'){
      sanitized.counts <- LaplaceMechanism(counts, eps, us)
    } else if (which.sensitivity=='both'){
      sanitized.counts.bounded <- LaplaceMechanism(counts, eps, bs)
      sanitized.counts.unbounded <- LaplaceMechanism(counts, eps, us)
    }
  } else if (mechanism=='Gaussian'){
    if (which.sensitivity=='bounded'){
      sanitized.counts <- GaussianMechanism(counts, eps, delta, bs, type.DP)
    } else if (which.sensitivity=='unbounded'){
      sanitized.counts <- GaussianMechanism(counts, eps, delta, us, type.DP)
    } else if (which.sensitivity=='both'){
      sanitized.counts.bounded <- GaussianMechanism(counts, eps, delta, bs,
                                                    type.DP)
      sanitized.counts.unbounded <- GaussianMechanism(counts, eps, delta, us,
                                                      type.DP)
    }
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    sanitized.hist <- list()

    sc <- sanitized.counts.bounded
    if (!allow.negative) sc[sc<0] <- 0
    if (normalize) {
      sc <- sc/sum(sc)
    } else{
      sc <- round(sc)
    }
    tv$counts <- sc
    sanitized.hist[["Bounded"]] <- tv

    sc <- sanitized.counts.unbounded
    if (!allow.negative) sc[sc<0] <- 0
    if (normalize) {
      sc <- sc/sum(sc)
    } else{
      sc <- round(sc)
    }
    tv$counts <- sc
    sanitized.hist[["Unbounded"]] <- tv

    class(sanitized.hist)<-"Sanitized Histogram"
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

  return(sanitized.hist)
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
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param allow.negative Logical value. If FALSE (default), any negative values
#'   in the sanitized table due to the added noise will be set to 0. If TRUE,
#'   the negative values (if any) will be returned.
#' @return Sanitized contingency table based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' x <- MASS::Cars93$Type
#' y <- MASS::Cars93$Origin
#' z <- MASS::Cars93$AirBags
#' tableDP(x,y,eps=1,which.sensitivity='bounded',mechanism='Laplace',
#'   type.DP='pDP')
#' tableDP(x,y,z,eps=.5,which.sensitivity='unbounded',mechanism='Gaussian',
#'   delta=0.01)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#' @export
tableDP <- function(..., eps=1, which.sensitivity='bounded',
                    mechanism='Laplace', delta=0, type.DP='aDP',
                    allow.negative=FALSE){
  #### INPUT CHECKING ####
  {
  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- tableDataAccess(..., mechanism=mechanism)
  tv <- results$True.Values
  bs <- results$Bounded.Sensitivities
  us <- results$Unbounded.Sensitivities

  # Flatten tv, but keep values to recreate later
  named.table <- tv - tv
  dims <- dim(tv)
  dim(tv) <- NULL
  ##########

  ########## Privacy layer
  if (mechanism=='Laplace'){
    if (which.sensitivity=='bounded'){
      sanitized.table <- LaplaceMechanism(tv, eps, bs)
    } else if (which.sensitivity=='unbounded'){
      sanitized.table <- LaplaceMechanism(tv, eps, us)
    } else if (which.sensitivity=='both'){
      sanitized.table.bounded <- LaplaceMechanism(tv, eps, bs)
      sanitized.table.unbounded <- LaplaceMechanism(tv, eps, us)
    }
  } else if (mechanism=='Gaussian'){
    if (which.sensitivity=='bounded'){
      sanitized.table <- GaussianMechanism(tv, eps, delta, bs, type.DP)
    } else if (which.sensitivity=='unbounded'){
      sanitized.table <- GaussianMechanism(tv, eps, delta, us, type.DP)
    } else if (which.sensitivity=='both'){
      sanitized.table.bounded <- GaussianMechanism(tv, eps, delta, bs,
                                                    type.DP)
      sanitized.table.unbounded <- GaussianMechanism(tv, eps, delta, us,
                                                      type.DP)
    }
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    out <- list()

    # Unflatten and round tables
    bounded.table <- sanitized.table.bounded
    dim(bounded.table) <- dims
    bounded.table <- named.table + bounded.table
    bounded.table <- round(bounded.table)
    if (!allow.negative) bounded.table[bounded.table<0] <- 0
    out[["Bounded"]] <- bounded.table

    unbounded.table <- sanitized.table.unbounded
    dim(unbounded.table) <- dims
    unbounded.table <- named.table + unbounded.table
    unbounded.table <- round(unbounded.table)
    if (!allow.negative) unbounded.table[unbounded.table<0] <- 0
    out[["Unbounded"]] <- unbounded.table

    class(out)<-"Sanitized Contingency Table"
  } else{
    dim(sanitized.table) <- dims
    sanitized.table <- named.table + sanitized.table
    sanitized.table <- round(sanitized.table)
    if (!allow.negative) sanitized.table[sanitized.table<0] <- 0
    out <- sanitized.table
  }

  return(out)
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
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param approx.n.max Logical indicating whether to approximate n.max (defined
#'   to be the length of the largest input vector) in the computation of the
#'   global sensitivity based on the upper and lower bounds of the data
#'   \insertCite{Liu2019b}{DPpack}. Approximation is best if n.max is very
#'   large.
#' @return Sanitized pooled variance based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' # Build datasets
#' D1 <- stats::rnorm(500, mean=3, sd=2)
#' D2 <- stats::rnorm(200, mean=3, sd=2)
#' D3 <- stats::rnorm(100, mean=3, sd=2)
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
#'                                   approx.n.max = TRUE)
#' private.pooled.var
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
pooledVarDP <- function(..., eps=1, lower.bound, upper.bound,
                        which.sensitivity='bounded', mechanism='Laplace',
                        delta=0, type.DP='aDP', approx.n.max=FALSE){
  samples <- list(...)
  #### INPUT CHECKING ####
  {
  J = length(samples)
  if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.")
  if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.")

  for (j in 1:J){
    samples[[j]][samples[[j]]<lower.bound] <- lower.bound
    samples[[j]][samples[[j]]>upper.bound] <- upper.bound
  }

  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- pooledVarDataAccess(samples, lower.bound=lower.bound,
                                 upper.bound=upper.bound,
                                 approx.n.max=approx.n.max)
  tv <- results$True.Values
  bs <- results$Bounded.Sensitivities
  us <- results$Unbounded.Sensitivities

  # Bounded/unbounded sensitivities are almost always equal for this statistic
  if (which.sensitivity=='both' & bs==us) {
    warning("Sensitivity based on bounded and unbounded differential privacy is identical for this statistic. Only one value will be returned.")
    which.sensitivity <- 'bounded'
  }
  ##########

  ########## Privacy layer
  if (mechanism=='Laplace'){
    if (which.sensitivity=='bounded'){
      sanitized.var <- -1
      while (sanitized.var<=0){ # Make sure variance is > 0 after noise
        sanitized.var <- LaplaceMechanism(tv, eps, bs)
      }
    } else if (which.sensitivity=='unbounded'){
      sanitized.var <- -1
      while (sanitized.var<=0){ # Make sure variance is > 0 after noise
        sanitized.var <- LaplaceMechanism(tv, eps, us)
      }
    } else if (which.sensitivity=='both'){
      sanitized.var.bounded <- -1
      while (sanitized.var.bounded<=0){ # Make sure variance is > 0 after noise
        sanitized.var.bounded <- LaplaceMechanism(tv, eps, bs)
      }
      sanitized.var.unbounded <- -1
      while (sanitized.var.unbounded<=0){ # Make sure variance is > 0 after noise
        sanitized.var.unbounded <- LaplaceMechanism(tv, eps, us)
      }
    }
  } else if (mechanism=='Gaussian'){
    if (which.sensitivity=='bounded'){
      sanitized.var <- -1
      while (sanitized.var<=0){ # Make sure variance is > 0 after noise
        sanitized.var <- GaussianMechanism(tv, eps, delta, bs, type.DP)
      }
    } else if (which.sensitivity=='unbounded'){
      sanitized.var <- -1
      while (sanitized.var<=0){ # Make sure variance is > 0 after noise
        sanitized.var <- GaussianMechanism(tv, eps, delta, us, type.DP)
      }
    } else if (which.sensitivity=='both'){
      sanitized.var.bounded <- -1
      while (sanitized.var.bounded<=0){ # Make sure variance is > 0 after noise
        sanitized.var.bounded <- GaussianMechanism(tv, eps, delta, bs, type.DP)
      }
      sanitized.var.unbounded <- -1
      while (sanitized.var.unbounded<=0){ # Make sure variance is > 0 after noise
        sanitized.var.unbounded <- GaussianMechanism(tv, eps, delta, us, type.DP)
      }
    }
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    sanitized.var <- list()
    sanitized.var[["Bounded"]] <- sanitized.var.bounded
    sanitized.var[["Unbounded"]] <- sanitized.var.unbounded
    class(sanitized.var) <- "Sanitized Pooled Variance"
  }
  return(sanitized.var)
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
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   delta)-differential privacy, but may not do so collectively and in
#'   composition. Care must be taken not to violate differential privacy in this
#'   case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'Laplace',
#'   'Gaussian'}. Default is Laplace. See \code{\link{LaplaceMechanism}} and
#'   \code{\link{GaussianMechanism}} for a description of the supported
#'   mechanisms.
#' @param delta Nonnegative real number defining the delta privacy parameter. If
#'   0 (default), reduces to eps-DP and the Laplace mechanism is used.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism (if selected). Can be either 'pDP' for probabilistic
#'   DP \insertCite{Machanavajjhala2008}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param approx.n.max Logical indicating whether to approximate n.max (defined
#'   to be the length of the largest input vector) in the computation of the
#'   global sensitivity based on the upper and lower bounds of the data
#'   \insertCite{Liu2019b}{DPpack}. Approximation is best if n.max is very
#'   large.
#' @return Sanitized pooled covariance based on the bounded and/or unbounded
#'   definitions of differential privacy.
#' @examples
#' # Build datasets
#' D1 <- sort(stats::rnorm(500, mean=3, sd=2))
#' D2 <- sort(stats::rnorm(500, mean=-1, sd=0.5))
#' D3 <- sort(stats::rnorm(200, mean=3, sd=2))
#' D4 <- sort(stats::rnorm(200, mean=-1, sd=0.5))
#' M1 <- matrix(c(D1, D2), ncol=2)
#' M2 <- matrix(c(D3, D4), ncol=2)
#'
#' lb1 <- -3 # 3 std devs below mean
#' lb2 <- -2.5 # 3 std devs below mean
#' ub1 <- 9 # 3 std devs above mean
#' ub2 <- .5 # 3 std devs above mean
#' # Pooled covariance satisfying (1,0)-differential privacy
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
#'   \insertRef{Machanavajjhala2008}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#'   \insertRef{Liu2019b}{DPpack}
#'
#' @export
pooledCovDP <- function(..., eps=1, lower.bound1, upper.bound1, lower.bound2,
                        upper.bound2, which.sensitivity='bounded',
                        mechanism='Laplace', delta=0, type.DP='aDP',
                        approx.n.max=FALSE){
  samples <- list(...)
  #### INPUT CHECKING ####
  {
  J = length(samples)
  if (length(lower.bound1)!=1) stop("Length of lower.bound1 must be 1.")
  if (length(lower.bound2)!=1) stop("Length of lower.bound2 must be 1.")
  if (length(upper.bound1)!=1) stop("Length of upper.bound1 must be 1.")
  if (length(upper.bound2)!=1) stop("Length of upper.bound2 must be 1.")

  for (j in 1:J){
    samples[[j]][samples[[j]][,1]<lower.bound1,1] <- lower.bound1
    samples[[j]][samples[[j]][,1]>upper.bound1,1] <- upper.bound1
    samples[[j]][samples[[j]][,2]<lower.bound2,2] <- lower.bound2
    samples[[j]][samples[[j]][,2]>upper.bound2,2] <- upper.bound2
  }

  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }

  if (delta==0 & mechanism=='Gaussian') mechanism <- 'Laplace'
  if (mechanism!='Laplace' & mechanism!='Gaussian'){
    stop("Mechanism must be one of {'Laplace', 'Gaussian'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- pooledCovDataAccess(samples, lower.bound1=lower.bound1,
                                 upper.bound1=upper.bound1,
                                 lower.bound2=lower.bound2,
                                 upper.bound2=upper.bound2,
                                 approx.n.max=approx.n.max)
  tv <- results$True.Values
  bs <- results$Bounded.Sensitivities
  us <- results$Unbounded.Sensitivities
  ##########

  ########## Privacy layer
  if (mechanism=='Laplace'){
    if (which.sensitivity=='bounded'){
      sanitized.cov <- LaplaceMechanism(tv, eps, bs)
    } else if (which.sensitivity=='unbounded'){
      sanitized.cov <- LaplaceMechanism(tv, eps, us)
    } else if (which.sensitivity=='both'){
      sanitized.cov.bounded <- LaplaceMechanism(tv, eps, bs)
      sanitized.cov.unbounded <- LaplaceMechanism(tv, eps, us)
    }
  } else if (mechanism=='Gaussian'){
    if (which.sensitivity=='bounded'){
      sanitized.cov <- GaussianMechanism(tv, eps, delta, bs, type.DP)
    } else if (which.sensitivity=='unbounded'){
      sanitized.cov <- GaussianMechanism(tv, eps, delta, us, type.DP)
    } else if (which.sensitivity=='both'){
      sanitized.cov.bounded <- GaussianMechanism(tv, eps, delta, bs, type.DP)
      sanitized.cov.unbounded <- GaussianMechanism(tv, eps, delta, us, type.DP)
    }
  }
  ##########

  ########## Postprocessing layer
  if (which.sensitivity=='both'){
    sanitized.cov <- list()
    sanitized.cov[["Bounded"]] <- sanitized.cov.bounded
    sanitized.cov[["Unbounded"]] <- sanitized.cov.unbounded
    class(sanitized.cov) <- "Sanitized Pooled Covariance"
  }
  return(sanitized.cov)
  ##########
}

#' Differentially Private Quantile
#'
#' This function computes the differentially private quantile of an input vector
#' at a user-specified privacy level of epsilon.
#'
#' @param x Numeric vector of which the quantile will be taken.
#' @param quant Real number between 0 and 1 indicating which quantile to return.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Real number giving the global or public lower bound of x.
#' @param upper.bound Real number giving the global or public upper bound of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   0)-differential privacy, but may not do so collectively and in composition.
#'   Care must be taken not to violate differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'exponential'}.
#'   See \code{\link{ExponentialMechanism}} for a description of the supported
#'   mechanisms.
#' @param uniform.sampling Boolean indicating whether to sample uniformly
#'   between sorted dataset values when returning the private quantile. If TRUE,
#'   it is possible for this function to return any number between lower.bound
#'   and upper.bound. If FALSE, only a value present in the dataset or the lower
#'   bound can be returned.
#' @return Sanitized quantile based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' D <- stats::rnorm(500)
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
  {
  if (length(upper.bound)!=1) stop("Length of upper.bound must be 1.")
  if (length(lower.bound)!=1) stop("Length of lower.bound must be 1.")
  x[x<lower.bound] <- lower.bound
  x[x>upper.bound] <- upper.bound

  if (quant<0 || quant>1) stop("quant must be between 0 and 1.")

  if (which.sensitivity!='bounded' & which.sensitivity!='unbounded' &
      which.sensitivity!='both'){
    stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}.")
  }
  if (which.sensitivity=='both'){
    warning("Sensitivity based on bounded and unbounded differential privacy is identical for this statistic. Only one value will be returned.")
  }

  if (mechanism!='exponential'){
    stop("Mechanism must be one of {'exponential'}.")
  }
  }
  ##########

  ########## Data access layer
  results <- quantileDataAccess(x, quant, lower.bound=lower.bound,
                                upper.bound=upper.bound)
  utility <- results$Utility
  sorted <- results$Sorted
  sens <- results$Sensitivity

  ##########

  ########## Privacy layer
  if (mechanism=='exponential'){
    sanitized.index <- ExponentialMechanism(utility, eps, sens,
                                            measure=diff(sorted))
    if (uniform.sampling){
      sanitized.quantile <- stats::runif(1)*
        (sorted[sanitized.index+1] - sorted[sanitized.index]) +
        sorted[sanitized.index]
    } else sanitized.quantile <- sorted[sanitized.index]
  }
  ##########

  ########## Postprocessing layer
  return(sanitized.quantile)
  ##########
}

#' Differentially Private Median
#'
#' This function computes the differentially private median of an input vector
#' at a user-specified privacy level of epsilon.
#'
#' @param x Numeric vector of which the median will be taken.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param lower.bound Real number giving the global or public lower bound of x.
#' @param upper.bound Real number giving the global or public upper bound of x.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result based on bounded definition for differential privacy. If
#'   'unbounded', returns result based on unbounded definition. If 'both',
#'   returns result based on both methods \insertCite{Kifer2011}{DPpack}. Note
#'   that if 'both' is chosen, each result individually satisfies (eps,
#'   0)-differential privacy, but may not do so collectively and in composition.
#'   Care must be taken not to violate differential privacy in this case.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. Currently the following mechanisms are supported: {'exponential'}.
#'   See \code{\link{ExponentialMechanism}} for a description of the supported
#'   mechanisms.
#' @param uniform.sampling Boolean indicating whether to sample uniformly
#'   between sorted dataset values when returning the private quantile. If TRUE,
#'   it is possible for this function to return any number between lower.bound
#'   and upper.bound. If FALSE, only a value present in the dataset or the lower
#'   bound can be returned.
#' @return Sanitized median based on the bounded and/or unbounded definitions
#'   of differential privacy.
#' @examples
#' D <- stats::rnorm(500)
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
  return(sanitized.median)
  ##########
}





