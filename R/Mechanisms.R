#' Laplace Mechanism
#'
#' This function implements the Laplace mechanism for differential privacy by
#' adding noise to the true value(s) of a function according to specified values
#' of epsilon and l1-global sensitivity(-ies). Global sensitivity calculated
#' based either on bounded or unbounded differential privacy can be used
#' \insertCite{Kifer2011}{DPpack}. If true.values is a vector, the provided
#' epsilon is divided such that epsilon-differential privacy is satisfied across
#' all function values. In the case that each element of true.values comes from
#' its own function with different corresponding sensitivities, a vector of
#' sensitivities may be provided. In this case, if desired, the user can specify
#' how to divide epsilon among the function values using alloc.proportions.
#'
#' @param true.values Real number or numeric vector corresponding to the true
#'   value(s) of the desired function.
#' @param eps Positive real number defining the epsilon privacy parameter.
#' @param sensitivities Real number or numeric vector corresponding to the
#'   l1-global sensitivity(-ies) of the function(s) generating true.values. This
#'   value must be of length 1 or of the same length as true.values. If it is of
#'   length 1 and true.values is a vector, this indicates that the given
#'   sensitivity applies simultaneously to all elements of true.values and that
#'   the privacy budget need not be allocated (alloc.proportions is unused in
#'   this case). If it is of the same length as true.values, this indicates that
#'   each element of true.values comes from its own function with different
#'   corresponding sensitivities. In this case, the l1-norm of the provided
#'   sensitivities is used to generate the Laplace noise.
#' @param alloc.proportions Optional numeric vector giving the allocation
#'   proportions of epsilon to the function values in the case of vector-valued
#'   sensitivities. For example, if sensitivities is of length two and
#'   alloc.proportions = c(.75, .25), then 75% of the privacy budget eps is
#'   allocated to the noise computation for the first element of true.values,
#'   and the remaining 25% is allocated to the noise computation for the second
#'   element of true.values. This ensures eps-level privacy across all
#'   computations. Input does not need to be normalized, meaning
#'   alloc.proportions = c(3,1) produces the same result as the example above.
#' @return Sanitized function values based on the bounded and/or unbounded
#'   definitions of differential privacy, sanitized via the Laplace mechanism.
#' @examples
#' # Simulate dataset
#' n <- 100
#' c0 <- 5 # Lower bound
#' c1 <- 10 # Upper bound
#' D1 <- stats::runif(n, c0, c1)
#' epsilon <- 1 # Privacy budget
#' sensitivity <- (c1-c0)/n
#'
#' private.mean <- LaplaceMechanism(mean(D1), epsilon, sensitivity)
#' private.mean
#'
#' # Simulate second dataset
#' d0 <- 3 # Lower bound
#' d1 <- 6 # Upper bound
#' D2 <- stats::runif(n, d0, d1)
#' D <- matrix(c(D1,D2),ncol=2)
#' sensitivities <- c((c1-c0)/n, (d1-d0)/n)
#' epsilon <- 1 # Total privacy budget for all means
#'
#' # Here, sensitivities are summed and the result is used to generate Laplace
#' # noise. This is essentially the same as allocating epsilon proportional to
#' # the corresponding sensitivity. The results satisfy 1-differential privacy.
#' private.means <- LaplaceMechanism(apply(D, 2, mean), epsilon, sensitivities)
#' private.means
#'
#' # Here, privacy budget is explicitly split so that 75% is given to the first
#' # vector element and 25% is given to the second.
#' private.means <- LaplaceMechanism(apply(D, 2, mean), epsilon, sensitivities,
#'                                   alloc.proportions = c(0.75, 0.25))
#' private.means
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#' @export
LaplaceMechanism <- function (true.values, eps, sensitivities,
                              alloc.proportions=NULL) {
  ### INPUT CHECKING ###
  {
  if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.")
  }
  if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0")
  if (length(sensitivities)!=length(true.values) & length(sensitivities)!=1){
    stop("Length of sensitivities must match length of true.values or be length 1.")
  }
  if (any(sensitivities<=0)) stop("Global sensitivities must be > 0.")

  if (!is.null(alloc.proportions)){
    if (length(alloc.proportions)!=length(sensitivities)) {
      stop("Length of alloc.proportions, if given, must match length of sensitivities.")
    }
    if (any(alloc.proportions<=0)){
      stop("Values in alloc.proportions, if given, must be > 0.")
    }
    alloc.proportions <- alloc.proportions/sum(alloc.proportions)
  }
  }
  ########
  n <- length(true.values)
  if (is.null(alloc.proportions)) {
    noise <- rmutil::rlaplace(n=n,s=sum(sensitivities)/eps)
  } else{
    noise <- double(n)
    for (i in 1:n){
      noise[i] <- rmutil::rlaplace(s=sensitivities[i]/(alloc.proportions[i]*eps))
    }
  }
  private.values <- true.values + noise
  return(private.values)
}

#' Gaussian Mechanism
#'
#' This function implements the Gaussian mechanism for differential privacy by
#' adding noise to the true value(s) of a function according to specified values
#' of epsilon, delta, and l2-global sensitivity(-ies). Global sensitivity
#' calculated based either on bounded or unbounded differential privacy can be
#' used \insertCite{Kifer2011}{DPpack}. If true.values is a vector, the provided
#' epsilon and delta are divided such that (epsilon, delta)-level differential
#' privacy is satisfied across all function values. In the case that each
#' element of true.values comes from its own function with different
#' corresponding sensitivities, a vector of sensitivities may be provided. In
#' this case, if desired, the user can specify how to divide epsilon and delta
#' among the function values using alloc.proportions.
#'
#' @param true.values Real number or numeric vector corresponding to the true
#'   value(s) of the desired function.
#' @param eps Positive real number defining the epsilon privacy parameter.
#' @param delta Positive real number defining the delta privacy parameter.
#' @param sensitivities Real number or numeric vector corresponding to the
#'   l2-global sensitivity(-ies) of the function(s) generating true.values. This
#'   value must be of length 1 or of the same length as true.values. If it is of
#'   length 1 and true.values is a vector, this indicates that the given
#'   sensitivity applies simultaneously to all elements of true.values and that
#'   the privacy budget need not be allocated (alloc.proportions is unused in
#'   this case). If it is of the same length as true.values, this indicates that
#'   each element of true.values comes from its own function with different
#'   corresponding sensitivities. In this case, the l2-norm of the provided
#'   sensitivities is used to generate the Gaussian noise.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either 'pDP' for probabilistic DP
#'   \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{Dwork2006b}{DPpack}. Note that if 'aDP' is chosen, epsilon must
#'   be strictly less than 1.
#' @param alloc.proportions Optional numeric vector giving the allocation
#'   proportions of epsilon and delta to the function values in the case of
#'   vector-valued sensitivities. For example, if sensitivities is of length two
#'   and alloc.proportions = c(.75, .25), then 75% of the privacy budget eps
#'   (and 75% of delta) is allocated to the noise computation for the first
#'   element of true.values, and the remaining 25% is allocated to the noise
#'   computation for the second element of true.values. This ensures (eps,
#'   delta)-level privacy across all computations. Input does not need to be
#'   normalized, meaning alloc.proportions = c(3,1) produces the same result as
#'   the example above.
#' @param analytic Indicates whether to use the analytic Gaussian mechanism to
#'   compute the noise scale \insertCite{Balle2018}{DPpack}. Defaults to FALSE.
#' @param tol Error tolerance for binary search used in determining the noise
#'   parameter for the analytic Gaussian mechanism. Unused if analytic is FALSE.
#'   Defaults to 1e-12.
#' @return Sanitized function values based on the bounded and/or unbounded
#'   definitions of differential privacy, sanitized via the Gaussian mechanism.
#' @examples
#' # Simulate dataset
#' n <- 100
#' c0 <- 5 # Lower bound
#' c1 <- 10 # Upper bound
#' D1 <- stats::runif(n, c0, c1)
#'
#' # Privacy budget
#' epsilon <- 0.9 # eps must be in (0, 1) for approximate differential privacy
#' delta <- 0.01
#' sensitivity <- (c1-c0)/n
#'
#' # Approximate differential privacy
#' private.mean.approx <- GaussianMechanism(mean(D1), epsilon, delta,
#'                                          sensitivity)
#' private.mean.approx
#'
#' # Probabilistic differential privacy
#' private.mean.prob <- GaussianMechanism(mean(D1), epsilon, delta, sensitivity,
#'                                        type.DP = 'pDP')
#' private.mean.prob
#'
#' # Analytic Gaussian mechanism
#' epsilon <- 1.1 # epsilon can be > 1 for analytic Gaussian mechanism
#' private.mean.analytic <- GaussianMechanism(mean(D1), epsilon, delta,
#'                                            sensitivity, analytic=TRUE)
#' private.mean.analytic
#'
#' # Simulate second dataset
#' d0 <- 3 # Lower bound
#' d1 <- 6 # Upper bound
#' D2 <- stats::runif(n, d0, d1)
#' D <- matrix(c(D1,D2),ncol=2)
#' sensitivities <- c((c1-c0)/n, (d1-d0)/n)
#' epsilon <- 0.9 # Total privacy budget for all means
#' delta <- 0.01
#'
#' # Here, sensitivities are summed and the result is used to generate Laplace
#' # noise. This is essentially the same as allocating epsilon proportional to
#' # the corresponding sensitivity. The results satisfy (0.9,0.01)-approximate
#' # differential privacy.
#' private.means <- GaussianMechanism(apply(D, 2, mean), epsilon, delta,
#'                                    sensitivities)
#' private.means
#'
#' # Here, privacy budget is explicitly split so that 75% is given to the first
#' # vector element and 25% is given to the second.
#' private.means <- GaussianMechanism(apply(D, 2, mean), epsilon, delta,
#'                                    sensitivities,
#'                                    alloc.proportions = c(0.75, 0.25))
#' private.means
#'
#' @importFrom Rdpack reprompt
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Balle2018}{DPpack}
#'
#'   \insertRef{Liu2019a}{DPpack}
#'
#'   \insertRef{Dwork2006b}{DPpack}
#'
#' @export
GaussianMechanism <- function (true.values, eps, delta, sensitivities,
                               type.DP='aDP', alloc.proportions=NULL,
                               analytic=FALSE, tol=1e-12){
  ### INPUT CHECKING ###
  {
  if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.")
  }
  if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0")
  if (!is.numeric(delta) || length(delta)>1 || delta<=0) stop("delta must be a scalar > 0")
  if (length(sensitivities)!=length(true.values) & length(sensitivities)!=1){
    stop("Length of sensitivities must match length of true.values or be length 1.")
  }
  if (any(sensitivities<=0)) stop("Global sensitivities must be > 0.")

  if (analytic && type.DP!='aDP')
    stop("Only approximate DP can be used when analytic is TRUE.")
  if (type.DP!='pDP' && type.DP!='aDP') stop("type.DP must be one of {'pDP', 'aDP'}.")
  if (!analytic && type.DP=='aDP' && eps>=1)
    stop("eps must be < 1 for aDP unless analytic is TRUE.")

  if (!is.null(alloc.proportions)){
    if (length(alloc.proportions)!=length(sensitivities)) {
      stop("Length of alloc.proportions, if given, must match length of sensitivities.")
    }
    if (any(alloc.proportions<=0)){
      stop("Values in alloc.proportions, if given, must be > 0.")
    }
    alloc.proportions <- alloc.proportions/sum(alloc.proportions)
  }
  }
  ########
  n <- length(true.values)
  if (is.null(alloc.proportions)){
    if (type.DP == 'pDP'){ # Equation 17 from Gaussian paper
      param <- sqrt(sum(sensitivities^2))*(sqrt(stats::qnorm(delta/2)^2+2*eps)-
                                             stats::qnorm(delta/2))/(2*eps)
      noise <- stats::rnorm(n, sd=param)
    } else if (type.DP == 'aDP'){ # Equation 18 from Gaussian paper
      if (analytic){
        param <- calibrateAnalyticGaussianMechanism(eps, delta,
                                                    sqrt(sum(sensitivities^2)),
                                                    tol)
      } else param <- sqrt(sum(sensitivities^2))*(sqrt(2*log(1.25/delta)))/eps
      noise <- stats::rnorm(n, sd=param)
    }
  } else{
    noise <- double(n)
    alloc.eps <- eps*alloc.proportions
    alloc.delta <- delta*alloc.proportions
    if (type.DP == 'pDP'){
      for (i in 1:n){
        param <- sensitivities[i]*
          (sqrt(stats::qnorm(alloc.delta[i]/2)^2+2*alloc.eps[i])-
             stats::qnorm(alloc.delta[i]/2))/(2*alloc.eps[i])
        noise[i] <- stats::rnorm(n=1, sd=param)
      }
    } else if (type.DP == 'aDP'){
      for (i in 1:n){
        if (analytic){
          param <- calibrateAnalyticGaussianMechanism(alloc.eps[i],
                                                      alloc.delta[i],
                                                      sensitivities[i], tol)
        } else{
          param <- sensitivities[i]*(sqrt(2*log(1.25/alloc.delta[i])))/alloc.eps[i]
        }
        noise[i] <- stats::rnorm(n=1, sd=param)
      }
    }
  }
  private.values <- true.values + noise
  return(private.values)
}

#' Exponential Mechanism
#'
#' This function implements the exponential mechanism for differential privacy
#' by selecting the index of a vector of candidates to return according to a
#' user-specified vector of utility function values, epsilon, and global
#' sensitivity. Sensitivity calculated based either on bounded or unbounded
#' differential privacy can be used \insertCite{Kifer2011}{DPpack}. If measure
#' is provided, the probabilities of selecting each value are scaled according
#' to the values in measure. If candidates is provided, the function returns the
#' value of candidates at the selected index, rather than the index itself.
#'
#' @param utility Numeric vector giving the utilities of the possible values.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param sensitivity Real number corresponding to the l1-global sensitivity of
#'   the function generating utility.
#' @param measure Optional numeric vector of scaling measures for the
#'   probabilities of selecting each value. Should be same size as utility.
#'   Defaults to uniform scaling.
#' @param candidates Optional vector of candidates of same size as utility. If
#'   given, the function returns the candidate at the selected index rather than
#'   the index itself.
#' @return Indices (or values if candidates given) selected by the mechanism
#'   based on the bounded and/or unbounded definitions of differential privacy.
#' @examples
#' candidates <- c('a','b','c','d','e','f','g')
#' # Release index
#' idx <- ExponentialMechanism(c(0,1,2,3,2,1,0), 1, 1)
#' candidates[idx] # Randomly chosen candidate
#'
#' # Release candidate
#' ExponentialMechanism(c(0,1,2,3,2,1,0), 1, .5, measure=c(1,1,2,1,2,1,1),
#'   candidates=candidates)
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{McSherry2007}{DPpack}
#'
#' @export
ExponentialMechanism <- function (utility, eps, sensitivity, measure=NULL,
                                  candidates=NULL){
  ### INPUT CHECKING ###
  {
  if (!is.numeric(utility) || !is.atomic(utility)){
    stop("utility must be a numeric atomic vector or scalar.")
  }
  if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0")
  if (length(sensitivity)!=1) stop("Length of sensitivity must be 1.")
  if (sensitivity<=0) stop("sensitivity must be > 0.")

  n <- length(utility)
  if (is.null(measure)) measure <- rep(1,n)
  if (length(measure)!=n) stop("Length of measure must match length of utility.")
  if (any(measure<0)) stop("Values in measure cannot be negative.")

  if (!is.null(candidates)){
    if (length(candidates)!=n) stop("Length of candidates must match length of utility")
  }
  }
  ########
  utility <- utility - max(utility)
  probabilities <- exp(eps*utility/(2*sensitivity))
  probabilities <- probabilities * measure
  probabilities <- probabilities/sum(probabilities)
  selected <- which.max(stats::runif(1)<=cumsum(probabilities))
  if (!is.null(candidates)) selected <- candidates[selected]
  return(selected)
}
