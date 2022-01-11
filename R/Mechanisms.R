#' Laplace Mechanism
#'
#' This function implements the Laplace mechanism for differential privacy by
#' adding noise to the true value of a statistic according to specified values
#' of epsilon and global sensitivity. Sensitivity calculated based either on
#' bounded or unbounded differential privacy can be used
#' \insertCite{Kifer2011}{DPpack}. If true.values is a vector, the provided
#' epsilon is divided such that epsilon-level differential privacy is satisfied
#' across all statistics. If desired, the user can specify how to divide epsilon
#' among the statistics using alloc.proportions.
#'
#' @param true.values Real number or numeric vector corresponding to the true
#'   value(s) of the desired statistic(s).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param bounded.sensitivities Real number or numeric vector corresponding to
#'   the global sensitivity(-ies) of the statistic based on bounded differential
#'   privacy \insertCite{Kifer2011}{DPpack}. This is defined to be the greatest
#'   amount by which the statistic could change in value by changing one entry
#'   in the dataset (i.e. the total number of elements in the dataset remain the
#'   same). This value can only be NULL if which.sensitivity is 'unbounded'.
#' @param unbounded.sensitivities Real number or numeric vector corresponding to
#'   the global sensitivity(-ies) of the statistic based on unbounded
#'   differential privacy \insertCite{Kifer2011}{DPpack}. This is defined to be
#'   the greatest amount by which the statistic could change in value by
#'   adding/removing one entry of the dataset (i.e. the total number of elements
#'   in the dataset increases/decreases by one). This value can only be NULL if
#'   which.sensitivity is 'bounded'.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively. Care must be taken not to violate differential privacy in
#'   this case.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta if relevant) to the statistics. For example, if
#'   true.values is of length two and alloc.proportions = c(.75, .25), then 75%
#'   of the privacy budget eps (and delta) is allocated to the noise computation
#'   for the first element of true.values, and the remaining 25% is allocated to
#'   the noise computation for the second element of true.values. This ensures
#'   (eps, delta)-level privacy across all computations. By default, it
#'   distributes eps and delta evenly among the calculations. Input does not need
#'   to be normalized, meaning alloc.proportions = c(3,1) produces the same
#'   result as the example above.
#' @return A list of the sanitized statistics based on the bounded and/or
#'   unbounded definitions of differential privacy, sanitized via the Laplace
#'   mechanism.
#' @examples
#' LaplaceMechanism(5, 1, bounded.sensitivities=0.5,
#'   which.sensitivity='bounded')
#' LaplaceMechanism(c(5,3), 1, unbounded.sensitivities=1,
#'   which.sensitivity='unbounded', alloc.proportions=c(.4,.6))
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#' \insertRef{Kifer2011}{DPpack}
#'
#' @export
LaplaceMechanism <- function (true.values, eps, bounded.sensitivities=NULL,
                              unbounded.sensitivities=NULL,
                              which.sensitivity='bounded',
                              alloc.proportions=NULL) {
  ### INPUT CHECKING ###
  {if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.");
  }
  if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps
                                                        must be a scalar > 0");
  if (which.sensitivity=='bounded') {
    out.bound = TRUE;
    out.unbound = FALSE;
    if (is.null(bounded.sensitivities)) {
      stop("Must provide bounded.sensitivities if
           which.sensitivity is 'bounded' or 'both'.")
    }
    if (length(bounded.sensitivities)!=length(true.values)){
      stop("Length of bounded.sensitivities must match length of true.values.");
    }
    if (any(bounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
  }
  else if (which.sensitivity=='unbounded') {
    out.unbound = TRUE;
    out.bound = FALSE;
    if (is.null(unbounded.sensitivities)) {
      stop("Must provide unbounded.sensitivities if
           which.sensitivity is 'unbounded' or 'both'.")
    }
    if (length(unbounded.sensitivities)!=length(true.values)){
      stop("Length of unbounded.sensitivities
           must match length of true.values.");
    }
    if (any(unbounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
  }
  else if (which.sensitivity=='both') {
    out.bound = TRUE;
    out.unbound = TRUE;
    if (is.null(bounded.sensitivities)) {
      stop("Must provide bounded.sensitivities if
           which.sensitivity is 'bounded' or 'both'.")
    }
    if (length(bounded.sensitivities)!=length(true.values)){
      stop("Length of bounded.sensitivities must match length of true.values.");
    }
    if (any(bounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
    if (is.null(unbounded.sensitivities)) {
      stop("Must provide unbounded.sensitivities if
           which.sensitivity is 'unbounded' or 'both'.")
    }
    if (length(unbounded.sensitivities)!=length(true.values)){
      stop("Length of unbounded.sensitivities
           must match length of true.values.");
    }
    if (any(unbounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
    if (all(bounded.sensitivities==unbounded.sensitivities)){
      message("Bounded and unbounded sensitivities are identical.
              Only bounded will be returned.")
      out.unbound = FALSE;
    }
  }
  else stop("which.sensitivity must be one of
            {'bounded', 'unbounded', 'both'}");
  n <- length(true.values);
  if (is.null(alloc.proportions)) alloc.proportions <- rep(1,n);
  if (length(alloc.proportions)!=length(true.values)) {
    stop("Length of alloc.proportions must match length of true.values.");
  }
  if (any(alloc.proportions<=0)){
    stop("Values in alloc.proportions must be > 0.");
  }
  alloc.proportions <- alloc.proportions/sum(alloc.proportions);
  }
  ########
  n <- length(true.values);
  out <- list();
  if (out.bound){
    bounded.noise <- double(n);
    for (i in 1:n){
      bounded.param <- bounded.sensitivities[i]/(alloc.proportions[i]*eps);
      bounded.noise[i] <- rmutil::rlaplace(s=bounded.param);
    }
    bounded.sanitized <- true.values + bounded.noise;
    out[["Bounded"]] <- bounded.sanitized;
  }
  if (out.unbound){
    unbounded.noise <- double(n);
    for (i in 1:n){
      unbounded.param <- unbounded.sensitivities[i]/(alloc.proportions[i]*eps);
      unbounded.noise[i] <- rmutil::rlaplace(s=unbounded.param);
    }
    unbounded.sanitized <- true.values + unbounded.noise;
    out[["Unbounded"]] <- unbounded.sanitized;
  }
  return(out);
}

#' Gaussian Mechanism
#'
#' This function implements the Gaussian mechanism for differential privacy by
#' adding noise to the true value of a statistic according to specified values
#' of epsilon, delta, and global sensitivity. Sensitivity calculated based
#' either on bounded or unbounded differential privacy can be used
#' \insertCite{Kifer2011}{DPpack}. If true.values is a vector, the provided
#' epsilon and delta are divided such that (epsilon, delta)-level differential
#' privacy is satisfied across all statistics. If desired, the user can specify
#' how to divide epsilon and delta among the statistics using alloc.proportions.
#'
#' @param true.values Real number or numeric vector corresponding to the true
#'   value(s) of the desired statistic(s).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param delta Positive real number defining the delta privacy parameter.
#' @param bounded.sensitivities Real number or numeric vector corresponding to
#'   the global sensitivity(-ies) of the statistic based on bounded differential
#'   privacy \insertCite{Kifer2011}{DPpack}. This is defined to be the greatest
#'   amount by which the statistic could change in value by changing one entry
#'   in the dataset (i.e. the total number of elements in the dataset remain the
#'   same). This value can only be NULL if which.sensitivity is 'unbounded'.
#' @param unbounded.sensitivities Real number or numeric vector corresponding to
#'   the global sensitivity(-ies) of the statistic based on unbounded
#'   differential privacy \insertCite{Kifer2011}{DPpack}. This is defined to be
#'   the greatest amount by which the statistic could change in value by
#'   adding/removing one entry of the dataset (i.e. the total number of elements
#'   in the dataset increases/decreases by one). This value can only be NULL if
#'   which.sensitivity is 'bounded'.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively. Care must be taken not to violate differential privacy in
#'   this case.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either 'pDP' for probabilistic DP
#'   \insertCite{Liu2019a}{DPpack} or 'aDP' for approximate DP
#'   \insertCite{DPtextbook}{DPpack}.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta if relevant) to the statistics. For example, if
#'   true.values is of length two and alloc.proportions = c(.75, .25), then 75%
#'   of the privacy budget eps (and delta) is allocated to the noise computation
#'   for the first element of true.values, and the remaining 25% is allocated to
#'   the noise computation for the second element of true.values. This ensures
#'   (eps, delta)-level privacy across all computations. By default, it
#'   distributes eps and delta evenly among the calculations. Input does not
#'   need to be normalized, meaning alloc.proportions = c(3,1) produces the same
#'   result as the example above.
#' @return A list of the sanitized statistics based on the bounded and/or
#'   unbounded definitions of differential privacy, sanitized via the Gaussian
#'   mechanism.
#' @examples
#' GaussianMechanism(5, 1, .5, bounded.sensitivities=0.5,
#'   which.sensitivity='bounded')
#' GaussianMechanism(c(5,3), 1, .5, unbounded.sensitivities=1,
#'   which.sensitivity='unbounded', type.DP='aDP', alloc.proportions=c(.4,.6))
#'
#' @importFrom Rdpack reprompt
#'
#' @references \insertRef{Dwork2006a}{DPpack}
#'
#' \insertRef{Kifer2011}{DPpack}
#'
#' \insertRef{Liu2019a}{DPpack}
#'
#' \insertRef{DPtextbook}{DPpack}
#'
#' @export
GaussianMechanism <- function (true.values, eps, delta, bounded.sensitivities=NULL,
                               unbounded.sensitivities=NULL,
                               which.sensitivity='bounded',
                               type.DP='pDP', alloc.proportions=NULL){
  ### INPUT CHECKING ###
  {if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.");
  }
    if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a
                                                          scalar > 0");
    if (!is.numeric(delta) || length(delta)>1 || delta<=0) stop("delta must be a
                                                             scalar > 0");
    if (which.sensitivity=='bounded') {
      out.bound = TRUE;
      out.unbound = FALSE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if
             which.sensitivity is 'bounded' or 'both'.")
      }
      if (length(bounded.sensitivities)!=length(true.values)){
        stop("Length of bounded.sensitivities must
             match length of true.values.");
      }
      if (any(bounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
    }
    else if (which.sensitivity=='unbounded') {
      out.unbound = TRUE;
      out.bound = FALSE;
      if (is.null(unbounded.sensitivities)) {
        stop("Must provide unbounded.sensitivities if
             which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=length(true.values)){
        stop("Length of unbounded.sensitivities must
             match length of true.values.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
    }
    else if (which.sensitivity=='both') {
      out.bound = TRUE;
      out.unbound = TRUE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if
             which.sensitivity is 'bounded' or 'both'.")
      }
      if (length(bounded.sensitivities)!=length(true.values)){
        stop("Length of bounded.sensitivities must
             match length of true.values.");
      }
      if (any(bounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (is.null(unbounded.sensitivities)) {
        stop("Must provide unbounded.sensitivities if
             which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=length(true.values)){
        stop("Length of unbounded.sensitivities must
             match length of true.values.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (all(bounded.sensitivities==unbounded.sensitivities)){
        message("Bounded and unbounded sensitivities are identical.
                Only bounded will be returned.")
        out.unbound = FALSE;
      }
    }
    else stop("which.sensitivity must be one of
              {'bounded', 'unbounded', 'both'}");
    if (type.DP!='pDP' && type.DP!='aDP') stop("type.DP must be one of
                                               {'pDP', 'aDP'}.");
    n <- length(true.values);
    if (is.null(alloc.proportions)) alloc.proportions <- rep(1,n);
    if (length(alloc.proportions)!=length(true.values)) {
      stop("Length of alloc.proportions must match length of true.values.");
    }
    if (any(alloc.proportions<=0)){
      stop("Values in alloc.proportions must be > 0.");
    }
    alloc.proportions <- alloc.proportions/sum(alloc.proportions);
  }
  ########
  n <- length(true.values);
  alloc.eps <- eps*alloc.proportions;
  alloc.delta <- delta*alloc.proportions;

  out <- list();
  if (out.bound){
    # bounded.l2.approx <- sqrt(sum(bounded.sensitivities^2)); # Lemma 2 from Gaussian paper
    bounded.noise <- double(n);
    if (type.DP == 'pDP'){ # Equation 17 from Gaussian paper
      for (i in 1:n){
        bounded.param <- bounded.sensitivities[i]*
          (sqrt(qnorm(alloc.delta[i]/2)^2+2*alloc.eps[i])-
             qnorm(alloc.delta[i]/2))/(2*alloc.eps[i]);
        bounded.noise[i] <- rnorm(1,sd=bounded.param);
      }
    } else if (type.DP == 'aDP'){ # Equation 18 from Gaussian paper
      for (i in 1:n){
        bounded.param <- bounded.sensitivities[i]*
          (sqrt(2*log(1.25/alloc.delta[i])))/alloc.eps[i];
        bounded.noise[i] <- rnorm(1,sd=bounded.param);
      }
    }
    bounded.sanitized <- true.values + bounded.noise;
    out[["Bounded"]] <- bounded.sanitized;
  }
  if (out.unbound){
    # bounded.l2.approx <- sqrt(sum(bounded.sensitivities^2)); # Lemma 2 from Gaussian paper
    unbounded.noise <- double(n);
    if (type.DP == 'pDP'){ # Equation 17 from Gaussian paper
      for (i in 1:n){
        unbounded.param <- unbounded.sensitivities[i]*
          (sqrt(qnorm(alloc.delta[i]/2)^2+2*alloc.eps[i])-
             qnorm(alloc.delta[i]/2))/(2*alloc.eps[i]);
        unbounded.noise[i] <- rnorm(1,sd=unbounded.param);
      }
    } else if (type.DP == 'aDP'){ # Equation 18 from Gaussian paper
      # unbounded.l2.approx <- sqrt(sum(unbounded.sensitivities^2)); # Lemma 2 from Gaussian paper
      for (i in 1:n){
        unbounded.param <- unbounded.sensitivities[i]*
          (sqrt(2*log(1.25/alloc.delta[i])))/alloc.eps[i];
        unbounded.noise[i] <- rnorm(1,sd=unbounded.param);
      }
    }
    unbounded.sanitized <- true.values + unbounded.noise;
    out[["Unbounded"]] <- unbounded.sanitized;
  }
  return(out);
}

#' Exponential Mechanism
#'
#' This function implements the Exponential mechanism for differential privacy
#' by selecting the index of a vector of candidates to return according to a
#' user-specified vector of utility function values, epsilon, and global
#' sensitivity. Sensitivity calculated based either on bounded or unbounded
#' differential privacy can be used \insertCite{Kifer2011}{DPpack}. If measure
#' is provided, the elements of the utility vector are scaled according to the
#' values in measure. If candidates is provided, the function returns the value
#' of candidates at the selected index, rather than the index itself.
#'
#' @param utility Numeric vector giving the utilities of the possible values.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param bounded.sensitivities Real number corresponding to the global
#'   sensitivity(-ies) of the statistic based on bounded differential privacy
#'   \insertCite{Kifer2011}{DPpack}. This is defined to be the greatest amount
#'   by which the statistic could change in value by changing one entry in the
#'   dataset (i.e. the total number of elements in the dataset remain the same).
#'   This value can only be NULL if which.sensitivity is 'unbounded'.
#' @param unbounded.sensitivities Real number corresponding to the global
#'   sensitivity(-ies) of the statistic based on unbounded differential privacy
#'   \insertCite{Kifer2011}{DPpack}. This is defined to be the greatest amount
#'   by which the statistic could change in value by adding/removing one entry
#'   of the dataset (i.e. the total number of elements in the dataset
#'   increases/decreases by one). This value can only be NULL if
#'   which.sensitivity is 'bounded'.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded definition for differential
#'   privacy. If 'unbounded', returns result plus noise based on unbounded
#'   definition. If 'both', returns result based on both methods
#'   \insertCite{Kifer2011}{DPpack}. Note that if 'both' is chosen, each result
#'   individually satisfies differential privacy at level eps, but may not do so
#'   collectively. Care must be taken not to violate differential privacy in
#'   this case.
#' @param measure Optional numeric vector of scaling measures for each element
#'   of utility. Should be same size as utility. Defaults to uniform scaling.
#' @param candidates Optional vector of candidates of same size as utility. If
#'   given, the function returns the candidate at the selected index rather than
#'   the index itself.
#' @return A list of indices (or values if candidates given) selected by the
#'   mechanism based on the bounded and/or unbounded definitions of differential
#'   privacy.
#' @examples
#' ExponentialMechanism(c(0,1,2,3,2,1,0), 1, bounded.sensitivities=1,
#'   which.sensitivity='bounded')
#' ExponentialMechanism(c(0,1,2,3,2,1,0), 1, unbounded.sensitivities=.5,
#'   which.sensitivity='unbounded', measure=c(1,1,2,1,2,1,1),
#'   candidates=c('a','b','c','d','e','f','g'))
#'
#' @references
#'   \insertRef{Dwork2006a}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#' @export
ExponentialMechanism <- function (utility, eps, bounded.sensitivities=NULL,
                               unbounded.sensitivities=NULL,
                               which.sensitivity='bounded',measure=NULL,
                               candidates=NULL){
  ### INPUT CHECKING ###
  {if (!is.numeric(utility) || !is.atomic(utility)){
    stop("utility must be a numeric atomic vector or scalar.");
  }
    if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0");
    if (which.sensitivity=='bounded') {
      out.bound = TRUE;
      out.unbound = FALSE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if
             which.sensitivity is 'bounded' or 'both'.")
      }
      if (length(bounded.sensitivities)!=1){
        stop("Length of bounded.sensitivities cannot be greater than 1.");
      }
      if (any(bounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
    }
    else if (which.sensitivity=='unbounded') {
      out.unbound = TRUE;
      out.bound = FALSE;
      if (is.null(unbounded.sensitivities)) {
        stop("Must provide unbounded.sensitivities if
             which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=1){
        stop("Length of unbounded.sensitivities cannot be greater than 1.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
    }
    else if (which.sensitivity=='both') {
      out.bound = TRUE;
      out.unbound = TRUE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if
             which.sensitivity is 'bounded' or 'both'.")
      }
      if (length(bounded.sensitivities)!=1){
        stop("Length of bounded.sensitivities cannot be greater than 1.");
      }
      if (any(bounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (is.null(unbounded.sensitivities)) {
        stop("Must provide unbounded.sensitivities if
             which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=1){
        stop("Length of unbounded.sensitivities cannot be greater than 1.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (all(bounded.sensitivities==unbounded.sensitivities)){
        message("Bounded and unbounded sensitivities are identical.
                Only bounded will be returned.")
        out.unbound = FALSE;
      }
    }
    else stop("which.sensitivity must be one of
              {'bounded', 'unbounded', 'both'}");
    n <- length(utility);
    if (is.null(measure)) measure <- rep(1,n);
    if (length(measure)!=n) stop("Length of measure must
                                 match length of utility.");
    if (any(measure<0)) stop("Values in measure cannot be negative.");

    if (!is.null(candidates)){
      if (length(candidates)!=n) stop("Length of candidates must
                                      match length of utility");
    }
  }
  ########
  out <- list();
  utility <- utility - max(utility);
  if (out.bound){
    bounded.probabilities <- exp(eps*utility/(2*bounded.sensitivities));
    bounded.probabilities <- bounded.probabilities * measure;
    bounded.probabilities <- bounded.probabilities/sum(bounded.probabilities);
    bounded.sanitized <- which.max(runif(1)<=cumsum(bounded.probabilities));
    if (!is.null(candidates)) bounded.sanitized <- candidates[bounded.sanitized];
    out[["Bounded"]] <- bounded.sanitized;
  }
  if (out.unbound){
    unbounded.probabilities <- exp(eps*utility/(2*unbounded.sensitivities));
    unbounded.probabilities <- unbounded.probabilities * measure;
    unbounded.probabilities <- unbounded.probabilities/sum(unbounded.probabilities);
    unbounded.sanitized <- which.max(runif(1)<=cumsum(unbounded.probabilities));
    if (!is.null(candidates)) unbounded.sanitized <- candidates[unbounded.sanitized];
    out[["Unbounded"]] <- unbounded.sanitized;
  }
  return(out);
}
