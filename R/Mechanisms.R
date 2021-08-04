#' Laplace Mechanism
#'
#' TODO: Add description.
#'
#' @param true.values Real number or numeric vector corresponding to the true
#'   value(s) of the desired statistic(s).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param bounded.sensitivities Real number or numeric vector corresponding to
#'   the bounded global sensitivity(-ies) of the statistic. This is defined to
#'   be the greatest amount by which the statistic could change in value by
#'   changing one entry in the dataset (i.e. the total number of elements in the
#'   dataset remain the same). This value can only be NULL if which.sensitivity
#'   is 'unbounded'.
#' @param unbounded.sensitivities Real number or numeric vector corresponding to
#'   the unbounded global sensitivity(-ies) of the statistic. This is defined to
#'   be the greatest amount by which the statistic could change in value by
#'   adding/removing one entry of the dataset (i.e. the total number of elements
#'   in the dataset increases/decreases by one). This value can only be NULL if
#'   which.sensitivity is 'bounded'.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta if relevant) to the statistics. For example, if
#'   true.values is of length two and alloc.proportions = c(.75, .25), then 75%
#'   of the privacy budget eps (and delt) is allocated to the noise computation
#'   for the first element of true.values, and the remaining 25% is allocated to
#'   the noise computation for the second element of true.values. This ensures
#'   (eps, delt)-level privacy across all computations. By default, it
#'   distributes eps and delt evenly among the calculations. Input does not need
#'   to be normalized, meaning alloc.proportions = c(3,1) produces the same
#'   result as the example above.
#' @return A list of the bounded and/or unbounded sanitized statistics,
#'   sanitized via the Laplace mechanism.
#' @examples
#' laplaceMechanism(5, 1, bounded.sensitivities=0.5,
#'   which.sensitivity='bounded')
#' laplaceMechanism(c(5,3), 1, unbounded.sensitivities=1,
#'   which.sensitivity='unbounded', alloc.proportions=c(.4,.6))
#'
#' @export
laplaceMechanism <- function (true.values, eps, bounded.sensitivities=NULL,
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
#' TODO: Add description.
#'
#' @param true.values Real number or numeric vector corresponding to the true
#'   value(s) of the desired statistic(s).
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param delt Positive real number defining the delta privacy parameter.
#' @param bounded.sensitivities Real number or numeric vector corresponding to
#'   the bounded global sensitivity(-ies) of the statistic. This is defined to
#'   be the greatest amount by which the statistic could change in value by
#'   changing one entry in the dataset (i.e. the total number of elements in the
#'   dataset remain the same). This value can only be NULL if which.sensitivity
#'   is 'unbounded'.
#' @param unbounded.sensitivities Real number or numeric vector corresponding to
#'   the unbounded global sensitivity(-ies) of the statistic. This is defined to
#'   be the greatest amount by which the statistic could change in value by
#'   adding/removing one entry of the dataset (i.e. the total number of elements
#'   in the dataset increases/decreases by one). This value can only be NULL if
#'   which.sensitivity is 'bounded'.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param type.DP String indicating the type of differential privacy desired for
#'   the Gaussian mechanism. Can be either probabilistic DP ('pDP') or
#'   approximate DP ('aDP'). If type.DP is 'pDP', noise is computed using the
#'   standard deviation given by [TODO:REFERENCE FORMULA FROM DESCRIPTION],
#'   while if type.DP is 'aDP', noise is computed using the standard deviation
#'   given by [TODO:REFERENCE FORMULAT FROM DESCRIPTION].
#' @param alloc.proportions Numeric vector giving the allocation proportions of
#'   epsilon (and delta if relevant) to the statistics. For example, if
#'   true.values is of length two and alloc.proportions = c(.75, .25), then 75%
#'   of the privacy budget eps (and delt) is allocated to the noise computation
#'   for the first element of true.values, and the remaining 25% is allocated to
#'   the noise computation for the second element of true.values. This ensures
#'   (eps, delt)-level privacy across all computations. By default, it
#'   distributes eps and delt evenly among the calculations. Input does not need
#'   to be normalized, meaning alloc.proportions = c(3,1) produces the same
#'   result as the example above.
#' @return A list of the bounded and/or unbounded sanitized statistics,
#'   sanitized via the Gaussian mechanism.
#' @examples
#' gaussianMechanism(5, 1, .5, bounded.sensitivities=0.5,
#'   which.sensitivity='bounded')
#' gaussianMechanism(c(5,3), 1, .5, unbounded.sensitivities=1,
#'   which.sensitivity='unbounded', type.DP='aDP', alloc.proportions=c(.4,.6))
#'
#' @export
gaussianMechanism <- function (true.values, eps, delt, bounded.sensitivities=NULL,
                               unbounded.sensitivities=NULL,
                               which.sensitivity='bounded',
                               type.DP='pDP', alloc.proportions=NULL){
  ### INPUT CHECKING ###
  {if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.");
  }
    if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a
                                                          scalar > 0");
    if (!is.numeric(delt) || length(delt)>1 || delt<=0) stop("delt must be a
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
  alloc.delt <- delt*alloc.proportions;

  out <- list();
  if (out.bound){
    # bounded.l2.approx <- sqrt(sum(bounded.sensitivities^2)); # Lemma 2 from Gaussian paper
    bounded.noise <- double(n);
    if (type.DP == 'pDP'){ # Equation 17 from Gaussian paper
      for (i in 1:n){
        bounded.param <- bounded.sensitivities[i]*
          (sqrt(qnorm(alloc.delt[i]/2)^2+2*alloc.eps[i])-
             qnorm(alloc.delt[i]/2))/(2*alloc.eps[i]);
        bounded.noise[i] <- rnorm(1,sd=bounded.param);
      }
    } else if (type.DP == 'aDP'){ # Equation 18 from Gaussian paper
      for (i in 1:n){
        bounded.param <- bounded.sensitivities[i]*
          (sqrt(2*log(1.25/alloc.delt[i])))/alloc.eps[i];
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
          (sqrt(qnorm(alloc.delt[i]/2)^2+2*alloc.eps[i])-
             qnorm(alloc.delt[i]/2))/(2*alloc.eps[i]);
        unbounded.noise[i] <- rnorm(1,sd=unbounded.param);
      }
    } else if (type.DP == 'aDP'){ # Equation 18 from Gaussian paper
      # unbounded.l2.approx <- sqrt(sum(unbounded.sensitivities^2)); # Lemma 2 from Gaussian paper
      for (i in 1:n){
        unbounded.param <- unbounded.sensitivities[i]*
          (sqrt(2*log(1.25/alloc.delt[i])))/alloc.eps[i];
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
#' TODO: Add description.
#'
#' @param utility Numeric vector giving the utilities of the possible values.
#' @param eps Positive real number defining the epsilon privacy budget.
#' @param bounded.sensitivities Real number corresponding to the bounded global
#'   sensitivity of the statistic. This is defined to be the greatest amount by
#'   which the statistic could change in value by changing one entry in the
#'   dataset (i.e. the total number of elements in the dataset remain the same).
#'   This value can only be NULL if which.sensitivity is 'unbounded'.
#' @param unbounded.sensitivities Real number corresponding to the unbounded
#'   global sensitivity of the statistic. This is defined to be the greatest
#'   amount by which the statistic could change in value by adding/removing one
#'   entry of the dataset (i.e. the total number of elements in the dataset
#'   increases/decreases by one). This value can only be NULL if
#'   which.sensitivity is 'bounded'.
#' @param which.sensitivity String indicating which type of sensitivity to use.
#'   Can be one of {'bounded', 'unbounded', 'both'}. If 'bounded' (default),
#'   returns result plus noise based on bounded sensitivities. If 'unbounded',
#'   returns result plus noise based on unbounded sensitivities. If 'both',
#'   returns result based on both methods. Note that if 'both' is chosen, each
#'   result individually satisfies differential privacy at level eps, but may
#'   not do so collectively. Care must be taken not to violate differential
#'   privacy in this case.
#' @param measure Optional numeric vector of scaling measures for each element
#'   of utility. Should be same size as utility. Defaults to uniform scaling.
#' @param candidates Optional vector of candidates of same size as utility. If
#'   given, the function returns the candidate at the selected index rather than
#'   the index itself.
#' @return List of indices (or values if candidates given) selected by the
#'   mechanism based on bounded and/or unbounded sensitivities.
#' @examples
#' exponentialMechanism(c(0,1,2,3,2,1,0), 1, bounded.sensitivities=1,
#'   which.sensitivity='bounded')
#' exponentialMechanism(c(0,1,2,3,2,1,0), 1, unbounded.sensitivities=.5,
#'   which.sensitivity='unbounded', measure=c(1,1,2,1,2,1,1),
#'   candidates=c('a','b','c','d','e','f','g'))
#'
#' @export
exponentialMechanism <- function (utility, eps, bounded.sensitivities=NULL,
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
