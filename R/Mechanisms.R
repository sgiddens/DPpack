#' @export
laplaceMechanism <- function (true.values, eps, bounded.sensitivities=NULL, 
                              unbounded.sensitivities=NULL,
                              which.sensitivity='bounded',
                              alloc.proportions=NULL) {
  # Adds noise to true statistic via the Laplace mechanism.
  #
  # true.values: The true value(s) of the statistic. 
  # eps: The epsilon privacy parameter.
  # bounded.sensitivities: The bounded (changing one entry) global 
  #       sensitivity(ies) of the statistic. At least one of this and 
  #       'unbounded.sensitivities' must be non-NULL.
  # unbounded.sensitivities: The unbounded (adding/removing one entry) global 
  #       sensitivity(ies) of the statistic. At least one of this and 
  #       'bounded.sensitivities' must be non-NULL.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus Laplace noise based on 
  #       bounded sensitivities. If 'unbounded', returns result plus Laplace 
  #       noise based on unbounded sensitivities.
  # alloc.proportions: This gives the proportional allocation of epsilon to the 
  #       statistics. By default, it distributes eps evenly among the 
  #       calculations. Input does not need to be normalized.
  #
  # Returns: A list of the bounded and/or unbounded sanitized statistics.
  ### INPUT CHECKING ###
  {if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.");
  } 
  if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0");
  if (which.sensitivity=='bounded') {
    out.bound = TRUE; 
    out.unbound = FALSE;
    if (is.null(bounded.sensitivities)) {
      stop("Must provide bounded.sensitivities if which.sensitivity is 'bounded' or 'both'.")
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
      stop("Must provide unbounded.sensitivities if which.sensitivity is 'unbounded' or 'both'.")
    }
    if (length(unbounded.sensitivities)!=length(true.values)){
      stop("Length of unbounded.sensitivities must match length of true.values.");
    }
    if (any(unbounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
  }
  else if (which.sensitivity=='both') {
    out.bound = TRUE; 
    out.unbound = TRUE;
    if (is.null(bounded.sensitivities)) {
      stop("Must provide bounded.sensitivities if which.sensitivity is 'bounded' or 'both'.")
    }
    if (length(bounded.sensitivities)!=length(true.values)){
      stop("Length of bounded.sensitivities must match length of true.values.");
    }
    if (any(bounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
    if (is.null(unbounded.sensitivities)) {
      stop("Must provide unbounded.sensitivities if which.sensitivity is 'unbounded' or 'both'.")
    }
    if (length(unbounded.sensitivities)!=length(true.values)){
      stop("Length of unbounded.sensitivities must match length of true.values.");
    }
    if (any(unbounded.sensitivities<=0)){
      stop("Sensitivities must be > 0.");
    }
    if (all(bounded.sensitivities==unbounded.sensitivities)){
      message("Bounded and unbounded sensitivities are identical. Only bounded will be returned.")
      out.unbound = FALSE;
    }
  }
  else stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}");
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

#' @export
gaussianMechanism <- function (true.values, eps, delt, bounded.sensitivities=NULL,
                               unbounded.sensitivities=NULL, 
                               which.sensitivity='bounded',
                               type.DP='pDP', alloc.proportions=NULL){
  # Adds noise to true statistic via the Gaussian mechanism.
  #
  # true.values: The true value(s) of the statistic.
  # eps: The epsilon privacy parameter.
  # delt: The delta privacy parameter.
  # bounded.sensitivities: The bounded (changing one entry) global 
  #       sensitivity(ies) of the statistic. At least one of this and 
  #       'unbounded.sensitivities' must be non-NULL.
  # unbounded.sensitivities: The unbounded (adding/removing one entry) global 
  #       sensitivity(ies) of the statistic. At least one of this and 
  #       'bounded.sensitivities' must be non-NULL.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result plus Gaussian noise based on 
  #       bounded sensitivities. If 'unbounded', returns result plus Gaussian 
  #       noise based on unbounded sensitivities.
  # type.DP: The type of differential privacy desired. Can be either 
  #       probabilistic DP ('pDP') or approximate DP ('aDP').
  # alloc.proportions: This gives the proportional allocation of epsilon to the 
  #       statistics. By default, it distributes eps and delt evenly among the 
  #       calculations. Input does not need to be normalized.
  #
  # Returns: A list of the bounded and/or unbounded sanitized statistics.
  ### INPUT CHECKING ###
  {if (!is.numeric(true.values) || !is.atomic(true.values)){
    stop("true.values must be numeric atomic vectors or scalars.");
  } 
    if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0");
    if (!is.numeric(delt) || length(delt)>1 || delt<=0) stop("delt must be a scalar > 0");
    if (which.sensitivity=='bounded') {
      out.bound = TRUE; 
      out.unbound = FALSE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if which.sensitivity is 'bounded' or 'both'.")
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
        stop("Must provide unbounded.sensitivities if which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=length(true.values)){
        stop("Length of unbounded.sensitivities must match length of true.values.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
    }
    else if (which.sensitivity=='both') {
      out.bound = TRUE; 
      out.unbound = TRUE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if which.sensitivity is 'bounded' or 'both'.")
      }
      if (length(bounded.sensitivities)!=length(true.values)){
        stop("Length of bounded.sensitivities must match length of true.values.");
      }
      if (any(bounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (is.null(unbounded.sensitivities)) {
        stop("Must provide unbounded.sensitivities if which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=length(true.values)){
        stop("Length of unbounded.sensitivities must match length of true.values.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (all(bounded.sensitivities==unbounded.sensitivities)){
        message("Bounded and unbounded sensitivities are identical. Only bounded will be returned.")
        out.unbound = FALSE;
      }
    }
    else stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}");
    if (type.DP!='pDP' && type.DP!='aDP') stop("type.DP must be one of {'pDP', 'aDP'}.");
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

#' @export
exponentialMechanism <- function (utility, eps, bounded.sensitivities=NULL,
                               unbounded.sensitivities=NULL, 
                               which.sensitivity='bounded',measure=NULL,
                               candidates=NULL){
  # Selects an index (or value from candidates) via the exponential mechanism.
  #
  # utility: A numeric vector with the utilities of the possible values.
  # eps: The epsilon privacy parameter.
  # bounded.sensitivities: The bounded (changing one entry) global 
  #       sensitivity(ies) of the statistic. At least one of this and 
  #       'unbounded.sensitivities' must be non-NULL.
  # unbounded.sensitivities: The unbounded (adding/removing one entry) global 
  #       sensitivity(ies) of the statistic. At least one of this and 
  #       'bounded.sensitivities' must be non-NULL.
  # which.sensitivity: Can be one of {'bounded', 'unbounded', 'both'}. If
  #       'bounded' (default), returns result  based on bounded sensitivities. 
  #       If 'unbounded', returns result based on unbounded sensitivities.
  # measure: Optional numeric vector of scaling measures for each utility. 
  #       Should be same size as utility. Defaults to uniform scaling.
  # candidates: Optional vector of candidates of same size as utility. 
  #       If given, the function returns the candidate at the selected index 
  #       rather than the index itself.
  #
  # Returns: A list of the chosen indices (or values) using bounded and/or
  #       unbounded sensitivities.
  ### INPUT CHECKING ###
  {if (!is.numeric(utility) || !is.atomic(utility)){
    stop("utility must be a numeric atomic vector or scalar.");
  } 
    if (!is.numeric(eps) || length(eps)>1 || eps<=0) stop("eps must be a scalar > 0");
    if (which.sensitivity=='bounded') {
      out.bound = TRUE; 
      out.unbound = FALSE;
      if (is.null(bounded.sensitivities)) {
        stop("Must provide bounded.sensitivities if which.sensitivity is 'bounded' or 'both'.")
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
        stop("Must provide unbounded.sensitivities if which.sensitivity is 'unbounded' or 'both'.")
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
        stop("Must provide bounded.sensitivities if which.sensitivity is 'bounded' or 'both'.")
      }
      if (length(bounded.sensitivities)!=1){
        stop("Length of bounded.sensitivities cannot be greater than 1.");
      }
      if (any(bounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (is.null(unbounded.sensitivities)) {
        stop("Must provide unbounded.sensitivities if which.sensitivity is 'unbounded' or 'both'.")
      }
      if (length(unbounded.sensitivities)!=1){
        stop("Length of unbounded.sensitivities cannot be greater than 1.");
      }
      if (any(unbounded.sensitivities<=0)){
        stop("Sensitivities must be > 0.");
      }
      if (all(bounded.sensitivities==unbounded.sensitivities)){
        message("Bounded and unbounded sensitivities are identical. Only bounded will be returned.")
        out.unbound = FALSE;
      }
    }
    else stop("which.sensitivity must be one of {'bounded', 'unbounded', 'both'}");
    n <- length(utility);
    if (is.null(measure)) measure <- rep(1,n);
    if (length(measure)!=n) stop("Length of measure must match length of utility.");
    if (any(measure<0)) stop("Values in measure cannot be negative.");
    
    if (!is.null(candidates)){
      if (length(candidates)!=n) stop("Length of candidates must match length of utility");
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
