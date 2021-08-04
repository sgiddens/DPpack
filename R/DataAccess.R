#' Differentially Private Mean Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private mean. The true values are computed using
#' \code{\link[base]{mean}}, while the bounded and unbounded sensitivities are
#' computed according to the theoretical values from [TODO: CITE PAPER HERE].
#'
#' @param x Numeric vector, matrix, or data frame. Means taken over columns
#'   (when applicable).
#' @param lower.bounds Numeric vector of lower bounds on each column of x. The
#'   length of lower.bounds must match the number of columns of x (length 1 if x
#'   is a vector).
#' @param upper.bounds Numeric vector of upper bounds on each column of x. The
#'   length of upper.bounds must match the number of columns of x (length 1 if x
#'   is a vector).
#' @return List of the true mean, and the bounded and unbounded sensitivities.
#' @examples
#' meanDataAccess(c(1,4,3,2), 0, 5)
#'
#' @keywords internal
meanDataAccess <- function (x, lower.bounds, upper.bounds){
  f <- mean; # Function to evaluate over x
  if (is.null(dim(x))){
    tv <- f(x);
    n <- length(x);
  } else{
    tv <- apply(x,2,f);
    n <- nrow(x);
  }
  bs <- (upper.bounds-lower.bounds)/n;
  us <- (upper.bounds-lower.bounds)/n;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Variance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private variance. The true values are computed using
#' \code{\link[stats]{var}}, while the bounded and unbounded sensitivities are
#' computed according to the theoretical values from [TODO: CITE PAPER HERE].
#'
#' @param x Numeric vector, matrix, or data frame. Variances taken over columns
#'   (when applicable).
#' @param lower.bounds Numeric vector of lower bounds on each column of x. The
#'   length of lower.bounds must match the number of columns of x (length 1 if x
#'   is a vector).
#' @param upper.bounds Numeric vector of upper bounds on each column of x. The
#'   length of upper.bounds must match the number of columns of x (length 1 if x
#'   is a vector).
#' @return List of the true variance, and the bounded and unbounded
#'   sensitivities.
#' @examples
#' varDataAccess(c(1,4,3,2), 0, 5)
#'
#' @keywords internal
varDataAccess <- function (x, lower.bounds, upper.bounds){
  f <- var; # Function to evaluate over x
  if (is.null(dim(x))){
    tv <- f(x);
    n <- length(x);
  } else{
    tv <- apply(x,2,f);
    n <- nrow(x);
  }
  bs <- (upper.bounds-lower.bounds)^2/n;
  us <- (upper.bounds-lower.bounds)^2/n;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Covariance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private covariance. The true values are computed using
#' \code{\link[stats]{cov}}, while the bounded and unbounded sensitivities are
#' computed according to the theoretical values from [TODO: CITE PAPER HERE].
#'
#' @param x1,x2 Numeric vectors.
#' @param lower.bound1,lower.bound2 Real numbers giving the lower bounds of x1
#'   and x2, respectively.
#' @param upper.bound1,upper.bound2 Real numbers giving the upper bounds of x1
#'   and x2, respectively.
#' @return List of the true covariance, and the bounded and unbounded
#'   sensitivities.
#' @examples
#' covDataAccess(c(1,4,3,2), c(-2,-3,-4,-1), 0, 5, -5, 0)
#'
#' @keywords internal
covDataAccess <- function (x1, x2, lower.bound1, upper.bound1,
                           lower.bound2, upper.bound2){
  tv <- cov(x1,x2);
  n <- length(x1);
  bs <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/n;
  us <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/n;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Histogram Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private histogram. The true values are computed using
#' \code{\link[graphics]{hist}}, while the bounded and unbounded sensitivities are
#' computed according to the theoretical values from [TODO: CITE PAPER HERE].
#'
#' @param x Numeric vector.
#' @param breaks Identical to the argument with the same name from
#'   \code{\link[graphics]{hist}}.
#' @return List of the true histogram, and the bounded and unbounded
#'   sensitivities.
#' @examples
#' histogramDataAccess(c(1,4,3,2,3), 'Sturges')
#'
#' @keywords internal
histogramDataAccess <- function (x, breaks){
  tv <- hist(x,breaks,plot=FALSE);
  tv$density <- NULL;
  bs <- 2;
  us <- 1;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Contingency Table Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private contingency table. The true values are computed using
#' \code{\link[base]{table}}, while the bounded and unbounded sensitivities are
#' computed according to the theoretical values from [TODO: CITE PAPER HERE].
#'
#' @param x,y Vectors of data from which to create the contingency table.
#' @return List of the true contingency table, and the bounded and unbounded
#'   sensitivities.
#' @examples
#' x <- MASS::Cars93$Type;
#' y <- MASS::Cars93$Origin;
#' tableDataAccess(x,y)
#'
#' @keywords internal
tableDataAccess <- function(x, y){
  tv <- table(x,y);
  bs <- 2;
  us <- 1;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Pooled Variance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private pooled variance. The true values are computed using
#' the theoretical formula and \code{\link[stats]{var}}, while the bounded and
#' unbounded sensitivities are computed according to the theoretical values from
#' [TODO: CITE PAPER HERE].
#'
#' @param samples List of vectors from which to compute the pooled variance.
#' @param lower.bound Real number giving the lower bound of the input data.
#' @param upper.bound Real number giving the upper bound of the input data.
#' @param approx.n.max Logical indicating whether to approximate n.max, which is
#'   defined to be the length of the largest input vector. Approximation is best
#'   if n.max is very large.
#' @return List of the true pooled variance, and the bounded and unbounded
#'   sensitivities.
#' @examples
#' pooledVarDataAccess(list(c(1,4,-2,8,-6),c(1,2),c(-5,-7)),-10,10,FALSE)
#'
#' @keywords internal
pooledVarDataAccess <- function(samples, lower.bound, upper.bound,approx.n.max){
  J <- length(samples);
  n <- 0;
  n.max <- 0;
  for (j in 1:J){
    nj <- length(samples[[j]]);
    n<-n+nj;
    if (n.max<nj){
      n.max <- nj;
    }
  }

  # Compute true value
  tv <- 0;
  for (j in 1:J){
    tv <- tv + (length(samples[[j]])-1)*var(samples[[j]]);
  }
  tv <- tv/(n-J);

  if (approx.n.max){
    # Make sure I can use max here this way instead of
    #   (max(upper_bounds)-min(lower_bounds))^2
    bs <- (upper.bound-lower.bound)^2/(n-J);
    if (n==2*J){
      us <- (upper.bound-lower.bound)^2*(n-1)/
        (n*(n-2));
    }else{
      us <- (upper.bound-lower.bound)^2/(n-J);
    }
  }else{
    bs <- (upper.bound-lower.bound)^2*(1-1/n.max)/(n-J);
    if (n==2*J){
      us <- (upper.bound-lower.bound)^2*(n-1)/(n*(n-2));
    }else{
      us <- (upper.bound-lower.bound)^2*(1-1/n.max)/(n-J);
    }
  }

  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Pooled Covariance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private pooled covariance. The true values are computed using
#' the theoretical formula and \code{\link[stats]{cov}}, while the bounded and
#' unbounded sensitivities are computed according to the theoretical values from
#' [TODO: CITE PAPER HERE].
#'
#' @param samples List of two-column matrices from which to compute the pooled
#'   covariance.
#' @param lower.bound1,lower.bound2 Real numbers giving the lower bounds of the
#'   first and second columns of samples, respectively.
#' @param upper.bound1,upper.bound2 Real numbers giving the upper bounds of the
#'   first and second columns of samples, respectively.
#' @param approx.n.max Logical indicating whether to approximate n.max, which is
#'   defined to be the length of the largest input vector. Approximation is best
#'   if n.max is very large.
#' @return List of the true pooled covariance, and the bounded and unbounded
#'   sensitivities.
#' @examples
#' x1 <- matrix(c(1,4,-2,8,-6,-3),ncol=2)
#' x2 <- matrix(c(1,2,-5,7),ncol=2)
#' pooledCovDataAccess(list(x1,x2),-10,10,-10,10,FALSE)
#'
#' @keywords internal
pooledCovDataAccess <- function(samples, lower.bound1, upper.bound1,
                                lower.bound2, upper.bound2,
                                approx.n.max){
  J <- length(samples);
  n <- 0;
  n.max <- 0;
  for (j in 1:J){
    nj <- nrow(samples[[j]]);
    n<-n+nj;
    if (n.max<nj){
      n.max <- nj;
    }
  }

  # Compute true value
  tv <- 0;
  for (j in 1:J){
    tv <- tv + (nrow(samples[[j]])-1)*cov(samples[[j]][,1], samples[[j]][,2]);
  }
  tv <- tv/(n-J);

  # Compute sensitivities
  if (approx.n.max){
    bs <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/(n-J);
  }else{
    bs <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)*
      (1-1/n.max)/(n-J);
  }
  us <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)*
    (1 + (n/(4*(n-1-J))) - (1/sqrt(n-J-1)))/
    (n-J);

  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

#' Differentially Private Quantile Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private quantile. The utility vector is computed as in
#' [TODO:CITE THIS], while the bounded and unbounded sensitivities are computed
#' according to the theoretical values from [TODO:PROVE THESE].
#'
#' @param x Numeric vector.
#' @param quant Real number between 0 and 1 indicating which quantile to return.
#' @param lower.bound Real number giving the lower bound of the input data.
#' @param upper.bound Real number giving the upper bound of the input data.
#' @return List of a vector corresponding to the utility function, the sorted
#'   and clipped vector of inputs, and the bounded and unbounded sensitivities.
#' @examples
#' quantileDataAccess(c(1,1,-2,8,-6),.25,-10,10)
#'
#' @keywords internal
quantileDataAccess <- function (x, quant, lower.bound, upper.bound){
  n <- length(x);
  utility <- -abs(0:n - quant*n);

  sorted <- sort(x);
  sorted[sorted<lower.bound] <- lower.bound;
  sorted[sorted>upper.bound] <- upper.bound;
  sorted <- c(lower.bound, sorted, upper.bound);

  bs <- 1;
  us <- 1;
  return(list("Utility"=utility, "Sorted"=sorted, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}
