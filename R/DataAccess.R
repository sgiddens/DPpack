#' Differentially Private Mean Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private mean. The true values are computed using
#' \code{\link[base]{mean}}, while the sensitivities are calculated based on
#' bounded and unbounded differential privacy \insertCite{Kifer2011}{DPpack}
#' according to the theoretical values \insertCite{Liu2019b}{DPpack}. For the
#' mean, the sensitivities based on bounded and unbounded differential privacy
#' are identical, so only one value is returned.
#'
#' @param x Dataset whose mean is desired.
#' @param lower.bound Scalar representing the global or public lower bound on
#'   values of x.
#' @param upper.bound Scalar representing the global or public upper bound on
#'   values of x.
#' @return List of the true mean and the sensitivity calculated based on
#'   theoretical values.
#' @examples
#' meanDataAccess(c(1,4,3,2), 0, 5)
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
meanDataAccess <- function (x, lower.bound, upper.bound){
  f <- mean # Function to evaluate over x
  tv <- f(x)
  n <- length(x)
  sens <- (upper.bound-lower.bound)/n
  return(list("True.Values"=tv, "Sensitivity"=sens))
}

#' Differentially Private Variance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private variance. The true values are computed using
#' \code{\link[stats]{var}}, while the sensitivities are calculated based on
#' bounded and unbounded differential privacy \insertCite{Kifer2011}{DPpack}
#' according to the theoretical values \insertCite{Liu2019b}{DPpack}. For the
#' variance, the sensitivities based on bounded and unbounded differential
#' privacy are identical, so only one value is returned.
#'
#' @param x Dataset whose variance is desired.
#' @param lower.bound Scalar representing the global or public lower bound on
#'   values of x.
#' @param upper.bound Scalar representing the global or public upper bound on
#'   values of x.
#' @return List of the true variance and the sensitivity calculated based on
#'   theoretical values.
#' @examples
#' varDataAccess(c(1,4,3,2), 0, 5)
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
varDataAccess <- function (x, lower.bound, upper.bound){
  f <- stats::var # Function to evaluate over x
  tv <- f(x)
  n <- length(x)
  sens <- (upper.bound-lower.bound)^2/n
  return(list("True.Values"=tv, "Sensitivity"=sens))
}

#' Differentially Private Covariance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private covariance. The true values are computed using
#' \code{\link[stats]{cov}}, while the sensitivities are calculated based on
#' bounded and unbounded differential privacy \insertCite{Kifer2011}{DPpack}
#' according to the theoretical values \insertCite{Liu2019b}{DPpack}. For the
#' covariance, the sensitivities based on bounded and unbounded differential
#' privacy are identical, so only one value is returned.
#'
#' @param x1,x2 Numeric vectors whose covariance is desired.
#' @param lower.bound1,lower.bound2 Real numbers giving the lower bounds of x1
#'   and x2, respectively.
#' @param upper.bound1,upper.bound2 Real numbers giving the upper bounds of x1
#'   and x2, respectively.
#' @return List of the true covariance and the sensitivity calculated based on
#'   theoretical values.
#' @examples
#' covDataAccess(c(1,4,3,2), c(-2,-3,-4,-1), 0, 5, -5, 0)
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#' \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
covDataAccess <- function (x1, x2, lower.bound1, upper.bound1,
                           lower.bound2, upper.bound2){
  tv <- stats::cov(x1,x2)
  n <- length(x1)
  sens <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/n
  return(list("True.Values"=tv, "Sensitivity"=sens))
}

#' Differentially Private Histogram Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private histogram. The true values are computed using
#' \code{\link[graphics]{hist}}, while the sensitivities are calculated based on
#' bounded and unbounded differential privacy \insertCite{Kifer2011}{DPpack}
#' according to the theoretical values \insertCite{Liu2019b}{DPpack}.
#'
#' @param x Numeric vector from which the histogram will be formed..
#' @param breaks Identical to the argument with the same name from
#'   \code{\link[graphics]{hist}}.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. If the 'Laplace' mechanism is chosen, l1 sensitivities are
#'   returned. If the 'Gaussian' mechanism is chosen, l2 sensitivities are
#'   returned.
#' @return List of the true histogram and the sensitivities calculated based on
#'   bounded and unbounded differential privacy.
#' @examples
#' histogramDataAccess(c(1,4,3,2,3), 'Sturges', 'Laplace')
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#'   \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
histogramDataAccess <- function (x, breaks, mechanism){
  tv <- graphics::hist(x, breaks, plot=FALSE)
  tv$density <- NULL
  if (mechanism=='Laplace'){
    bs <- 2
    us <- 1
  } else if (mechanism=='Gaussian'){
    bs <- sqrt(2)
    us <- 1
  }
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us))
}

#' Differentially Private Contingency Table Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private contingency table. The true values are computed using
#' \code{\link[base]{table}},while the sensitivities are calculated based on
#' bounded and unbounded differential privacy \insertCite{Kifer2011}{DPpack}
#' according to the theoretical values \insertCite{Liu2019b}{DPpack}.
#'
#' @param ... Vectors of data from which to create the contingency table.
#' @param mechanism String indicating which mechanism to use for differential
#'   privacy. If the 'Laplace' mechanism is chosen, l1 sensitivities are
#'   returned. If the 'Gaussian' mechanism is chosen, l2 sensitivities are
#'   returned.
#' @return List of the true contingency table and the sensitivities calculated
#'   based on bounded and unbounded differential privacy.
#' @examples
#' x <- MASS::Cars93$Type;
#' y <- MASS::Cars93$Origin;
#' tableDataAccess(x, y, mechanism='Laplace')
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#' \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
tableDataAccess <- function(..., mechanism='Laplace'){
  tv <- table(...)
  if (mechanism=='Laplace'){
    bs <- 2
    us <- 1
  } else if (mechanism=='Gaussian'){
    bs <- sqrt(2)
    us <- 1
  }
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us))
}

#' Differentially Private Pooled Variance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private pooled variance. The true values are computed using
#' the theoretical formula and \code{\link[stats]{var}}, while the sensitivities
#' are calculated based on bounded and unbounded differential privacy
#' \insertCite{Kifer2011}{DPpack} according to the theoretical values
#' \insertCite{Liu2019b}{DPpack}.
#'
#' @param samples List of vectors from which to compute the pooled variance.
#' @param lower.bound Real number giving the lower bound of the input data.
#' @param upper.bound Real number giving the upper bound of the input data.
#' @param approx.n.max Logical indicating whether to approximate n.max, which is
#'   defined to be the length of the largest input vector. Approximation is best
#'   if n.max is very large.
#' @return List of the true pooled variance and the sensitivities calculated
#'   based on bounded and unbounded differential privacy.
#' @examples
#' pooledVarDataAccess(list(c(1,4,-2,8,-6),c(1,2),c(-5,-7)),-10,10,FALSE)
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#' \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
pooledVarDataAccess <- function(samples, lower.bound, upper.bound,approx.n.max){
  J <- length(samples)
  n <- 0
  n.max <- 0
  for (j in 1:J){
    nj <- length(samples[[j]])
    n <- n+nj
    if (n.max<nj){
      n.max <- nj
    }
  }

  # Compute true value
  tv <- 0
  for (j in 1:J){
    tv <- tv + (length(samples[[j]])-1)*stats::var(samples[[j]])
  }
  tv <- tv/(n-J)

  if (approx.n.max){
    bs <- (upper.bound-lower.bound)^2/(n-J)
    if (n==2*J){
      us <- (upper.bound-lower.bound)^2*(n-1)/
        (n*(n-2))
    }else{
      us <- (upper.bound-lower.bound)^2/(n-J)
    }
  }else{
    bs <- (upper.bound-lower.bound)^2*(1-1/n.max)/(n-J)
    if (n==2*J){
      us <- (upper.bound-lower.bound)^2*(n-1)/(n*(n-2))
    }else{
      us <- (upper.bound-lower.bound)^2*(1-1/n.max)/(n-J)
    }
  }

  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us))
}

#' Differentially Private Pooled Covariance Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private pooled covariance. The true values are computed using
#' the theoretical formula and \code{\link[stats]{cov}}, while the sensitivities
#' are calculated based on bounded and unbounded differential privacy
#' \insertCite{Kifer2011}{DPpack} according to the theoretical values
#' \insertCite{Liu2019b}{DPpack}.
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
#' @return List of the true pooled covariance and the sensitivities calculated
#'   based on bounded and unbounded differential privacy.
#' @examples
#' x1 <- matrix(c(1,4,-2,8,-6,-3),ncol=2)
#' x2 <- matrix(c(1,2,-5,7),ncol=2)
#' pooledCovDataAccess(list(x1,x2),-10,10,-10,10,FALSE)
#'
#' @references \insertRef{Liu2019b}{DPpack}
#'
#' \insertRef{Kifer2011}{DPpack}
#'
#' @keywords internal
#'
#' @export
pooledCovDataAccess <- function(samples, lower.bound1, upper.bound1,
                                lower.bound2, upper.bound2,
                                approx.n.max){
  J <- length(samples)
  n <- 0
  n.max <- 0
  for (j in 1:J){
    nj <- nrow(samples[[j]])
    n<-n+nj
    if (n.max<nj){
      n.max <- nj
    }
  }

  # Compute true value
  tv <- 0
  for (j in 1:J){
    tv <- tv + (nrow(samples[[j]])-1)*stats::cov(samples[[j]][,1], samples[[j]][,2])
  }
  tv <- tv/(n-J)

  # Compute sensitivities
  if (approx.n.max){
    bs <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/(n-J)
  }else{
    bs <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)*
      (1-1/n.max)/(n-J)
  }
  us <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)*
    (1 + (n/(4*(n-1-J))) - (1/sqrt(n-J-1)))/
    (n-J)

  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us))
}

#' Differentially Private Quantile Data Access Function
#'
#' This function performs the data access step in the computation of a
#' differentially private quantile. The utility vector is computed as in
#' \insertCite{Smith2011a;textual}{DPpack}, while the sensitivities are
#' calculated based on bounded and unbounded differential privacy
#' \insertCite{Kifer2011}{DPpack} according to the theoretical values
#' \insertCite{Gillenwater2021}{DPpack}.
#'
#' @param x Numeric vector.
#' @param quant Real number between 0 and 1 indicating which quantile to return.
#' @param lower.bound Real number giving the lower bound of the input data.
#' @param upper.bound Real number giving the upper bound of the input data.
#' @return List of a vector corresponding to the utility function, the sorted
#'   and clipped vector of inputs and the sensitivity calculated based on
#'   theoretical values.
#' @examples
#' quantileDataAccess(c(1,1,-2,8,-6),.25,-10,10)
#'
#' @references \insertRef{Kifer2011}{DPpack}
#'
#'   \insertRef{Smith2011a}{DPpack}
#'
#'   \insertRef{Gillenwater2021}{DPpack}
#'
#' @keywords internal
#'
#' @export
quantileDataAccess <- function (x, quant, lower.bound, upper.bound){
  n <- length(x)
  utility <- -abs(0:n - quant*n)

  sorted <- sort(x)
  sorted[sorted<lower.bound] <- lower.bound
  sorted[sorted>upper.bound] <- upper.bound
  sorted <- c(lower.bound, sorted, upper.bound)

  sens <- 1
  return(list("Utility"=utility, "Sorted"=sorted, "Sensitivity"=sens))
}
