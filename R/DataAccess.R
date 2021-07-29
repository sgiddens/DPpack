meanDataAccess <- function (x, lower.bounds, upper.bounds){
  # Performs the data access step for differentially private mean.
  #
  # x: The dataset over which the mean will be taken.
  #     It could be a subset of the full dataset. Means taken over columns.
  # lower.bounds: The list of lower bounds (minimums) of the full dataset from
  #     which each column of x was taken.
  # upper.bounds: The list of upper bounds (maximums) of the full dataset from
  #     which each column of x was taken.
  #
  # Returns:
  #     List containing the true means, the bounded sensitivities, and the
  #         unbounded sensitivities for each column of x.
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

varDataAccess <- function (x, lower.bounds, upper.bounds){
  # Performs the data access step for differentially private variance.
  #
  # x: The dataset over which the variance will be taken.
  #     It could be a subset of the full dataset. Means taken over columns.
  # lower.bounds: The list of lower bounds (minimums) of the full dataset from
  #     which each column of x was taken.
  # upper.bounds: The list of upper bounds (maximums) of the full dataset from
  #     which each column of x was taken.
  #
  # Returns:
  #     List containing the true variances, the bounded sensitivities, and the
  #         unbounded sensitivities for each column of x.
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

covDataAccess <- function (x1, x2, lower.bound1, upper.bound1,
                           lower.bound2, upper.bound2){
  # Performs the data access step for differentially private covariance.
  #
  # x1, x2: the datasets between which the covariance will be taken.
  #     Each could be a subset of the full dataset. Must be single dimensional.
  # lower.bound1: The lower bound (minimum) of the full dataset from
  #     which x1 was taken.
  # upper.bound1: The upper bound (maximum) of the full dataset from
  #     which x1 was taken.
  # lower.bound2: The lower bound (minimum) of the full dataset from
  #     which x2 was taken.
  # upper.bound2: The upper bound (maximum) of the full dataset from
  #     which x2 was taken.
  #
  # Returns:
  #     List containing the true covariance, the bounded sensitivity, and the
  #         unbounded sensitivity between x1, and x2.
  tv <- cov(x1,x2);
  n <- length(x1);
  bs <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/n;
  us <- (upper.bound1-lower.bound1)*(upper.bound2-lower.bound2)/n;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

histogramDataAccess <- function (x, breaks){
  # Performs the data access step for differentially private histogram.
  #
  # x: The dataset over which the histogram will be taken.
  #     It could be a subset of the full dataset.
  # breaks: same argument as with hist function.
  #
  # Returns:
  #     List containing the true histogram, the bounded sensitivity, and the
  #         unbounded sensitivity x.
  tv <- hist(x,breaks,plot=FALSE);
  tv$density <- NULL;
  bs <- 2;
  us <- 1;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

tableDataAccess <- function(x, y){
  # Performs the data access step for a differential private contingency table.
  #
  # x: the first set of data over which to create the contingency table.
  # y: the second set of data over which to create the contingency table.
  #
  # Returns:
  #     List containing the true histogram, the bounded sensitivity, and the
  #         unbounded sensitivity x.
  tv <- table(x,y);
  bs <- 2;
  us <- 1;
  return(list("True.Values"=tv, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}

pooledVarDataAccess <- function(samples, lower.bound, upper.bound,approx.n.max){
  # Performs the data access step of pooled variance.
  #
  # samples: the collection of samples from which to compute the pooled variance.
  # lower.bound: the lower bound (minimum) of the collections.
  # upper.bound: the upper bound (maximum) of the collections.
  # approx.n.max: Boolean indicating whether to approximate n.max. Approximation
  #     is best if n.max is very large.
  #
  # Returns:
  #     List containing the true histogram, the bounded sensitivity, and the
  #         unbounded sensitivity x.
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

pooledCovDataAccess <- function(samples, lower.bound1, upper.bound1,
                                lower.bound2, upper.bound2,
                                approx.n.max){
  # Performs the data access step of pooled variance.
  #
  # samples: the collection of samples from which to compute the pooled variance.
  # lower_bounds1: the lower bounds (minimums) of each of the x1 collections.
  #     If NULL, it is computed to be the min values of each of them.
  # upper_bounds1: the upper bounds (maximums) of each of the x2 collections.
  #     If NULL, it is computed to be the max values of each of them.
  # lower_bounds2: the lower bounds (minimums) of each of the x1 collections.
  #     If NULL, it is computed to be the min values of each of them.
  # upper_bounds2: the upper bounds (maximums) of each of the x2 collections.
  #     If NULL, it is computed to be the max values of each of them.
  # approx.n.max: Boolean indicating whether to approximate n.max. Approximation
  #     is best if n.max is very large.
  #
  # Returns:
  #     List containing the true histogram, the bounded sensitivity, and the
  #         unbounded sensitivity x.
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

quantileDataAccess <- function (x, quant, lower.bound, upper.bound){
  # Performs the data access step for differentially private quantile.
  #
  # x: The dataset over which the quantile will be taken.
  #     It could be a subset of the full dataset.
  # quant: Real number between 0 and 1 indicating which quantile to return.
  # lower.bound: The lower bounds (minimums) of the full dataset.
  # upper.bound: The upper bounds (maximums) of the full dataset.
  #
  # Returns:
  #     List containing the ordered data clipped at bounds, the bounded
  #         sensitivities, and the unbounded sensitivities for x.
  n <- length(x);
  utility <- -abs(0:n - quant*n);

  sorted <- sort(x);
  sorted[sorted<lower.bound] <- lower.bound;
  sorted[sorted>upper.bound] <- upper.bound;
  sorted <- c(lower.bound, sorted, upper.bound);

  bs <- 1;
  us <- 1/2;
  return(list("Utility"=utility, "Sorted"=sorted, "Bounded.Sensitivities"=bs,
              "Unbounded.Sensitivities"=us));
}
