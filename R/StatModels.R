EmpiricalRiskMinimizationDP <- R6::R6Class("EmpiricalRiskMinimizationDP", list(
    loss = NULL,
    regularizer = NULL,
    lambda = NULL,
    eps = NULL,
    c = NULL,
    weights = NULL,
    initialize = function(loss, regularizer, lambda, eps, c,
                          loss.gr = NULL, regularizer.gr = NULL){
      self$loss <- loss
      self$regularizer <- regularizer
      self$lambda <- lambda
      self$eps <- eps
      self$c <- c
      self$loss.gr <- loss.gr
      self$regularizer.gr <- regularizer.gr
    }
    train = function(X, y){
      n <- length(y)
      d <- ncol(X)
      eps.prime <- eps - log(1 + 2*c/(n*lambda) + c^2/(n^2*lambda^2))
      if (eps.prime > 0) {
        Delta <- 0
      }
      else {
        Delta <- c/(n*(exp(eps/4) - 1)) - lambda
        eps.prime <- eps/2
      }
      beta <- eps.prime/2
      norm.b <- rgamma(1, d, rate=beta)
      direction.b <- rnorm(d);
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b * direction.b
      self$weights <- self$optimize_weights(X, y, Delta, b)
        # minimizer for loss(weights, X, y) + Delta*regularizer(weights) + b%*%weights/n
      invisible(self)
    },
    optimize_weights = function(X, y, Delta, b){


    },
    predict = function(X){
      # Predict based on trained weights
    }
  )
)
