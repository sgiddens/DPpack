sigmoid <- function(z){1/(1+exp(-z))}

h.sigmoid <- function(X, coeff){
  sigmoid(X%*%coeff)
}

h.gr.sigmoid <- function(X, coeff){
  as.numeric(sigmoid(X%*%coeff)*(1-sigmoid(X%*%coeff)))*t(X)
}

h.linear <- function(X, coeff){
  X%*%coeff
}

h.gr.linear <- function(X, coeff){
  t(X)
}

loss.cross.entropy <- function(y.hat,y){
  -(y*log(y.hat) + (1-y)*log(1-y.hat))
}

loss.gr.cross.entropy <- function(y.hat,y){
  -y/y.hat + (1-y)/(1-y.hat)
}

loss.squared.error <- function(y.hat,y){
  (y.hat-y)^2/2
}

loss.gr.squared.error <- function(y.hat,y){
  y.hat-y
}

regularizer.l2 <- function(coeff){
  coeff%*%coeff/2
}

regularizer.gr.l2 <- function(coeff){
  coeff
}

domain.linear <- list("constraints"=function(coeff) coeff%*%coeff-length(coeff),
                      "jacobian"=function(coeff) 2*coeff)


EmpiricalRiskMinimizationDP.CMS <- R6::R6Class("EmpiricalRiskMinimizationDP.CMS",
  public=list(
  h = NULL,
  h.gr = NULL,
  loss = NULL,
  loss.gr = NULL,
  regularizer = NULL,
  regularizer.gr = NULL,
  lambda = NULL,
  eps = NULL,
  c = NULL,
  coeff = NULL,
  initialize = function(h, loss, regularizer, eps, lambda, c, h.gr = NULL,
                        loss.gr = NULL, regularizer.gr = NULL){
    self$h <- h # Must be of form h(X, coeff)
    self$h.gr <- h.gr
    self$loss <- loss # Must be of form loss(y.hat,y)
    self$loss.gr <- loss.gr
    if (is.character(regularizer)){
      if (regularizer == 'l2') {
        self$regularizer <- regularizer.l2
        self$regularizer.gr <- regularizer.gr.l2
      }
      else stop("regularizer must be one of {'l2'}, or a function.")
    }
    else {
      self$regularizer <- regularizer # Must be of form regularizer(coeff)
      self$regularizer.gr <- regularizer.gr
    }
    self$eps <- eps
    self$lambda <- lambda
    self$c <- c
  },
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
    # Make sure values are correct
    if (is.null(upper.bounds)){
      warning(paste("Upper bounds missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      upper.bounds <- c(apply(X,2,max),max(y));
    } else{
      if (length(upper.bounds)!=(ncol(X)+1)) stop("Length of upper.bounds
                  must be equal to the number of columns of X plus 1 (for y).");
    }
    if (is.null(lower.bounds)){
      warning(paste("Lower bounds missing and will be calculated from the data.",
                    "This may represent additional privacy loss."));
      lower.bounds <- c(apply(X,2,min),min(y));
    } else{
      if (length(lower.bounds)!=(ncol(X)+1)) stop("Length of lower.bounds
                  must be equal to the number of columns of X plus 1 (for y).");
    }

    # Add bias if needed
    if (add.bias){
      X <- dplyr::mutate(X, bias =1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
      upper.bounds <- c(1,upper.bounds)
      lower.bounds <- c(1,lower.bounds)
    }

    # Make matrices for multiplication purposes
    X <- as.matrix(X,nrow=length(y))
    y <- as.matrix(y,ncol=1)

    # Clip based on provided bounds
    for (i in 1:(length(upper.bounds)-1)){
      X[X[,i]>upper.bounds[i],i] <- upper.bounds[i]
      X[X[,i]<lower.bounds[i],i] <- lower.bounds[i]
    }
    y[y>upper.bounds[length(upper.bounds)]] <- upper.bounds[length(upper.bounds)]
    y[y<lower.bounds[length(lower.bounds)]] <- lower.bounds[length(lower.bounds)]
    list(X=X, y=y, upper.bounds=upper.bounds, lower.bounds=lower.bounds)
  },
  postprocess_coeff = function(coeff, preprocess){
    coeff
  },
  fit = function(X, y, upper.bounds=NULL, lower.bounds=NULL, add.bias=FALSE){
    preprocess <- self$preprocess_data(X,y,upper.bounds,lower.bounds,add.bias)
    X <- preprocess$X
    y <- preprocess$y

    n <- length(y)
    d <- ncol(X)
    if (!is.infinite(eps)){
      eps.prime <- self$eps - log(1 + 2*self$c/(n*self$lambda) +
                                    self$c^2/(n^2*self$lambda^2))
      if (eps.prime > 0) {
        Delta <- 0
      }
      else {
        Delta <- self$c/(n*(exp(self$eps/4) - 1)) - self$lambda
        eps.prime <- self$eps/2
      }
      beta <- eps.prime/2
      norm.b <- rgamma(1, d, rate=beta)
      direction.b <- rnorm(d);
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b * direction.b
    } else {
      Delta <- 0
      b <- numeric(d)
    }

    tmp.coeff <- self$optimize_coeff(X, y, self$lambda, Delta, b)
    self$coeff <- self$postprocess_coeff(tmp.coeff, preprocess)
    invisible(self)
  },
  predict = function(X, add.bias=FALSE){
    if (add.bias){
      X <- dplyr::mutate(X, bias =1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
    }
    self$h(as.matrix(X), self$coeff)
  },
  optimize_coeff=function(X, y, lambda, Delta, b){
    n <- length(y)
    d <- ncol(X)

    # Get objective function
    objective <- function(par, X, y, lambda, Delta, b){
      as.numeric(sum(self$loss(self$h(X,par),y))/n +
              lambda*self$regularizer(par)/n + t(b)%*%par/n + Delta*par%*%par/2)
    }

    # Get gradient function
    if (!is.null(self$h.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      objective.gr <- function(par, X, y, lambda, Delta, b){
        as.numeric(self$h.gr(X,par)%*%self$loss.gr(self$h(X,par),y)/n +
                     lambda*self$regularizer.gr(par)/n + b/n + Delta*par)
      }
    }
    else objective.gr <- NULL

    # Run optimization
    coeff0 <- numeric(ncol(X))
    opt.res <- optim(coeff0, fn=objective, gr=objective.gr, #method="BFGS",
                     X=X, y=y, lambda=lambda, Delta=Delta, b=b)
    opt.res$par
  }
))

LogisticRegressionDP <- R6::R6Class("LogisticRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.CMS,
  public=list(
    initialize = function(regularizer, eps, lambda, regularizer.gr = NULL){
      super$initialize(h.sigmoid, loss.cross.entropy, regularizer, eps,
                      lambda, 1/4, h.gr.sigmoid, loss.gr.cross.entropy,
                      regularizer.gr)
    },
    preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
      res <- super$preprocess_data(X,y,upper.bounds,lower.bounds, add.bias)
      X <- res$X
      y <- res$y
      upper.bounds <- res$upper.bounds
      lower.bounds <- res$lower.bounds

      # Process X
      lb <- lower.bounds[1:(length(lower.bounds)-1)]
      ub <- upper.bounds[1:(length(upper.bounds)-1)]
      tmp <- matrix(c(lb,ub),ncol=2)
      scale1 <- apply(abs(tmp),1,max)
      scale2 <- sqrt(sum(scale1^2))
      X.norm <- t(t(X)/scale1)/scale2

      # Process y (not necessary since y expected to be 0 or 1)

      list(X=X.norm, y=y, scale1=scale1, scale2=scale2)
    },
    postprocess_coeff = function(coeff, preprocess){
      coeff/(preprocess$scale1*preprocess$scale2)
    },
    predict = function(X, add.bias=FALSE){
      round(super$predict(X, add.bias))
    },
    optimize_coeff = function(X, y, lambda, Delta, b){
      # Repeated implementation to avoid numerical gradient issues
      n <- length(y)
      d <- ncol(X)

      # Get objective function
      objective <- function(par, X, y, lambda, Delta, b){
        as.numeric(sum(self$loss(self$h(X,par),y))/n +
                     lambda*self$regularizer(par)/n + t(b)%*%par/n +
                     Delta*par%*%par/2)
      }

      # Get gradient function
      if (!is.null(self$h.gr) && !is.null(self$loss.gr) &&
          !is.null(self$regularizer.gr)) {
        objective.gr <- function(par, X, y, lambda, Delta, b){
          as.numeric(t(X)%*%(self$h(X, par)-y)/n +
                       lambda*self$regularizer.gr(par)/n + b/n + Delta*par)
        }
      }
      else objective.gr <- NULL

      # Run optimization
      coeff0 <- numeric(ncol(X))
      opt.res <- optim(coeff0, fn=objective, gr=objective.gr, method="BFGS",
                       X=X, y=y, lambda=lambda, Delta=Delta, b=b)
      opt.res$par
    }
))

SupportVectorMachineDP <- R6::R6Class("SupportVectorMachineDP",
  inherit=EmpiricalRiskMinimizationDP.CMS,
  public=list(

  ))

EmpiricalRiskMinimizationDP.KST <- R6::R6Class("EmpiricalRiskMinimizationDP.KST",
  public=list(
    h = NULL,
    h.gr = NULL,
    loss = NULL,
    loss.gr = NULL,
    regularizer = NULL,
    regularizer.gr = NULL,
    eps = NULL,
    delt = NULL,
    domain = NULL,
    zeta = NULL,
    lambda = NULL,
    gamma = NULL,
    coeff = NULL,
    initialize = function(h, loss, regularizer, eps, delt, domain, zeta, lambda,
                          gamma, h.gr=NULL, loss.gr=NULL, regularizer.gr=NULL){
      self$h <- h # Must be of form h(X, coeff)
      self$h.gr <- h.gr
      self$loss <- loss # Must be of form loss(y.hat,y)
      self$loss.gr <- loss.gr
      if (is.character(regularizer)){
        if (regularizer == 'l2') {
          self$regularizer <- regularizer.l2
          self$regularizer.gr <- regularizer.gr.l2
        }
        else stop("regularizer must be one of {'l2'}, or a function.")
      }
      else {
        self$regularizer <- regularizer # Must be of form regularizer(coeff)
        self$regularizer.gr <- regularizer.gr
      }
      self$eps <- eps
      self$delt <- delt
      self$domain <- domain # of form list("constraints"=g(x) (<=0), "jacobian"=g'(x))
      self$zeta <- zeta
      self$lambda <- lambda
      self$gamma <- gamma
    },
    preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
      # Make sure values are correct
      if (is.null(upper.bounds)){
        warning(paste("Upper bounds missing and will be calculated from the data.",
                      "This may represent additional privacy loss."));
        upper.bounds <- c(apply(X,2,max),max(y));
      } else{
        if (length(upper.bounds)!=(ncol(X)+1)) stop("Length of upper.bounds
                  must be equal to the number of columns of X plus 1 (for y).");
      }
      if (is.null(lower.bounds)){
        warning(paste("Lower bounds missing and will be calculated from the data.",
                      "This may represent additional privacy loss."));
        lower.bounds <- c(apply(X,2,min),min(y));
      } else{
        if (length(lower.bounds)!=(ncol(X)+1)) stop("Length of lower.bounds
                  must be equal to the number of columns of X plus 1 (for y).");
      }

      # Add bias if needed
      if (add.bias){
        X <- dplyr::mutate(X, bias =1)
        X <- X[, c(ncol(X), 1:(ncol(X)-1))]
        upper.bounds <- c(1,upper.bounds)
        lower.bounds <- c(1,lower.bounds)
      }

      # Make matrices for multiplication purposes
      X <- as.matrix(X,nrow=length(y))
      y <- as.matrix(y,ncol=1)

      # Clip based on provided bounds
      for (i in 1:(length(upper.bounds)-1)){
        X[X[,i]>upper.bounds[i],i] <- upper.bounds[i]
        X[X[,i]<lower.bounds[i],i] <- lower.bounds[i]
      }
      y[y>upper.bounds[length(upper.bounds)]] <- upper.bounds[length(upper.bounds)]
      y[y<lower.bounds[length(lower.bounds)]] <- lower.bounds[length(lower.bounds)]
      list(X=X, y=y, upper.bounds=upper.bounds, lower.bounds=lower.bounds)
    },
    postprocess_coeff = function(coeff, preprocess){
      coeff
    },
    fit = function(X, y, upper.bounds=NULL, lower.bounds=NULL, add.bias=FALSE){
      # Assumptions:
      # If assumptions are not met in child class, implement
      #     preprocess/postprocess data functions so they are
      preprocess <- self$preprocess_data(X,y,upper.bounds,lower.bounds,add.bias)
      X <- preprocess$X
      y <- preprocess$y

      n <- length(y)
      p <- ncol(X)
      if (is.null(self$eps) || is.infinite(self$eps)) Delta <- 0
      else Delta <- 2*self$lambda/self$eps
      if (is.null(self$delt) || self$delt==0){
        norm.b <- rgamma(1, p, rate=self$eps/(2*self$zeta))
        direction.b <- rnorm(p)
        direction.b <- direction.b/sqrt(sum(direction.b^2))
        b <- norm.b*direction.b
      } else{
        mu <- numeric(p)
        Sigma <- ((self$zeta^2*(8*log(2/self$delt)+4*self$eps))/
                    (self$eps^2))*diag(p)
        b <- MASS::mvrnorm(n=1,mu,Sigma)
      }
      if (is.null(self$eps) || is.infinite(self$eps)) b <- numeric(p)
      b <- as.matrix(b)

     tmp.coeff <- self$optimize_coeff(X, y, Delta, b)
     self$coeff <- self$postprocess_coeff(tmp.coeff, preprocess)
     invisible(self)
    },
    predict = function(X, add.bias=FALSE){
      if (add.bias){
        X <- dplyr::mutate(X, bias =1)
        X <- X[, c(ncol(X), 1:(ncol(X)-1))]
      }
      self$h(as.matrix(X), self$coeff)
    },
    optimize_coeff=function(X, y, Delta, b){
     n <- length(y)
     p <- ncol(X)

     # Get objective function
     objective <- function(par, X, y, Delta, b){
       as.numeric(sum(self$loss(self$h(X,par),y))/n +
         self$gamma*self$regularizer(par)/n + Delta*par%*%par/(2*n) +
         t(b)%*%par/n)
     }

     # g <- function(par, X, y, Delta, b) self$domain$constraints(par) # For new opt method

     # Get gradient function
     if (!is.null(self$h.gr) && !is.null(self$loss.gr) &&
         !is.null(self$regularizer.gr)) {
       objective.gr <- function(par, X, y, Delta, b){
         as.numeric(self$h.gr(X,par)%*%self$loss.gr(self$h(X,par),y)/n +
           self$gamma*self$regularizer.gr(par)/n + b/n + Delta*par/n)
       }
       # alg <- "NLOPT_LD_MMA" # For new opt method
       # g_jac <- function(par, X, y, Delta, b) self$domain$jacobian(par) # For new opt method
     }
     else {
       objective.gr <- NULL
       # alg <- "NLOPT_LN_COBYLA" # For new opt method
       # g_jac <- NULL # For new opt method
     }

     # Run optimization
     coeff0 <- numeric(ncol(X))

     # For old opt method
     opt.res <- optim(coeff0, fn=objective, gr=objective.gr, method="CG",
                      X=X, y=y, Delta=Delta, b=b)
     opt.res$par

     # For new opt method
     # opts <- list("algorithm"=alg, xtol_rel = 1e-4)
     # opt.res <- nloptr::nloptr(x0=coeff0, eval_f=objective,
     #                           eval_grad_f=objective.gr, eval_g_ineq=g,
     #                           eval_jac_g_ineq=g_jac, opts=opts,
     #                           X=X, y=y, Delta=Delta, b=b)
     # opt.res$solution
    }
  ))

LinearRegressionDP <- R6::R6Class("LinearRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.KST,
  public=list(
    initialize = function(regularizer, eps, delt, gamma, regularizer.gr=NULL){
      super$initialize(h.linear, loss.squared.error, regularizer, eps, delt,
                       domain.linear, NULL, NULL, gamma,
                       h.gr=h.gr.linear, loss.gr=loss.gr.squared.error,
                       regularizer.gr=regularizer.gr)
    },
    preprocess_data=function(X, y, upper.bounds, lower.bounds, add.bias){
      res <- super$preprocess_data(X,y,upper.bounds,lower.bounds, add.bias)
      X <- res$X
      y <- res$y
      upper.bounds <- res$upper.bounds
      lower.bounds <- res$lower.bounds

      # Set necessary values for linear regression
      p <- ncol(X)
      self$zeta <- 2*p^(3/2)
      self$lambda <- p

      # Process X
      lb <- lower.bounds[1:(length(lower.bounds)-1)]
      ub <- upper.bounds[1:(length(upper.bounds)-1)]
      tmp <- matrix(c(lb,ub),ncol=2)
      scale1.X <- apply(abs(tmp),1,max)
      scale2.X <- sqrt(sum(scale1.X^2))/sqrt(p)
      X.norm <- t(t(X)/scale1.X)/scale2.X

      # Process y
      lb.y <- lower.bounds[length(lower.bounds)]
      ub.y <- upper.bounds[length(upper.bounds)]
      shift.y <- lb.y + (ub.y - lb.y)/2 # Subtracting this centers at 0
      scale.y <- (ub.y-lb.y)/(2*p) # Dividing by this scales to max size p in abs()
      y.norm <- (y-shift.y)/scale.y

      list(X=X.norm,y=y.norm, scale1.X=scale1.X, scale2.X=scale2.X,
           shift.y=shift.y, scale.y=scale.y)
    },
    postprocess_coeff=function(coeff, preprocess){
      # Undo X.norm
      coeff <- coeff/(preprocess$scale1.X*preprocess$scale2.X)

      # Undo y.norm
      coeff <- coeff*preprocess$scale.y
      coeff[1] <- coeff[1] + preprocess$shift.y # Assumes coeff[1] is bias term
      coeff
    }
  ))

