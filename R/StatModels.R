# Use keywords internal for these auxiliary functions
#' @keywords internal
#' @export
h.sigmoid <- function(X, coeff){
  e1071::sigmoid(X%*%coeff)
}

#' @export
h.gr.sigmoid <- function(X, coeff){
  as.numeric(e1071::dsigmoid(X%*%coeff))*t(X)
}

#' @export
h.linear <- function(X, coeff){
  X%*%coeff
}

#' @export
h.gr.linear <- function(X, coeff){
  t(X)
}

#' @export
loss.cross.entropy <- function(y.hat,y){
  -(y*log(y.hat) + (1-y)*log(1-y.hat))
}

#' @export
loss.gr.cross.entropy <- function(y.hat,y){
  -y/y.hat + (1-y)/(1-y.hat)
}

#' @export
loss.squared.error <- function(y.hat,y){
  (y.hat-y)^2/2
}

#' @export
loss.gr.squared.error <- function(y.hat,y){
  y.hat-y
}

#' @export
generate.loss.huber <- function(h){
  function(y.hat, y){
    z <- y.hat*y
    mask1 <- z>(1+h)
    mask2 <- abs(1-z)<=h
    mask3 <- z<(1-h)
    z[mask1] <- 0
    z[mask2] <- (1+h-z[mask2])^2/(4*h)
    z[mask3] <- 1-z[mask3]
    z
  }
}

#' @export
generate.loss.gr.huber <- function(h){
  function(y.hat, y){
    z <- y.hat*y
    mask1 <- z>(1+h)
    mask2 <- abs(1-z)<=h
    mask3 <- z<(1-h)
    z[mask1] <- 0
    z[mask2] <- -y[mask2]*(1+h-z[mask2])/(2*h)
    z[mask3] <- -y[mask3]
    z
  }
}

#' @export
regularizer.l2 <- function(coeff){
  coeff%*%coeff/2
}

#' @export
regularizer.gr.l2 <- function(coeff){
  coeff
}

#' Tune Model
#'
#' @export
tune_model<- function(models, X, y, upper.bounds, lower.bounds, add.bias=FALSE){
  m <- length(models)
  n <- length(y)
  # Split data into m+1 groups
  validateX <- X[seq(m+1,n,m+1),]
  validatey <- y[seq(m+1,n,m+1)]
  z <- numeric(m)
  for (i in 1:m){
    subX <- X[seq(i,n,m+1),]
    suby <- y[seq(i,n,m+1)]
    models[[i]]$fit(subX,suby,upper.bounds,lower.bounds,add.bias)
    validatey.hat <- models[[i]]$predict(validateX,add.bias,raw.value=FALSE)
    # Ensure predicted values match validation values (SVM uses (1,-1), LR uses (1,0))
    validatey[validatey>0] <- 1
    validatey[validatey<=0] <- 0
    validatey.hat[validatey.hat>0] <- 1
    validatey.hat[validatey.hat<=0] <- 0
    z[i] <- sum(validatey!=validatey.hat)
  }
  res <- exponentialMechanism(-z, models[[i]]$eps, 1, candidates=models)
  res$Bounded[[1]]
}

#' Privacy-preserving Empirical Risk Minimization for Binary Classification
#'
#' @description This class implements differentially private empirical risk
#'   minimization using the objective perturbation technique
#'   \insertCite{chaudhuri2011}{DPpack}. It is intended to be a framework for
#'   building more specific models via inheritance. See
#'   \code{\link{LogisticRegressionDP}} for an example of this type of
#'   structure.
#'
#' @details A new model object of this class accepts as inputs a hypothesis
#'   function, a loss function, a regularizer, an epsilon value for differential
#'   privacy, a lambda value that scales the regularizer, and a constant c
#'   meeting certain constraints related to the loss function. The model can
#'   then be fit with a dataset X (given as a data.frame), a set of binary
#'   labels y for each row of X, as well as upper and lower bounds on the
#'   possible values for each column of X and for y. In fitting, the model
#'   stores a vector of coefficients coeff which satisfy epsilon-level
#'   differential privacy. These can be released directly, or used in
#'   conjunction with the predict method to predict the label of new datapoints.
#'
#'   Note that in order to guarantee epsilon-level privacy for the empirical
#'   risk minimization model, certain constraints must be satisfied for the
#'   values used to construct the object, as well as for the data used to fit.
#'   Specifically, the provided loss function must be convex and doubly
#'   differentiable w.r.t. y.hat with |loss'(y.hat,y)|<=1 and
#'   |loss''(y.hat,y)|<=c for some constant c and for all possible values of
#'   y.hat and y, where y.hat is the predicted label and y is the true label.
#'   Additionally, it is assumed that if x represents a single row of the
#'   dataset X, then ||x||<=1 for all x. Note that because of this, a bias term
#'   cannot be included without appropriate scaling/preprocessing of the
#'   dataset. To ensure privacy, the add.bias argument in the $fit and $predict
#'   methods should only be utilized in subclasses within this package, not in
#'   this class.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
EmpiricalRiskMinimizationDP.CMS <- R6::R6Class("EmpiricalRiskMinimizationDP.CMS",
  public=list(
  #' @field h Hypothesis function of the form h(X, coeff), where X is a matrix
  #'   and coeff is a vector or matrix, that returns a column matrix of
  #'   predicted labels for each row of X.
  h = NULL,
  #' @field h.gr Function representing the gradient of the hypothesis function
  #'   with respect to the values in coeff and of the same form as h. Should be
  #'   given such that the ith row of the output represents the gradient of h
  #'   with respect to the ith coefficient.
  h.gr = NULL,
  #' @field loss Loss function of the form loss(y.hat, y), where y.hat and y are
  #'   matrices, that returns a matrix of the same shape as y.hat and y of loss
  #'   function values for the empirical risk minimization model with predicted
  #'   labels y.hat and true labels y.
  loss = NULL,
  #' @field loss.gr Function representing the gradient of the loss function with
  #'   respect to y.hat and of the same form as loss. Should be given such that
  #'   the ith row of the output represents the gradient of loss at the ith set
  #'   of input values.
  loss.gr = NULL,
  #' @field regularizer Regularization function. Must be of the form
  #'   regularizer(coeff), where coeff is a vector or matrix, that returns the
  #'   value of the regularizer at coeff.
  regularizer = NULL,
  #' @field regularizer.gr Function representing the gradient of the
  #'   regularization function with respect to coeff and of the same form as
  #'   regularizer. Should return a vector. If regularizer is a string, this
  #'   value is ignored.
  regularizer.gr = NULL,
  #' @field lambda Nonnegative real number representing the regularization
  #'   constant.
  lambda = NULL,
  #' @field eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  eps = NULL,
  #' @field c Positive real number denoting the upper bound on the absolute
  #'   value of the second derivative of the loss function, as required to
  #'   ensure differential privacy.
  c = NULL,
  #' @field coeff Numeric vector of coefficients for the model.
  coeff = NULL,
  # Below this is for KernelSVM
  D = NULL,
  sampling=NULL,
  phi=NULL,
  gamma=NULL,
  prefilter=NULL,
  #' @description Create a new EmpiricalRiskMinimizationDP.CMS object.
  #' @param h Hypothesis function. Must have form as given in h field
  #'   description.
  #' @param loss Loss function. Must have form as given in loss field
  #'   description. Additionally, in order to ensure differential privacy, the
  #'   function must be convex and doubly differentiable w.r.t. y.hat with
  #'   |loss'(y.hat, y)|<=1 and |loss''(y.hat, y)|<=c for some constant c and
  #'   for all possible values of y.hat and y.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   as given in regularizer field description. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param lambda Nonnegative real number representing the regularization
  #'   constant.
  #' @param c Positive real number such that |loss''(y.hat,y)|<=c for all
  #'   possible values of y.hat and y.
  #' @param h.gr Optional function representing the gradient of the hypothesis
  #'   function with respect to the values in coeff. Must have form as given in
  #'   h.gr field description. If not given, gradients are not used to compute
  #'   the coefficient values in fitting the model.
  #' @param loss.gr Optional function representing the gradient of the loss
  #'   function with respect to y.hat. Must have form as given in loss.gr field
  #'   description. If not given, gradients are not used to compute the
  #'   coefficient values in fitting the model.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularizer function function with respect to coeff. Must have form as
  #'   given in regularizer.gr field description. If not given, gradients are
  #'   not used to compute the coefficient values in fitting the model.
  #'
  #' @examples
  #' # Construct object for logistic regression
  #' h <- function(X, coeff) e1071::sigmoid(X%*%coeff)
  #' # Cross entropy loss
  #' loss <- function(y.hat,y) -(y*log(y.hat) + (1-y)*log(1-y.hat))
  #' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
  #' eps <- 1
  #' lambda <- 0.1
  #' c <- 1/4 # Required value for logistic regression
  #' h.gr <- function(X, coeff) as.numeric(e1071::dsigmoid(X%*%coeff))*t(X)
  #' loss.gr <- function(y.hat, y) -y/y.hat + (1-y)/(1-y.hat)
  #' regularizer.gr <- function(coeff) coeff
  #' ermdp <- EmpiricalRiskMinimizationDP.CMS$new(h, loss, regularizer, eps,
  #'                                              lambda, c, h.gr, loss.gr,
  #'                                              regularizer.gr)
  #'
  #' @return A new EmpiricalRiskMinimizationDP.CMS object.
  initialize = function(h, loss, regularizer, eps, lambda, c, h.gr = NULL,
                        loss.gr = NULL, regularizer.gr = NULL){
    self$h <- h
    self$h.gr <- h.gr
    self$loss <- loss
    self$loss.gr <- loss.gr
    if (is.character(regularizer)){
      if (regularizer == 'l2') {
        self$regularizer <- regularizer.l2
        self$regularizer.gr <- regularizer.gr.l2
      }
      else stop("regularizer must be one of {'l2'}, or a function.")
    }
    else {
      self$regularizer <- regularizer
      self$regularizer.gr <- regularizer.gr
    }
    self$eps <- eps
    self$lambda <- lambda
    self$c <- c
  },
  #' @description Fit the differentially private emprirical risk minimization
  #'   model. The function runs the objective perturbation algorithm
  #'   \insertCite{chaudhuri2011}{DPpack} to generate an objective function. A
  #'   numerical optimization method is then run to find optimal coefficients
  #'   for fitting the model given the training data and hyperparameters. The
  #'   built-in \code{\link{optim}} function using the "BFGS" optimization
  #'   method is used. If h.gr, loss.gr, and regularizer.gr are all given in the
  #'   construction of the object, the gradient of the objective function is
  #'   utilized by optim as well. The resulting privacy-preserving coefficients
  #'   are stored in coeff.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true labels for each row of X.
  #' @param upper.bounds Optional vector of length ncol(X)+1 giving upper bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the upper bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. If NULL (default), the values are computed to be the maximum values
  #'   in the data, which results in additional privacy loss. Any value in the
  #'   columns of X and y larger than the corresponding upper bound is clipped
  #'   at the bound.
  #' @param lower.bounds Optional vector of length ncol(X)+1 giving lower bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the lower bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. If NULL (default), the values are computed to be the minimum values
  #'   in the data, which results in additional privacy loss. Any value in the
  #'   columns of X and y smaller than the corresponding lower bound is clipped
  #'   at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE.
  #'
  #' @examples
  #' # Assume X is dataframe meeting assumptions for privacy
  #' # Assume 2 columns of X each bounded between -1 and 1
  #' # Assume y is 0 or 1 labels for each row of X
  #' # Assume ermdp is previously constructed object as in $new example
  #' upper.bounds <- c( 1, 1, 1) # Bounds for X and y
  #' lower.bounds <- c(-1,-1, 0) # Bounds for X and y
  #' ermdp$fit(X, y, upper.bounds, lower.bounds)
  #' ermdp$coeff # Gets private coefficients
  #'
  fit = function(X, y, upper.bounds=NULL, lower.bounds=NULL, add.bias=FALSE){
    preprocess <- private$preprocess_data(X,y,upper.bounds,lower.bounds,add.bias)
    X <- preprocess$X
    y <- preprocess$y

    n <- length(y)
    d <- ncol(X)
    if (!is.infinite(self$eps)){
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

    tmp.coeff <- private$optimize_coeff(X, y, Delta, b)
    self$coeff <- private$postprocess_coeff(tmp.coeff, preprocess)
    invisible(self)
  },
  #' @description Predict label(s) for given X using the fitted coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as X used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #'
  #' @examples
  #' # Assume Xtest is a new dataframe of the same form as X from fit
  #' # method example, with true labels ytest
  #' # Also assume ermdp$fit() has already been run on training data
  #' predicted.y <- ermdp$predict(Xtest) # Note these values need to be rounded
  #' n.errors <- sum(abs(round(predicted.y)-ytest))
  #'
  #' @return Matrix of predicted labels corresponding to each row of X.
  predict = function(X, add.bias=FALSE){
    if (add.bias){
      X <- dplyr::mutate(X, bias =1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
    }
    self$h(as.matrix(X), self$coeff)
  }
), private = list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Vector of length ncol(X)+1 giving upper bounds on the
  #   values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns
  #   of X. If NULL, the values are computed to be the maximum values in the
  #   data, which results in additional privacy loss. Any value in the columns
  #   of X and y larger than the corresponding upper bound is clipped at the
  #   bound.
  # param lower.bounds Vector of length ncol(X)+1 giving lower bounds on the
  #   values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns
  #   of X. If NULL, the values are computed to be the minimum values in the
  #   data, which results in additional privacy loss. Any value in the columns
  #   of X and y smaller than the corresponding lower bound is clipped at the
  #   bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving empirical risk
  #   minimization algorithm.
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

    # Add bias if needed (highly recommended to not do this due to unwanted
    #       regularization of bias term)
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
  # description Postprocess coefficients obtained by fitting the model using
  #   differential privacy to ensure they match original inputs X and y.
  #   Effectively undoes the processing done in preprocess_data.
  # param coeff Vector of coefficients obtained by fitting the model.
  # param preprocess List of values returned by preprocess_data and used to
  #   process coeff.
  # return The processed coefficient vector.
  postprocess_coeff = function(coeff, preprocess){
    coeff
  },
  # description Run numerical optimization method to find optimal coefficients
  #   for fitting model given training data and hyperparameters. This function
  #   builds the objective function based on the training data, and runs the
  #   built-in optim function using the "BFGS" optimization method. If h.gr,
  #   loss.gr, and regularizer.gr are all given in the construction of the
  #   object, the gradient of the objective function is utilized by optim as
  #   well.
  # param X Matrix of data to be fit.
  # param y Vector or matrix of true labels for each row of X.
  # param Delta Nonnegative real number set by the differentially private
  #   algorithm and required for constructing the objective function.
  # param b Numeric perturbation vector randomly drawn by the differentially
  #   private algorithm and required for constructing the objective function.
  # return Vector of fitted coefficients.
  optimize_coeff=function(X, y, Delta, b){
    n <- length(y)
    d <- ncol(X)

    # Get objective function
    objective <- function(par, X, y, Delta, b){
      as.numeric(sum(self$loss(self$h(X,par),y))/n +
                   self$lambda*self$regularizer(par)/n + t(b)%*%par/n +
                   Delta*par%*%par/2)
    }

    # Get gradient function
    if (!is.null(self$h.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      objective.gr <- function(par, X, y, Delta, b){
        as.numeric(self$h.gr(X,par)%*%self$loss.gr(self$h(X,par),y)/n +
                     self$lambda*self$regularizer.gr(par)/n + b/n + Delta*par)
      }
    }
    else objective.gr <- NULL

    # Run optimization
    coeff0 <- numeric(ncol(X))
    opt.res <- optim(coeff0, fn=objective, gr=objective.gr, method="BFGS",
                     X=X, y=y, Delta=Delta, b=b)
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
    # Specifically, it is
    #   assumed that each row of X satisfies ||X[i,]||<=1. In order to guarantee
    #   this assumption, the data are preprocessed and scaled using the upper and
    #   lower bounds on the columns of X provided.
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

      # Process y.
      # We need labels to be 0 or 1. Any nonnegative value is considered 1,
      #     while any negative or 0 value is considered 0.
      if (any((y!=1) & (y!=0))) warning("Positive values in y will be coerced
                                        to 1, and non-positive values in y
                                        will be coerced to 0.")
      y[y>0] <- 1
      y[y<=0] <- 0

      list(X=X.norm, y=y, scale1=scale1, scale2=scale2)
    },
    postprocess_coeff = function(coeff, preprocess){
      coeff/(preprocess$scale1*preprocess$scale2)
    },
    predict = function(X, add.bias=FALSE, raw.value=FALSE){
      if (raw.value) super$predict(X, add.bias)
      else round(super$predict(X, add.bias))
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
    initialize = function(regularizer, eps, lambda, regularizer.gr=NULL, huber.h=1){
      super$initialize(h.linear, generate.loss.huber(huber.h), regularizer, eps,
                       lambda, 1/(2*huber.h), h.gr.linear,
                       generate.loss.gr.huber(huber.h), regularizer.gr)
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

      # Process y.
      # We need labels to be -1 or 1. Any nonnegative value is considered 1,
      #     while any negative or 0 value is considered -1.
      if (any((y!=1) & (y!=-1))) warning("Positive values in y will be coerced
                                        to 1, and non-positive values in y
                                        will be coerced to -1.")
      y[y>0] <- 1
      y[y<=0] <- -1

      list(X=X.norm, y=y, scale1=scale1, scale2=scale2)
    },
    postprocess_coeff = function(coeff, preprocess){
      coeff/(preprocess$scale1*preprocess$scale2)
    },
    predict = function(X, add.bias=FALSE, raw.value=FALSE){
      if (raw.value) super$predict(X, add.bias)
      else sign(super$predict(X, add.bias))
    }
  ))

generate.sampling <- function(gamma){
  function(d){
    omega <- rnorm(d,sd=sqrt(2*gamma))
    phi <- runif(1,-pi,pi)
    c(omega,phi)
  }
}

phi.gaussian <- function(x, theta){
  d <- length(x)
  cos(x%*%theta[1:d] + theta[d+1])
}

KernelSupportVectorMachineDP <- R6::R6Class("KernelSupportVectorMachineDP",
  inherit=SupportVectorMachineDP,
  public=list(
    initialize = function(regularizer, eps, lambda, D, sampling="Gaussian",
                          regularizer.gr=NULL, huber.h=1, phi=NULL, gamma=1){
      super$initialize(regularizer, eps, lambda, regularizer.gr, huber.h)
      self$D <- D
      if (is.character(sampling)){
        if (sampling=="Gaussian"){
          self$sampling <- generate.sampling(gamma)
          self$phi <- phi.gaussian
        } else stop("sampling must be one of {'Gaussian'}, or a function.")
      } else {
        if (is.null(phi)) stop("phi must be specified if samping is a function.")
        self$sampling <- sampling
        self$phi <- phi
      }
    },
    preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
      # Add bias if needed (highly recommended to not do this due to unwanted
      #       regularization of bias term)
      if (add.bias){
        X <- dplyr::mutate(X, bias =1)
        X <- X[, c(ncol(X), 1:(ncol(X)-1))]
        upper.bounds <- c(1,upper.bounds)
        lower.bounds <- c(1,lower.bounds)
      }

      # Make matrices for multiplication purposes
      X <- as.matrix(X,nrow=length(y))
      y <- as.matrix(y,ncol=1)

      # Get prefilter
      d <- ncol(X)
      prefilter <- matrix(NaN, nrow=self$D, ncol=(d+1))
      for (i in 1:self$D){
        prefilter[i,] <- self$sampling(d)
      }
      self$prefilter <- prefilter

      V <- self$XtoV(X)

      lbs <- c(numeric(ncol(V)) - sqrt(1/D), -1)
      ubs <- c(numeric(ncol(V)) + sqrt(1/D), 1)
      super$preprocess_data(V, y, ubs, lbs, add.bias=FALSE)
    },
    XtoV = function(X){
      # Get V (higher dimensional X)
      V <- matrix(NaN, nrow=nrow(X), ncol=self$D)
      for (i in 1:nrow(X)){
        for (j in 1:self$D){
          V[i,j] <- sqrt(1/self$D)*self$phi(X[i,],self$prefilter[j,])
        }
      }
      V
    },
    predict = function(X, add.bias=FALSE, raw.value=FALSE){
      if (add.bias){
        X <- dplyr::mutate(X, bias =1)
        X <- X[, c(ncol(X), 1:(ncol(X)-1))]
      }
      V <- self$XtoV(as.matrix(X))
      super$predict(V, add.bias=FALSE, raw.value)
    }
  ))

#' Privacy-preserving Empirical Risk Minimization for Regression
#'
#' @description This class implements differentially private empirical risk
#'   minimization using the objective perturbation technique
#'   \insertCite{Kifer2012}{DPpack}. It is intended to be a framework for
#'   building more specific models via inheritance. See
#'   \code{\link{LinearRegressionDP}} for an example of this type of structure.
#'
#' @details A new model object of this class accepts as inputs various functions
#'   and hyperparameters related to the specific regression problem. The model
#'   can then be fit with a dataset X (given as a data.frame), a set of binary
#'   labels y for each row of X, as well as upper and lower bounds on the
#'   possible values for each column of X and for y. In fitting, the model
#'   stores a vector of coefficients coeff which satisfy epsilon-level
#'   differential privacy. These can be released directly, or used in
#'   conjunction with the predict method to predict the label of new datapoints.
#'
#'   Note that in order to guarantee epsilon-level privacy for the empirical
#'   risk minimization model, certain constraints must be satisfied for the
#'   values used to construct the object, as well as for the data used to fit.
#'   Specifically, the following constraints must be met. First, the provided
#'   domain must be a closed convex subset of R^p, where p is the number of
#'   columns of X. Second, if \eqn{L(\theta,X,y) = (1/n)\sum l(\theta,x_i,y_i)}
#'   is the loss function, then \eqn{L} must be convex with a continuous
#'   Hessian. Third, \eqn{||l(\theta,x_i,y_i)||\le} zeta for some constant zeta
#'   for all \eqn{x_i, y_i, \theta}. Fourth, for all \eqn{x_i, y_i, \theta} the
#'   Hessian of \eqn{l(\theta,x_i,y_i)} must be of rank at most one and its
#'   Eigenvalues must be bounded above by some value lambda. Finally, the
#'   regularizer must be convex.
#'
#' @references \insertRef{Kifer2012}{DPpack}
#'
#' @export
EmpiricalRiskMinimizationDP.KST <- R6::R6Class("EmpiricalRiskMinimizationDP.KST",
  public=list(
  #' @field h Hypothesis function of the form h(X, coeff), where X is a matrix
  #'   and coeff is a vector or matrix, that returns a column matrix of
  #'   predicted labels for each row of X.
  h = NULL,
  #' @field h.gr Function representing the gradient of the hypothesis function
  #'   with respect to the values in coeff and of the same form as h. Should
  #'   be given such that the ith row of the output represents the gradient of
  #'   h with respect to the ith coefficient.
  h.gr = NULL,
  #' @field loss Loss function of the form loss(y.hat, y), where y.hat and y
  #'   are matrices, that returns a matrix of the same shape as y.hat and y of
  #'   loss function values for the empirical risk minimization model with
  #'   predicted labels y.hat and true labels y.
  loss = NULL,
  #' @field loss.gr Function representing the gradient of the loss function
  #'   with respect to y.hat and of the same form as loss. Should be given
  #'   such that the ith row of the output represents the gradient of loss at
  #'   the ith set of input values.
  loss.gr = NULL,
  #' @field regularizer Regularization function. Must be of the form
  #'   regularizer(coeff), where coeff is a vector or matrix, that returns the
  #'   value of the regularizer at coeff.
  regularizer = NULL,
  #' @field regularizer.gr Function representing the gradient of the
  #'   regularization function with respect to coeff and of the same form as
  #'   regularizer. Should return a vector. If regularizer is a string, this
  #'   value is ignored.
  regularizer.gr = NULL,
  #' @field eps Positive real number defining the epsilon privacy budget. If
  #'   set to Inf, runs algorithm without differential privacy.
  eps = NULL,
  #' @field delt Nonnegative real number defining the delta parameter for
  #'   approximate differential privacy. If set to 0, pure differential
  #'   privacy is used.
  delt = NULL,
  #' @field domain List of functions representing the constraints on the
  #'   search space for the objective perturbation algorithm. Must contain two
  #'   function, labeled "constraints" and "jacobian", respectively. The
  #'   "constraints" function accepts a vector of coefficients from the search
  #'   space and returns a value such that the value is \eqn{\le 0} if and
  #'   only if the input coefficient vector is within the constrained search
  #'   space. The "jacobian" function also accepts a vector of coefficients
  #'   and returns the Jacobian of the constraint function. For example, in
  #'   linear regression, the coefficient vector \eqn{\theta} is assumed to
  #'   satisfy \eqn{||\theta||_2 \le sqrt(p)}, where \eqn{p} is the length of
  #'   \eqn{\theta} \insertCite{Kifer2012}{DPpack}. So, domain could be
  #'   defined as `domain <- list("constraints"=function(coeff)
  #'   coeff%*%coeff-length(coeff), "jacobian"=function(coeff) 2*coeff)`.
  domain = NULL,
  #' @field zeta Positive real number corresponding to the upper bound of the
  #'   2-norm of the gradient of the loss function with respect to the
  #'   coefficient vector, i.e. \eqn{||l(\theta,x_i,y_i)||\le} zeta for some
  #'   constant zeta for all \eqn{x_i, y_i, \theta}.
  zeta = NULL,
  #' @field lambda Positive real number corresponding to the upper bound of
  #'   the Eigenvalues of the Hessian of \eqn{l(\theta,x_i,y_i)} for all
  #'   \eqn{x_i, y_i, \theta}.
  lambda = NULL,
  #' @field gamma Nonnegative real number representing the regularization
  #'   constant.
  gamma = NULL,
  #' @field coeff Numeric vector of coefficients for the model.
  coeff = NULL,
  #' @description Create a new EmpiricalRiskMinimizationDP.KST object.
  #' @param h Hypothesis function. Must have form as given in h field
  #'   description.
  #' @param loss Loss function. Must have form as given in loss field
  #'   description. Additionally, to ensure differential privacy, it must be
  #'   convex, \eqn{||l(\theta,x_i,y_i)||\le} zeta for some constant zeta for
  #'   all \eqn{x_i, y_i, \theta}, and for all \eqn{x_i, y_i, \theta} the
  #'   Hessian of \eqn{l(\theta,x_i,y_i)} must be of rank at most one and its
  #'   Eigenvalues must be bounded above by some value lambda.
  #' @param regularizer String or regularization function. If a string, must
  #'   be 'l2', indicating to use l2 regularization. If a function, must have
  #'   form as given in regularizer field description. Additionally, in order
  #'   to ensure differential privacy, the function must be 1-strongly convex
  #'   and doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If
  #'   set to Inf, runs algorithm without differential privacy.
  #' @param delt Nonnegative real number defining the delta parameter for
  #'   approximate differential privacy. If set to 0, pure differential
  #'   privacy is used.
  #' @param domain List of functions representing the constraints on the
  #'   search space for the objective perturbation algorithm of form as given
  #'   in domain field description.
  #' @param zeta Positive real number corresponding to the upper bound of the
  #'   2-norm of the gradient of the loss function with respect to the
  #'   coefficient vector, i.e. \eqn{||l(\theta,x_i,y_i)||\le} zeta for some
  #'   constant zeta for all \eqn{x_i, y_i, \theta}.
  #' @param lambda Positive real number corresponding to the upper bound of
  #'   the Eigenvalues of the Hessian of \eqn{l(\theta,x_i,y_i)} for all
  #'   \eqn{x_i, y_i, \theta}.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param h.gr Optional function representing the gradient of the hypothesis
  #'   function with respect to the values in coeff. Must have form as given
  #'   in h.gr field description. If not given, gradients are not used to
  #'   compute the coefficient values in fitting the model.
  #' @param loss.gr Optional function representing the gradient of the loss
  #'   function with respect to y.hat. Must have form as given in loss.gr
  #'   field description. If not given, gradients are not used to compute the
  #'   coefficient values in fitting the model.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularizer function function with respect to coeff. Must have form as
  #'   given in regularizer.gr field description. If not given, gradients are
  #'   not used to compute the coefficient values in fitting the model.
  #'
  #' @examples
  #' # Construct object for linear regression
  #' h <- function(X, coeff) X%*%coeff
  #' loss <- function(y.hat, y) (y.hat-y)^2/2
  #' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
  #' eps <- 1
  #' delt <- 1
  #' domain <- list("constraints"=function(coeff) coeff%*%coeff-length(coeff),
  #'   "jacobian"=function(coeff) 2*coeff)
  #' # Set p to be the number of predictors desired including intercept term (length of coeff)
  #' zeta <- 2*p^(3/2)
  #' lambda <- p
  #' gamma <- 1
  #' h.gr <- function(X, coeff) t(X)
  #' loss.gr <- function(y.hat, y) y.hat-y
  #' regularizer.gr <- function(coeff) coeff
  #'
  #' ermdp <- EmpiricalRiskMinimizationDP.KST$new(h, loss, 'l2', eps, delt,
  #'                                              domain, zeta, lambda,
  #'                                              gamma, h.gr, loss.gr,
  #'                                              regularizer.gr)
  #'
  #' @return A new EmpiricalRiskMinimizationDP.KST object.
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
  #' @description Fit the differentially private emprirical risk minimization
  #'   model. The function runs the objective perturbation algorithm
  #'   \insertCite{Kifer2012}{DPpack} to generate an objective function. A
  #'   numerical constrained optimization method is then run to find optimal
  #'   coefficients for fitting the model given the training data and
  #'   hyperparameters. The \code{\link{nloptr}} function is used. If h.gr,
  #'   loss.gr, and regularizer.gr are all given in the construction of the
  #'   object, the gradient of the objective function and the Jacobian of the
  #'   constraint function are utilized for the algorithm, and the NLOPT_LD_MMA
  #'   method is used. If one or more of these gradient functions are not given,
  #'   the NLOPT_LN_COBYLA method is used. The resulting privacy-preserving
  #'   coefficients are stored in coeff.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true labels for each row of X.
  #' @param upper.bounds Optional vector of length ncol(X)+1 giving upper bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the upper bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. If NULL (default), the values are computed to be the maximum values
  #'   in the data, which results in additional privacy loss. Any value in the
  #'   columns of X and y larger than the corresponding upper bound is clipped
  #'   at the bound.
  #' @param lower.bounds Optional vector of length ncol(X)+1 giving lower bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the lower bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. If NULL (default), the values are computed to be the minimum values
  #'   in the data, which results in additional privacy loss. Any value in the
  #'   columns of X and y smaller than the corresponding lower bound is clipped
  #'   at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE.
  #'
  #' @examples
  #' # Assume X is dataframe meeting assumptions for privacy
  #' # Assume 2 columns of X, with the first being all 1 (for intercept), and
  #' #   the second being between -1 and 1
  #' # Assume y is a matrix with values between -2 and 2
  #' # Assume ermdp is previously constructed object as in $new example
  #' upper.bounds <- c( 1, 1, 2) # Bounds for X and y
  #' lower.bounds <- c(1,-1, -2) # Bounds for X and y
  #' ermdp$fit(X, y, upper.bounds, lower.bounds)
  #' ermdp$coeff # Gets private coefficients
  #'
  fit = function(X, y, upper.bounds=NULL, lower.bounds=NULL, add.bias=FALSE){
    # Assumptions:
    # If assumptions are not met in child class, implement
    #     preprocess/postprocess data functions so they are
    preprocess <- private$preprocess_data(X,y,upper.bounds,lower.bounds,add.bias)
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

   tmp.coeff <- private$optimize_coeff(X, y, Delta, b)
   self$coeff <- private$postprocess_coeff(tmp.coeff, preprocess)
   invisible(self)
  },
  #' @description Predict y values for given X using the fitted coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as X used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #'
  #' @examples
  #' # Assume Xtest is a new dataframe of the same form as X from fit
  #' # method example, with true labels ytest
  #' # Also assume ermdp$fit() has already been run on training data
  #' predicted.y <- ermdp$predict(Xtest)
  #'
  #' @return Matrix of predicted y values corresponding to each row of X.
  predict = function(X, add.bias=FALSE){
    if (add.bias){
      X <- dplyr::mutate(X, bias =1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
    }
    self$h(as.matrix(X), self$coeff)
  }
), private=list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Vector of length ncol(X)+1 giving upper bounds on the
  #   values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns
  #   of X. If NULL, the values are computed to be the maximum values in the
  #   data, which results in additional privacy loss. Any value in the columns
  #   of X and y larger than the corresponding upper bound is clipped at the
  #   bound.
  # param lower.bounds Vector of length ncol(X)+1 giving lower bounds on the
  #   values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns
  #   of X. If NULL, the values are computed to be the minimum values in the
  #   data, which results in additional privacy loss. Any value in the columns
  #   of X and y smaller than the corresponding lower bound is clipped at the
  #   bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving empirical risk
  #   minimization algorithm.
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
  # description Postprocess coefficients obtained by fitting the model using
  #   differential privacy to ensure they match original inputs X and y.
  #   Effectively undoes the processing done in preprocess_data.
  # param coeff Vector of coefficients obtained by fitting the model.
  # param preprocess List of values returned by preprocess_data and used to
  #   process coeff.
  # return The processed coefficient vector.
  postprocess_coeff = function(coeff, preprocess){
    coeff
  },
  # description Run numerical optimization method to find optimal coefficients
  #   for fitting model given training data and hyperparameters. This function
  #   builds the objective function based on the training data and the provided
  #   search space domain, and runs the constrained optimization algorithm
  #   nloptr from the nloptr package. If h.gr, loss.gr, and regularizer.gr are
  #   all given in the construction of the object, the gradient of the objective
  #   function and the Jacobian of the constraint function are utilized for the
  #   algorithm, and the NLOPT_LD_MMA method is used. If one or more of these
  #   gradient functions are not given, the NLOPT_LN_COBYLA method is used.
  # param X Matrix of data to be fit.
  # param y Vector or matrix of true labels for each row of X.
  # param Delta Nonnegative real number set by the differentially private
  #   algorithm and required for constructing the objective function.
  # param b Numeric perturbation vector randomly drawn by the differentially
  #   private algorithm and required for constructing the objective function.
  # return Vector of fitted coefficients.
  optimize_coeff=function(X, y, Delta, b){
    n <- length(y)
    p <- ncol(X)

    # Get objective function
    objective <- function(par, X, y, Delta, b){
      as.numeric(sum(self$loss(self$h(X,par),y))/n +
                   self$gamma*self$regularizer(par)/n + Delta*par%*%par/(2*n) +
                   t(b)%*%par/n)
    }

    g <- function(par, X, y, Delta, b) self$domain$constraints(par) # For new opt method

    # Get gradient function
    if (!is.null(self$h.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      objective.gr <- function(par, X, y, Delta, b){
        as.numeric(self$h.gr(X,par)%*%self$loss.gr(self$h(X,par),y)/n +
                     self$gamma*self$regularizer.gr(par)/n + b/n + Delta*par/n)
      }
      alg <- "NLOPT_LD_MMA" # For new opt method
      g_jac <- function(par, X, y, Delta, b) self$domain$jacobian(par) # For new opt method
    }
    else {
      objective.gr <- NULL
      alg <- "NLOPT_LN_COBYLA" # For new opt method
      g_jac <- NULL # For new opt method
    }

    # Run optimization
    coeff0 <- numeric(ncol(X))

    # For old opt method
    # opt.res <- optim(coeff0, fn=objective, gr=objective.gr, method="CG",
    #                  X=X, y=y, Delta=Delta, b=b)
    # opt.res$par

    # For new opt method
    opts <- list("algorithm"=alg, xtol_rel = 1e-4)
    opt.res <- nloptr::nloptr(x0=coeff0, eval_f=objective,
                              eval_grad_f=objective.gr, eval_g_ineq=g,
                              eval_jac_g_ineq=g_jac, opts=opts,
                              X=X, y=y, Delta=Delta, b=b)
    opt.res$solution
  }
))

LinearRegressionDP <- R6::R6Class("LinearRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.KST,
  public=list(
    initialize = function(regularizer, eps, delt, gamma, regularizer.gr=NULL){
      domain.linear <- list("constraints"=function(coeff) coeff%*%coeff-length(coeff),
                            "jacobian"=function(coeff) 2*coeff)
      super$initialize(h.linear, loss.squared.error, regularizer, eps, delt,
                       domain.linear, NULL, NULL, gamma,
                       h.gr=h.gr.linear, loss.gr=loss.gr.squared.error,
                       regularizer.gr=regularizer.gr)
    },
    preprocess_data=function(X, y, upper.bounds, lower.bounds, add.bias=TRUE){
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
