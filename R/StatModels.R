#' Sigmoid Hypothesis Function
#'
#' This function implements the sigmoid hypothesis function used for logistic
#' regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the sigmoid function corresponding to each row of
#'   X.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   h.sigmoid(X,coeff)
#'
#' @keywords internal
#'
#' @export
h.sigmoid <- function(X, coeff) e1071::sigmoid(X%*%coeff)

#' Sigmoid Hypothesis Function Gradient
#'
#' This function implements the gradient of the sigmoid hypothesis function with
#' respect to coeff used for logistic regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the gradient of the sigmoid function with respect
#'   to each value of coeff.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   h.gr.sigmoid(X,coeff)
#'
#' @keywords internal
#'
#' @export
h.gr.sigmoid <- function(X, coeff) as.numeric(e1071::dsigmoid(X%*%coeff))*t(X)

#' Linear Hypothesis Function
#'
#' This function implements the linear hypothesis function used for linear SVM
#' and linear regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}} and
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the linear hypothesis function corresponding to
#'   each row of X.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   h.linear(X,coeff)
#'
#' @keywords internal
#'
#' @export
h.linear <- function(X, coeff) X%*%coeff

#' Linear Hypothesis Function Gradient
#'
#' This function implements the gradient of the linear hypothesis function with
#' respect to coeff used for linear SVM and linear regression in the form
#' required by \code{\link{EmpiricalRiskMinimizationDP.CMS}} and
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the gradient of the linear hypothesis function
#'   with respect to each value of coeff.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   h.gr.linear(X,coeff)
#'
#' @keywords internal
#'
#' @export
h.gr.linear <- function(X, coeff) t(X)

#' Cross Entropy Loss Function
#'
#' This function implements cross entropy loss used for logistic regression in
#' the form required by \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param y.hat Vector or matrix of estimated labels.
#' @param y Vector or matrix of true labels.
#' @return Vector or matrix of the cross entropy loss for each element of y.hat
#'   and y.
#' @examples
#'   y.hat <- c(0.1, 0.88, 0.02)
#'   y <- c(0, 1, 0)
#'   loss.cross.entropy(y.hat,y)
#'
#' @keywords internal
#'
#' @export
loss.cross.entropy <- function(y.hat,y) -(y*log(y.hat) + (1-y)*log(1-y.hat))

#' Cross Entropy Loss Function Gradient
#'
#' This function implements cross entropy loss gradient with respect to y.hat
#' used for logistic regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param y.hat Vector or matrix of estimated labels.
#' @param y Vector or matrix of true labels.
#' @return Vector or matrix of the cross entropy loss gradient for each element
#'   of y.hat and y.
#' @examples
#'   y.hat <- c(0.1, 0.88, 0.02)
#'   y <- c(0, 1, 0)
#'   loss.gr.cross.entropy(y.hat,y)
#'
#' @keywords internal
#'
#' @export
loss.gr.cross.entropy <- function(y.hat,y) -y/y.hat + (1-y)/(1-y.hat)

#' Squared Error Loss Function
#'
#' This function implements squared error loss used for linear regression in
#' the form required by \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param y.hat Vector or matrix of estimated labels.
#' @param y Vector or matrix of true labels.
#' @return Vector or matrix of the squared error loss for each element of y.hat
#'   and y.
#' @examples
#'   y.hat <- c(0.1, 0.88, 0.02)
#'   y <- c(0, 1, 0)
#'   loss.squared.error(y.hat,y)
#'
#' @keywords internal
#'
#' @export
loss.squared.error <- function(y.hat,y) (y.hat-y)^2/2

#' Squared error Loss Function Gradient
#'
#' This function implements the squared error loss gradient with respect to
#' y.hat used for linear regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param y.hat Vector or matrix of estimated values.
#' @param y Vector or matrix of true values.
#' @return Vector or matrix of the squared error loss gradient for each element
#'   of y.hat and y.
#' @examples
#'   y.hat <- c(0.1, 0.88, 0.02)
#'   y <- c(-0.1, 1, .2)
#'   loss.gr.squared.error(y.hat,y)
#'
#' @keywords internal
#'
#' @export
loss.gr.squared.error <- function(y.hat,y) y.hat-y

#' Generator for Huber Loss Function
#'
#' This function generates and returns the Huber loss function used for
#' privacy-preserving SVM at the specified value of h in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param h Positive real number for the Huber loss parameter. Lower values more
#'   closely approximate hinge loss. Higher values produce smoother Huber loss
#'   functions.
#' @return Huber loss function with parameter h in the form required by
#'   \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#' @examples
#'   h <- 1
#'   huber <- generate.loss.huber(h)
#'   y.hat <- c(-.5, 1.2, -0.9)
#'   y <- c(-1, 1, -1)
#'   huber(y.hat,y)
#'
#' @keywords internal
#'
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

#' Generator for Huber Loss Function Gradient
#'
#' This function generates and returns the Huber loss function gradient used for
#' privacy-preserving SVM at the specified value of h in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param h Positive real number for the Huber loss parameter. Lower values more
#'   closely approximate hinge loss. Higher values produce smoother Huber loss
#'   functions.
#' @return Huber loss function gradient with parameter h in the form required by
#'   \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#' @examples
#'   h <- 1
#'   huber <- generate.loss.gr.huber(h)
#'   y.hat <- c(-.5, 1.2, -0.9)
#'   y <- c(-1, 1, -1)
#'   huber(y.hat,y)
#'
#' @keywords internal
#'
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

#' l2 Regularizer
#'
#' This function implements the l2 regularizer in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}} and
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Regularizer value at the given coefficient vector.
#' @examples
#'   coeff <- c(0.5,-1,2)
#'   regularizer.l2(coeff)
#'
#' @keywords internal
#'
#' @export
regularizer.l2 <- function(coeff) coeff%*%coeff/2

#' l2 Regularizer Gradient
#'
#' This function implements the l2 regularizer gradient in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}} and
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Regularizer gradient value at the given coefficient vector.
#' @examples
#'   coeff <- c(0.5,-1,2)
#'   regularizer.gr.l2(coeff)
#'
#' @keywords internal
#'
#' @export
regularizer.gr.l2 <- function(coeff) coeff

#' Privacy-preserving Hyperparameter Tuning
#'
#' This function implements the privacy-preserving hyperparameter tuning
#' function \insertCite{chaudhuri2011}{DPpack} using the exponential mechanism.
#' It accepts a list of models with various chosen hyperparameters, a dataset X
#' with corresponding labels y, upper and lower bounds on the columns of X and
#' the values of y, and a boolean indicating whether to add bias in the
#' construction of each of the models. The data are split into m+1 equal groups,
#' where m is the number of models being compared. One group is set aside as the
#' validation group, and each of the other m groups are used to train each of
#' the given m models. The number of errors on the validation set is counted for
#' each model and used as the utility values in the exponential mechanism
#' (\code{\link{ExponentialMechanism}}) to select tuned model in a
#' privacy-preserving way.
#'
#' @param models Vector of model objects, each initialized with a different
#'   combination of hyperparameter values from the search space for tuning.
#'   Currently, only binary classification models are supported. Each model
#'   should be initialized with the same epsilon privacy parameter value eps.
#'   The tuned model satisfies eps-level differential privacy.
#' @param X Dataframe of data to be used in tuning the model. Note it is assumed
#'   the data rows and corresponding labels are randomly shuffled.
#' @param y Vector or matrix of true labels for each row of X.
#' @param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
#'   the values in each column of X and the values in y. The last value in the
#'   vector is assumed to be the upper bound on y, while the first ncol(X)
#'   values are assumed to be in the same order as the corresponding columns of
#'   X. Any value in the columns of X and y larger than the corresponding upper
#'   bound is clipped at the bound.
#' @param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
#'   the values in each column of X and the values in y. The last value in the
#'   vector is assumed to be the lower bound on y, while the first ncol(X)
#'   values are assumed to be in the same order as the corresponding columns of
#'   X. Any value in the columns of X and y smaller than the corresponding lower
#'   bound is clipped at the bound.
#' @param add.bias Boolean indicating whether to add a bias term to X. Defaults
#'   to FALSE.
#' @return Single model object selected from the input list models with tuned
#'   parameters.
#' @examples
#' # Assume X is dataframe meeting assumptions for privacy
#' # Assume 2 columns of X each bounded between -1 and 1
#' # Assume y is 0 or 1 labels for each row of X
#' upper.bounds <- c( 1, 1, 1) # Bounds for X and y
#' lower.bounds <- c(-1,-1, 0) # Bounds for X and y
#' eps <- 1
#'
#' # Grid of possible lambda values for tuning logistric regression model
#' grid.search <- c(100, 1, .0001)
#'
#' lrdp1 <- LogisticRegressionDP$new("l2", eps, grid.search[1])
#' lrdp2 <- LogisticRegressionDP$new("l2", eps, grid.search[2])
#' lrdp3 <- LogisticRegressionDP$new("l2", eps, grid.search[3])
#' models <- c(lrdp1, lrdp2, lrdp3)
#' tuned.model <- tune_model(models, X, y, upper.bounds, lower.bounds)
#' tuned.model$lambda # Gives resulting selected hyperparameter
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
tune_model<- function(models, X, y, upper.bounds=NULL, lower.bounds=NULL,
                      add.bias=FALSE){
  # Make sure values are correct
  if (is.null(upper.bounds)){
    warning(paste("Upper bounds missing and will be calculated from the data.",
                  "This may represent additional privacy loss."));
    upper.bounds <- c(apply(X,2,max),max(y));
  }
  if (is.null(lower.bounds)){
    warning(paste("Lower bounds missing and will be calculated from the data.",
                  "This may represent additional privacy loss."));
    lower.bounds <- c(apply(X,2,min),min(y));
  }
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
  res <- ExponentialMechanism(-z, models[[1]]$eps, 1, candidates=models)
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
#'   The regularizer must be 1-strongly convex and doubly differentiable.
#'   Additionally, it is assumed that if x represents a single row of the
#'   dataset X, then ||x||<=1 for all x. Note that because of this, a bias term
#'   cannot be included without appropriate scaling/preprocessing of the
#'   dataset. To ensure privacy, the add.bias argument in the $fit and $predict
#'   methods should only be utilized in subclasses within this package, not in
#'   this class.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @keywords internal
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
  #' @field D Value only used in child class
  #'   \code{\link{KernelSupportVectorMachineDP}}. Nonnegative integer
  #'   indicating the dimensionality of the transform space approximating the
  #'   kernel. Higher values of D provide better kernel approximations at a cost
  #'   of computational efficiency.
  D = NULL,
  #' @field sampling Value only used in child class
  #'   \code{\link{KernelSupportVectorMachineDP}}. String or sampling function.
  #'   If a string, must be 'Gaussian' (default), indicating to use the sampling
  #'   function corresponding to the Gaussian (radial) kernel approximation. If
  #'   a function, must be of the form sampling(d), where d is the input
  #'   dimension, and return a (d+1)-dimensional vector of samples corresponding
  #'   to the Fourier transform of the kernel to be approximated.
  sampling=NULL,
  #' @field phi Value only used in child class
  #'   \code{\link{KernelSupportVectorMachineDP}}. Function or NULL (default).
  #'   If sampling is given as one of the predefined strings, this input is
  #'   unnecessary. If sampling is a function, this should also be a function of
  #'   the form phi(x, theta), where x is an individual row of of the original
  #'   dataset, and theta is a (d+1)-dimensional vector sampled from the Fourier
  #'   transform of the kernel to be approximated, where d is the dimension of
  #'   x. The function then returns a numeric scalar corresponding to the
  #'   pre-filtered value at the given row with the given sampled vector.
  phi=NULL,
  #' @field gamma Value only used in child class
  #'   \code{\link{KernelSupportVectorMachineDP}}. Positive real number
  #'   corresponding to the Gaussian kernel parameter.
  gamma=NULL,
  #' @field prefilter Value only used in child class
  #'   \code{\link{KernelSupportVectorMachineDP}}. Matrix of pre-filter values
  #'   used in converting data into transform space.
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
  #' @param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the upper bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. Any value in the columns of X and y larger than the corresponding
  #'   upper bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the lower bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. Any value in the columns of X and y smaller than the corresponding
  #'   lower bound is clipped at the bound.
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
  # param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y larger than the corresponding
  #   upper bound is clipped at the bound.
  # param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y smaller than the corresponding
  #   lower bound is clipped at the bound.
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

#' Privacy-preserving Logistic Regression
#'
#' @description This class implements differentially private logistic regression
#'   using the objective perturbation technique
#'   \insertCite{chaudhuri2011}{DPpack}.
#'
#' @details A new model object of this class accepts as inputs a regularizer, an
#'   epsilon value for differential privacy, and a lambda value that scales the
#'   regularizer. The model can then be fit with a dataset X (given as a
#'   data.frame), a set of binary labels y for each row of X, as well as upper
#'   and lower bounds on the possible values for each column of X and for y. In
#'   fitting, the model stores a vector of coefficients coeff which satisfy
#'   epsilon-level differential privacy. These can be released directly, or used
#'   in conjunction with the predict method to predict the label of new
#'   datapoints.
#'
#'   Note that in order to guarantee epsilon-level privacy for the empirical
#'   risk minimization model, certain constraints must be satisfied for the
#'   values used to construct the object, as well as for the data used to fit.
#'   The regularizer must be 1-strongly convex and doubly differentiable. Also,
#'   it is assumed that if x represents a single row of the dataset X, then
#'   \eqn{||x||\le 1} for all \eqn{x}. In order to ensure this constraint is
#'   satisfied, the dataset is preprocessed using provided upper and lower
#'   bounds on the columns of X to scale the values in such a way that this
#'   constraint is met. After the private coefficients are generated, these are
#'   then postprocessed and un-scaled so that the stored coefficients correspond
#'   to the original data. This does not result in additional privacy loss as
#'   long as the upper and lower bounds provided when fitting the model do not
#'   depend directly on the data. Due to this constraint on \eqn{x}, it is best
#'   to avoid using a bias term in the model whenever possible. If a bias term
#'   must be used, the issue can be partially circumvented by adding a constant
#'   column to X before fitting the model, which will be scaled along with the
#'   rest of X. The \code{fit} method contains functionality to add a column of
#'   constant 1s to X before scaling, if desired.
#'
#'   The preprocessing of X is done as follows. First, the largest in absolute
#'   value of the upper and lower bounds on each column are used to scale each
#'   column individually such that the largest value in each column is at most 1
#'   in absolute value. Second, each value in X is divided by the square root of
#'   the number of predictors of X (including bias term). These two scalings
#'   ensure that each row of X satisfies the necessary constraints for
#'   differential privacy. Additionally, the labels y are assumed to be either 0
#'   or 1. If different values are provided, they are coerced to be either 0 or
#'   1 prior to fitting the model. Values in y that are \eqn{\le} 0 are assigned
#'   to be 0, while values in y \eqn{>} 0 are assigned to be 1. Accordingly, new
#'   predicted labels are output as either 0 or 1.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
LogisticRegressionDP <- R6::R6Class("LogisticRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.CMS,
  public=list(
  #' @description Create a new LogisticRegressionDP object.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   as given in regularizer field description. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param lambda Nonnegative real number representing the regularization
  #'   constant.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularizer function function with respect to coeff. Must have form as
  #'   given in regularizer.gr field description for parent class. If not given,
  #'   gradients are not used to compute the coefficient values in fitting the
  #'   model.
  #'
  #' @examples
  #' # Construct object for logistic regression
  #' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
  #' eps <- 1
  #' lambda <- 0.1
  #' regularizer.gr <- function(coeff) coeff # If function given for regularizer
  #' lrdp <- LogisticRegressioinDP$new(regularizer, eps, lambda,
  #'                                   regularizer.gr)
  #'
  #' @return A new LogisticRegressionDP object.
  initialize = function(regularizer, eps, lambda, regularizer.gr = NULL){
    super$initialize(h.sigmoid, loss.cross.entropy, regularizer, eps,
                    lambda, 1/4, h.gr.sigmoid, loss.gr.cross.entropy,
                    regularizer.gr)
  },
  #' @description Predict label(s) for given X using the fitted coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as X used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #' @param raw.value Boolean indicating whether to return the raw predicted
  #'   value or the rounded class label. If FALSE (default), rounds the values
  #'   to 0 or 1. If TRUE, returns the raw score from the logistic regression.
  #'
  #' @examples
  #' # Assume Xtest is a new dataframe of the same form as X from fit
  #' # method example, with true labels ytest
  #' # Also assume lrdp$fit() has already been run on training data
  #' predicted.y <- lrdp$predict(Xtest)
  #' n.errors <- sum(predicted.y!=ytest)
  #'
  #' @return Matrix of predicted labels corresponding to each row of X.
  predict = function(X, add.bias=FALSE, raw.value=FALSE){
      if (raw.value) super$predict(X, add.bias)
      else round(super$predict(X, add.bias))
    }
), private=list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y larger than the corresponding
  #   upper bound is clipped at the bound.
  # param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y smaller than the corresponding
  #   lower bound is clipped at the bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving logistic regression
  #   algorithm.
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
    res <- super$preprocess_data(X,y,upper.bounds,lower.bounds, add.bias)
    X <- res$X
    y <- res$y
    upper.bounds <- res$upper.bounds
    lower.bounds <- res$lower.bounds

    # Process X
    p <- length(lower.bounds)-1
    lb <- lower.bounds[1:p]
    ub <- upper.bounds[1:p]
    tmp <- matrix(c(lb,ub),ncol=2)
    scale1 <- apply(abs(tmp),1,max)
    scale2 <- sqrt(p)
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
  # description Postprocess coefficients obtained by fitting the model using
  #   differential privacy to ensure they match original inputs X and y.
  #   Effectively undoes the processing done in preprocess_data.
  # param coeff Vector of coefficients obtained by fitting the model.
  # param preprocess List of values returned by preprocess_data and used to
  #   process coeff.
  # return The processed coefficient vector.
  postprocess_coeff = function(coeff, preprocess){
    coeff/(preprocess$scale1*preprocess$scale2)
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
  optimize_coeff = function(X, y, Delta, b){
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

#' Generator for Sampling Distribution Function for Gaussian Kernel
#'
#' This function generates and returns a sampling function corresponding to the
#' Fourier transform of a Gaussian kernel with parameter gamma
#' \insertCite{chaudhuri2011}{DPpack} of form needed for
#' \code{\link{KernelSupportVectorMachineDP}}.
#'
#' @param gamma Positive real number for the Gaussian (radial) kernel parameter.
#' @return Sampling function for the Gaussian kernel of form required by
#'   \code{\link{KernelSupportVectorMachineDP}}.
#' @examples
#'   gamma <- 1
#'   sample <- generate.sampling(gamma)
#'   d <- 5
#'   sample(d)
#'
#' @keywords internal
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
generate.sampling <- function(gamma){
  function(d){
    omega <- rnorm(d,sd=sqrt(2*gamma))
    phi <- runif(1,-pi,pi)
    c(omega,phi)
  }
}

#' Transform Function for Gaussian Kernel Approximation
#'
#' This function maps an input data row x with a given prefilter to an output
#' value in such a way as to approximate the Gaussian kernel
#' \insertCite{chaudhuri2011}{DPpack}.
#'
#' @param x Vector or matrix corresponding to one row of the dataset X.
#' @param theta Randomly sampled prefilter vector of length n+1, where n is the
#'   length of x.
#' @return Mapped value corresponding to one element of the transformed space.
#' @examples
#'   x <- c(1,2,3)
#'   theta <- c(0.1, 1.1, -0.8, 3)
#'   phi.gaussian(x, theta)
#'
#' @keywords internal
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
phi.gaussian <- function(x, theta){
  d <- length(x)
  cos(x%*%theta[1:d] + theta[d+1])
}

#' Privacy-preserving Linear Support Vector Machine
#'
#' @description This class implements a differentially private support vector
#'   machine (SVM) using the objective perturbation technique
#'   \insertCite{chaudhuri2011}{DPpack}.
#'
#' @details A new model object of this class accepts as inputs a regularizer, an
#'   epsilon value for differential privacy, and a lambda value that scales the
#'   regularizer. The model can then be fit with a dataset X (given as a
#'   data.frame), a set of binary labels y for each row of X, as well as upper
#'   and lower bounds on the possible values for each column of X and for y. In
#'   fitting, the model stores a vector of coefficients coeff which satisfy
#'   epsilon-level differential privacy. These can be released directly, or used
#'   in conjunction with the predict method to predict the label of new
#'   datapoints.
#'
#'   Note that in order to guarantee epsilon-level privacy for the empirical
#'   risk minimization model, certain constraints must be satisfied for the
#'   values used to construct the object, as well as for the data used to fit.
#'   First, the loss function is assumed to be doubly differentiable. The hinge
#'   loss, which is typically used for linear support vector machines, is not
#'   differentiable at 1. Thus, to satisfy this constraint, this class utilizes
#'   the Huber loss, a smooth approximation to the hinge loss. The level of
#'   approximation to the hinge loss is determined by a user-specified constant,
#'   h, which defaults to 1. Additionally, the regularizer must be 1-strongly
#'   convex and doubly differentiable. Finally, it is assumed that if x
#'   represents a single row of the dataset X, then \eqn{||x||\le 1} for all
#'   \eqn{x}. In order to ensure this constraint is satisfied, the dataset is
#'   preprocessed using provided upper and lower bounds on the columns of X to
#'   scale the values in such a way that this constraint is met. After the
#'   private coefficients are generated, these are then postprocessed and
#'   un-scaled so that the stored coefficients correspond to the original data.
#'   This does not result in additional privacy loss as long as the upper and
#'   lower bounds provided when fitting the model do not depend directly on the
#'   data. Due to this constraint on \eqn{x}, it is best to avoid using a bias
#'   term in the model whenever possible. If a bias term must be used, the issue
#'   can be partially circumvented by adding a constant column to X before
#'   fitting the model, which will be scaled along with the rest of X. The
#'   \code{fit} method contains functionality to add a column of constant 1s to
#'   X before scaling, if desired.
#'
#'   The preprocessing of X is done as follows. First, the largest in absolute
#'   value of the upper and lower bounds on each column are used to scale each
#'   column individually such that the largest value in each column is at most 1
#'   in absolute value. Second, each value in X is divided by the square root of
#'   the number of predictors of X (including bias term). These two scalings
#'   ensure that each row of X satisfies the necessary constraints for
#'   differential privacy. Additionally, the labels y are assumed to be either
#'   -1 or 1. If different values are provided, they are coerced to be either -1
#'   or 1 prior to fitting the model. Values in y that are \eqn{\le} 0 are
#'   assigned to be -1, while values in y \eqn{>} 0 are assigned to be 1.
#'   Accordingly, new predicted labels are output as either -1 or 1.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
SupportVectorMachineDP <- R6::R6Class("SupportVectorMachineDP",
  inherit=EmpiricalRiskMinimizationDP.CMS,
  public=list(
  #' @description Create a new SupportVectorMachineDP object.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   as given in regularizer field description. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param lambda Nonnegative real number representing the regularization
  #'   constant.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularizer function function with respect to coeff. Must have form as
  #'   given in regularizer.gr field description for parent class. If not given,
  #'   gradients are not used to compute the coefficient values in fitting the
  #'   model.
  #' @param huber.h Positive real number indicating the degree to which the
  #'   Huber loss approximates the hinge loss. A smaller value indicates closer
  #'   resemblance to the hinge loss, but a larger value of the constant c used
  #'   in the objective perturbation algorithm, meaning more noise needs to be
  #'   added to ensure differential privacy. Conversely, larger values for this
  #'   parameter represent looser approximations to the hinge loss, but smaller
  #'   c values and less noise to ensure privacy. Defaults to 1.
  #'
  #' @examples
  #' # Construct object for SVM
  #' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
  #' eps <- 1
  #' lambda <- 0.1
  #' regularizer.gr <- function(coeff) coeff # If function given for regularizer
  #' huber.h <- 1
  #' svmdp <- SupportVectorMachineDP$new(regularizer, eps, lambda,
  #'                                   regularizer.gr, huber.h)
  #'
  #' @return A new SupportVectorMachineDP object.
  initialize = function(regularizer, eps, lambda, regularizer.gr=NULL, huber.h=1){
    super$initialize(h.linear, generate.loss.huber(huber.h), regularizer, eps,
                     lambda, 1/(2*huber.h), h.gr.linear,
                     generate.loss.gr.huber(huber.h), regularizer.gr)
  },
  #' @description Predict label(s) for given X using the fitted coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as X used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #' @param raw.value Boolean indicating whether to return the raw predicted
  #'   value or the rounded class label. If FALSE (default), rounds the values
  #'   to -1 or 1. If TRUE, returns the raw score from the SVM model.
  #'
  #' @examples
  #' # Assume Xtest is a new dataframe of the same form as X from fit
  #' # method example, with true labels ytest
  #' # Also assume svmdp$fit() has already been run on training data
  #' predicted.y <- svmdp$predict(Xtest)
  #' n.errors <- sum(predicted.y!=ytest)
  #'
  #' @return Matrix of predicted labels corresponding to each row of X.
  predict = function(X, add.bias=FALSE, raw.value=FALSE){
      if (raw.value) super$predict(X, add.bias)
      else sign(super$predict(X, add.bias))
    }
), private=list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y larger than the corresponding
  #   upper bound is clipped at the bound.
  # param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y smaller than the corresponding
  #   lower bound is clipped at the bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving SVM algorithm.
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
    res <- super$preprocess_data(X,y,upper.bounds,lower.bounds, add.bias)
    X <- res$X
    y <- res$y
    upper.bounds <- res$upper.bounds
    lower.bounds <- res$lower.bounds

    # Process X
    p <- length(lower.bounds)-1
    lb <- lower.bounds[1:p]
    ub <- upper.bounds[1:p]
    tmp <- matrix(c(lb,ub),ncol=2)
    scale1 <- apply(abs(tmp),1,max)
    scale2 <- sqrt(p)
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
  # description Postprocess coefficients obtained by fitting the model using
  #   differential privacy to ensure they match original inputs X and y.
  #   Effectively undoes the processing done in preprocess_data.
  # param coeff Vector of coefficients obtained by fitting the model.
  # param preprocess List of values returned by preprocess_data and used to
  #   process coeff.
  # return The processed coefficient vector.
  postprocess_coeff = function(coeff, preprocess){
    coeff/(preprocess$scale1*preprocess$scale2)
  }
))

#' Privacy-preserving Nonlinear Kernel Support Vector Machine
#'
#' @description This class implements a differentially private support vector
#'   machine with a nonlinear kernel using the objective perturbation technique
#'   \insertCite{chaudhuri2011}{DPpack}.
#'
#' @details A new model object of this class accepts as inputs various functions
#'   and hyperparameters related to the specific regression problem. The model
#'   can then be fit with a dataset X (given as a data.frame), a set of binary
#'   labels y for each row of X, as well as upper and lower bounds on the
#'   possible values for each column of X and for y. The algorithm first
#'   transforms the data using a kernel approximation method, then runs the
#'   privacy-preserving linear SVM model on the results
#'   (\code{\link{SupportVectorMachineDP}}). In fitting the model, the class
#'   stores a method \code{XtoV} that transforms additional input data X into
#'   the new dimension, as well as a vector of coefficients coeff for the
#'   transformed data V which satisfy epsilon-level differential privacy. These
#'   can be released directly, or used in conjunction with the predict method to
#'   predict the label of new datapoints.
#'
#'   Note that in order to guarantee epsilon-level privacy for the empirical
#'   risk minimization model, certain constraints must be satisfied for the
#'   values used to construct the object, as well as for the data used to fit.
#'   First, the loss function is assumed to be doubly differentiable. The hinge
#'   loss, which is typically used for linear support vector machines, is not
#'   differentiable at 1. Thus, to satisfy this constraint, this class utilizes
#'   the Huber loss, a smooth approximation to the hinge loss. The level of
#'   approximation to the hinge loss is determined by a user-specified constant,
#'   h, which defaults to 1. Additionally, the regularizer must be 1-strongly
#'   convex and doubly differentiable.
#'
#'   The labels y are assumed to be either -1 or 1. If different values are
#'   provided, they are coerced to be either -1 or 1 prior to fitting the model.
#'   Values in y that are \eqn{\le} 0 are assigned to be -1, while values in y
#'   \eqn{>} 0 are assigned to be 1. Accordingly, new predicted labels are
#'   output as either -1 or 1.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
KernelSupportVectorMachineDP <- R6::R6Class("KernelSupportVectorMachineDP",
  inherit=SupportVectorMachineDP,
  public=list(
  #' @description Create a new KernelSupportVectorMachineDP object.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   as given in regularizer field description. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param lambda Nonnegative real number representing the regularization
  #'   constant.
  #' @param D Nonnegative integer indicating the dimensionality of the transform
  #'   space approximating the kernel. Higher values of D provide better kernel
  #'   approximations at a cost of computational efficiency.
  #' @param sampling String or sampling function. If a string, must be
  #'   'Gaussian' (default), indicating to use the sampling function
  #'   corresponding to the Gaussian (radial) kernel approximation. If a
  #'   function, must take as input a dimension d and return a (d+1)-dimensional
  #'   vector of samples corresponding to the Fourier transform of the kernel to
  #'   be approximated.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularizer function function with respect to coeff. Must have form as
  #'   given in regularizer.gr field description for parent class. If not given,
  #'   gradients are not used to compute the coefficient values in fitting the
  #'   model.
  #' @param huber.h Positive real number indicating the degree to which the
  #'   Huber loss approximates the hinge loss. A smaller value indicates closer
  #'   resemblance to the hinge loss, but a larger value of the constant c used
  #'   in the objective perturbation algorithm, meaning more noise needs to be
  #'   added to ensure differential privacy. Conversely, larger values for this
  #'   parameter represent looser approximations to the hinge loss, but smaller
  #'   c values and less noise to ensure privacy. Defaults to 1.
  #' @param phi Function or NULL (default). If sampling is given as one of the
  #'   predefined strings, this input is unnecessary. If sampling is a function,
  #'   this should also be a function that takes as inputs an individual row of
  #'   of the original dataset x, and a (d+1)-dimensional vector theta sampled
  #'   from the Fourier transform of the kernel to be approximated, where d is
  #'   the dimension of x. The function then outputs a numeric scalar
  #'   corresponding to the pre-filtered value at the given row with the given
  #'   sampled vector.
  #' @param gamma Positive real number corresponding to the Gaussian kernel
  #'   parameter.
  #'
  #' @examples
  #' # Construct object for nonlinear SVM
  #' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
  #' eps <- 1
  #' lambda <- 0.1
  #' D <- 20
  #' sampling <- "Gaussian"
  #' regularizer.gr <- function(coeff) coeff # If function given for regularizer
  #' huber.h <- 1
  #' ksvmdp <- KernelSupportVectorMachineDP$new(regularizer, eps, lambda, D,
  #'                                   sampling, regularizer.gr, huber.h)
  #'
  #' @return A new KernelSupportVectorMachineDP object.
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
  #' @description Convert input data X into transformed data V. Uses sampled
  #'   pre-filter values and the provided mapping function to produce
  #'   D-dimensional data V on which to train the model or predict future
  #'   values.
  #' @param X Matrix corresponding to the original dataset.
  #' @return Matrix V of size n by D representing the transformed dataset, where
  #'   n is the number of rows of X, and D is the provided transformed space
  #'   dimension.
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
  #' @description Predict label(s) for given X using the XtoV method and fitted
  #'   coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as X used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to X.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #' @param raw.value Boolean indicating whether to return the raw predicted
  #'   value or the rounded class label. If FALSE (default), rounds the values
  #'   to -1 or 1. If TRUE, returns the raw score from the SVM model.
  #'
  #' @examples
  #' # Assume Xtest is a new dataframe of the same form as X from fit
  #' # method example, with true labels ytest
  #' # Also assume ksvmdp$fit() has already been run on training data
  #' predicted.y <- ksvmdp$predict(Xtest)
  #' n.errors <- sum(predicted.y!=ytest)
  #'
  #' @return Matrix of predicted labels corresponding to each row of X.
  predict = function(X, add.bias=FALSE, raw.value=FALSE){
      if (add.bias){
        X <- dplyr::mutate(X, bias =1)
        X <- X[, c(ncol(X), 1:(ncol(X)-1))]
      }
      V <- self$XtoV(as.matrix(X))
      super$predict(V, add.bias=FALSE, raw.value)
    }
), private=list(
  # description Preprocess input data and bounds to convert to kernel
  #   approximation dimension. Converts original dataset X to transformed datset
  #   V, then preprocesses V and y. Additionally, ensures V and y meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y larger than the corresponding
  #   upper bound is clipped at the bound.
  # param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y smaller than the corresponding
  #   lower bound is clipped at the bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for V, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving nonlinear SVM algorithm.
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
#'   can then be fit with a dataset X (given as a data.frame), a set of true
#'   values y for each row of X, as well as upper and lower bounds on the
#'   possible values for each column of X and for y. In fitting, the model
#'   stores a vector of coefficients coeff which satisfy epsilon-level
#'   differential privacy. These can be released directly, or used in
#'   conjunction with the predict method to predict the value of new datapoints.
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
#' @keywords internal
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
  #' @field delta Nonnegative real number defining the delta parameter for
  #'   approximate differential privacy. If set to 0, pure differential
  #'   privacy is used.
  delta = NULL,
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
  #'   to ensure differential privacy, the function must be convex.
  #' @param eps Positive real number defining the epsilon privacy budget. If
  #'   set to Inf, runs algorithm without differential privacy.
  #' @param delta Nonnegative real number defining the delta parameter for
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
  #' delta <- 1
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
  #' ermdp <- EmpiricalRiskMinimizationDP.KST$new(h, loss, 'l2', eps, delta,
  #'                                              domain, zeta, lambda,
  #'                                              gamma, h.gr, loss.gr,
  #'                                              regularizer.gr)
  #'
  #' @return A new EmpiricalRiskMinimizationDP.KST object.
  initialize = function(h, loss, regularizer, eps, delta, domain, zeta, lambda,
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
    self$delta <- delta
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
  #' @param y Vector or matrix of true values for each row of X.
  #' @param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the upper bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. Any value in the columns of X and y larger than the corresponding
  #'   upper bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds
  #'   on the values in each column of X and the values in y. The last value in
  #'   the vector is assumed to be the lower bound on y, while the first ncol(X)
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of X. Any value in the columns of X and y smaller than the corresponding
  #'   lower bound is clipped at the bound.
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
    if (is.null(self$delta) || self$delta==0){
      norm.b <- rgamma(1, p, rate=self$eps/(2*self$zeta))
      direction.b <- rnorm(p)
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b*direction.b
    } else{
      mu <- numeric(p)
      Sigma <- ((self$zeta^2*(8*log(2/self$delta)+4*self$eps))/
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
  # param y Vector or matrix of true values for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y larger than the corresponding
  #   upper bound is clipped at the bound.
  # param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y smaller than the corresponding
  #   lower bound is clipped at the bound.
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

#' Privacy-preserving Linear Regression
#'
#' @description This class implements differentially private linear regression
#'   using the objective perturbation technique \insertCite{Kifer2012}{DPpack}.
#'
#' @details A new model object of this class accepts as inputs functions and
#'   hyperparameters indicating the privacy levels and regularization. The model
#'   can then be fit with a dataset X (given as a data.frame), a set of true
#'   values y for each row of X, as well as upper and lower bounds on the
#'   possible values for each column of X and for y. In fitting, the model
#'   stores a vector of coefficients coeff which satisfy epsilon-level
#'   differential privacy. These can be released directly, or used in
#'   conjunction with the predict method to predict the label of new datapoints.
#'
#'   Note that in order to guarantee epsilon-level privacy for the empirical
#'   risk minimization model, certain constraints must be satisfied for the
#'   values used to construct the object, as well as for the data used to fit.
#'   The regularization function must be convex. Also, it is assumed that if x
#'   represents a single row of the dataset X, then \eqn{||x||\le p} for all
#'   \eqn{x}, where \eqn{p} is the number of predictors (including any possible
#'   intercept term). In order to ensure this constraint is satisfied, the
#'   dataset is preprocessed using provided upper and lower bounds on the
#'   columns of X to scale the values in such a way that this constraint is met.
#'   After the private coefficients are generated, these are then postprocessed
#'   and un-scaled so that the stored coefficients correspond to the original
#'   data. This does not result in additional privacy loss as long as the upper
#'   and lower bounds provided when fitting the model do not depend directly on
#'   the data. Due to this constraint on \eqn{x}, it is best to avoid using an
#'   intercept term in the model whenever possible. If an intercept term must be
#'   used, the issue can be partially circumvented by adding a constant column
#'   to X before fitting the model, which will be scaled along with the rest of
#'   X. The \code{fit} method contains functionality to add a column of constant
#'   1s to X before scaling, if desired.
#'
#'   The preprocessing of X is done as follows. First, the largest in absolute
#'   value of the upper and lower bounds on each column are used to scale each
#'   column individually such that the largest value in each column is at most 1
#'   in absolute value. Second, each value in X is divided by the square root of
#'   the number of predictors of X (including bias term). These two scalings
#'   ensure that each row of X satisfies the necessary constraints for
#'   differential privacy. Additionally, the true values y are assumed to be
#'   bounded between \eqn{-p} and \eqn{p}. To accomodate this assumption, the
#'   provided upper and lower bounds on y are used to shift the true values to
#'   be centered around 0, then the bounds are used to scale the true values to
#'   be within this range.
#'
#' @references \insertRef{Kifer2012}{DPpack}
#'
#' @export
LinearRegressionDP <- R6::R6Class("LinearRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.KST,
  public=list(
  #' @description Create a new LinearRegressionDP object.
  #' @param regularizer String or regularization function. If a string, must
  #'   be 'l2', indicating to use l2 regularization. If a function, must have
  #'   form as given in regularizer field description. Additionally, in order
  #'   to ensure differential privacy, the function must be convex.
  #' @param eps Positive real number defining the epsilon privacy budget. If
  #'   set to Inf, runs algorithm without differential privacy.
  #' @param delta Nonnegative real number defining the delta parameter for
  #'   approximate differential privacy. If set to 0, pure differential
  #'   privacy is used.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularizer function function with respect to coeff. Must have form as
  #'   given in regularizer.gr field description. If not given, gradients are
  #'   not used to compute the coefficient values in fitting the model.
  #'
  #' @examples
  #' # Construct object for linear regression
  #' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
  #' eps <- 1
  #' delta <- 1
  #' gamma <- 1
  #' regularizer.gr <- function(coeff) coeff
  #'
  #' lrdp <- LinearRegressionDP$new('l2', eps, delta, gamma,
  #'                                              regularizer.gr)
  #'
  #' @return A new LinearRegressionDP object.
  initialize = function(regularizer, eps, delta, gamma, regularizer.gr=NULL){
      domain.linear <- list("constraints"=function(coeff) coeff%*%coeff-length(coeff),
                            "jacobian"=function(coeff) 2*coeff)
      super$initialize(h.linear, loss.squared.error, regularizer, eps, delta,
                       domain.linear, NULL, NULL, gamma,
                       h.gr=h.gr.linear, loss.gr=loss.gr.squared.error,
                       regularizer.gr=regularizer.gr)
    }
), private=list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true values for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X)+1 giving upper bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the upper bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y larger than the corresponding
  #   upper bound is clipped at the bound.
  # param lower.bounds Numeric vector of length ncol(X)+1 giving lower bounds on
  #   the values in each column of X and the values in y. The last value in the
  #   vector is assumed to be the lower bound on y, while the first ncol(X)
  #   values are assumed to be in the same order as the corresponding columns of
  #   X. Any value in the columns of X and y smaller than the corresponding
  #   lower bound is clipped at the bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving empirical risk
  #   minimization algorithm.
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
    scale.X <- apply(abs(tmp),1,max)
    # Need each row at most norm sqrt(p), already achieved by first scaling
    X.norm <- t(t(X)/scale.X)

    # Process y
    lb.y <- lower.bounds[length(lower.bounds)]
    ub.y <- upper.bounds[length(upper.bounds)]
    shift.y <- lb.y + (ub.y - lb.y)/2 # Subtracting this centers at 0
    scale.y <- (ub.y-lb.y)/(2*p) # Dividing by this scales to max size p in abs()
    y.norm <- (y-shift.y)/scale.y

    list(X=X.norm,y=y.norm, scale.X=scale.X, shift.y=shift.y, scale.y=scale.y)
  },
  # description Postprocess coefficients obtained by fitting the model using
  #   differential privacy to ensure they match original inputs X and y.
  #   Effectively undoes the processing done in preprocess_data.
  # param coeff Vector of coefficients obtained by fitting the model.
  # param preprocess List of values returned by preprocess_data and used to
  #   process coeff.
  # return The processed coefficient vector.
  postprocess_coeff=function(coeff, preprocess){
    # Undo X.norm
    coeff <- coeff/preprocess$scale.X

    # Undo y.norm
    coeff <- coeff*preprocess$scale.y
    coeff[1] <- coeff[1] + preprocess$shift.y # Assumes coeff[1] is bias term
    coeff
  }
))
