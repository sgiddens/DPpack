#' Sigmoid Map Function
#'
#' This function implements the sigmoid map function from data X to labels y
#' used for logistic regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the sigmoid function corresponding to each row of
#'   X.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   mapXy.sigmoid(X,coeff)
#'
#' @keywords internal
#'
#' @export
mapXy.sigmoid <- function(X, coeff) e1071::sigmoid(X%*%coeff)

#' Sigmoid Map Function Gradient
#'
#' This function implements the gradient of the sigmoid map function with
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
#'   mapXy.gr.sigmoid(X,coeff)
#'
#' @keywords internal
#'
#' @export
mapXy.gr.sigmoid <- function(X, coeff) as.numeric(e1071::dsigmoid(X%*%coeff))*t(X)

#' Linear Map Function
#'
#' This function implements the linear map function from data X to labels y used
#' for linear SVM and linear regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}} and
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the linear map function corresponding to each row
#'   of X.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   mapXy.linear(X,coeff)
#'
#' @keywords internal
#'
#' @export
mapXy.linear <- function(X, coeff) X%*%coeff

#' Linear Map Function Gradient
#'
#' This function implements the gradient of the linear map function with respect
#' to coeff used for linear SVM and linear regression in the form required by
#' \code{\link{EmpiricalRiskMinimizationDP.CMS}} and
#' \code{\link{EmpiricalRiskMinimizationDP.KST}}.
#'
#' @param X Matrix of data.
#' @param coeff Vector or matrix of coefficients or weights.
#' @return Matrix of values of the gradient of the linear map function with
#'   respect to each value of coeff.
#' @examples
#'   X <- matrix(c(1,2,3,4,5,6),nrow=2)
#'   coeff <- c(0.5,-1,2)
#'   mapXy.gr.linear(X,coeff)
#'
#' @keywords internal
#'
#' @export
mapXy.gr.linear <- function(X, coeff) t(X)

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
#'   h <- 0.5
#'   huber <- generate.loss.huber(h)
#'   y.hat <- c(-.5, 1.2, -0.9)
#'   y <- c(-1, 1, -1)
#'   huber(y.hat,y)
#'   huber(y.hat, y, w=c(0.1, 0.5, 1)) # Weights observation-level loss
#'
#' @keywords internal
#'
#' @export
generate.loss.huber <- function(h){
  # Weight vector w included for WeightedERMDP.CMS
  function(y.hat, y, w = NULL){
    if (is.null(w)) w <- rep(1, length(y))
    if (length(w)!=length(y)) stop("Weight vector w must be same length as label vector y.")
    z <- y.hat*y
    mask1 <- z>(1+h)
    mask2 <- abs(1-z)<=h
    mask3 <- z<(1-h)
    z[mask1] <- 0
    z[mask2] <- w[mask2]*(1+h-z[mask2])^2/(4*h)
    z[mask3] <- w[mask3]*(1-z[mask3])
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
#'   huber(y.hat, y, w=c(0.1, 0.5, 1)) # Weights observation-level loss gradient
#'
#' @keywords internal
#'
#' @export
generate.loss.gr.huber <- function(h){
  function(y.hat, y, w = NULL){
    # Weight vector w included for WeightedERMDP.CMS
    if (is.null(w)) w <- rep(1, length(y))
    if (length(w)!=length(y)) stop("Weight vector w must be same length as label vector y.")
    z <- y.hat*y
    mask1 <- z>(1+h)
    mask2 <- abs(1-z)<=h
    mask3 <- z<(1-h)
    z[mask1] <- 0
    z[mask2] <- -w[mask2]*y[mask2]*(1+h-z[mask2])/(2*h)
    z[mask3] <- -w[mask3]*y[mask3]
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

#' Privacy-preserving Hyperparameter Tuning for Binary Classification Models
#'
#' This function implements the privacy-preserving hyperparameter tuning
#' function for binary classification \insertCite{chaudhuri2011}{DPpack} using
#' the exponential mechanism. It accepts a list of models with various chosen
#' hyperparameters, a dataset X with corresponding labels y, upper and lower
#' bounds on the columns of X, and a boolean indicating whether to add bias in
#' the construction of each of the models. The data are split into m+1 equal
#' groups, where m is the number of models being compared. One group is set
#' aside as the validation group, and each of the other m groups are used to
#' train each of the given m models. The number of errors on the validation set
#' is counted for each model and used as the utility values in the exponential
#' mechanism (\code{\link{ExponentialMechanism}}) to select a tuned model in a
#' privacy-preserving way.
#'
#' @param models Vector of binary classification model objects, each initialized
#'   with a different combination of hyperparameter values from the search space
#'   for tuning. Each model should be initialized with the same epsilon privacy
#'   parameter value eps. The tuned model satisfies eps-level differential
#'   privacy.
#' @param X Dataframe of data to be used in tuning the model. Note it is assumed
#'   the data rows and corresponding labels are randomly shuffled.
#' @param y Vector or matrix of true labels for each row of X.
#' @param upper.bounds Numeric vector giving upper bounds on the values in each
#'   column of X. Should be of length ncol(X). The values are assumed to be in
#'   the same order as the corresponding columns of X. Any value in the columns
#'   of X larger than the corresponding upper bound is clipped at the bound.
#' @param lower.bounds Numeric vector giving lower bounds on the values in each
#'   column of X. Should be of length ncol(X). The values are assumed to be in
#'   the same order as the corresponding columns of X. Any value in the columns
#'   of X smaller than the corresponding lower bound is clipped at the bound.
#' @param add.bias Boolean indicating whether to add a bias term to X. Defaults
#'   to FALSE.
#' @param weights Numeric vector of observation weights of the same length as
#'   \code{y}.
#' @param weights.upper.bound Numeric value representing the global or public
#'   upper bound on the weights.
#' @return Single model object selected from the input list models with tuned
#'   parameters.
#' @examples
#' # Build train dataset X and y, and test dataset Xtest and ytest
#' N <- 200
#' K <- 2
#' X <- data.frame()
#' y <- data.frame()
#' for (j in (1:K)){
#'   t <- seq(-.25,.25,length.out = N)
#'   if (j==1) m <- stats::rnorm(N,-.2,.1)
#'   if (j==2) m <- stats::rnorm(N, .2,.1)
#'   Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
#'   ytemp <- data.frame(matrix(j-1, N, 1))
#'   X <- rbind(X, Xtemp)
#'   y <- rbind(y, ytemp)
#' }
#' Xtest <- X[seq(1,(N*K),10),]
#' ytest <- y[seq(1,(N*K),10),,drop=FALSE]
#' X <- X[-seq(1,(N*K),10),]
#' y <- y[-seq(1,(N*K),10),,drop=FALSE]
#' y <- as.matrix(y)
#' weights <- rep(1, nrow(y)) # Uniform weighting
#' weights[nrow(y)] <- 0.5 # half weight for last observation
#' wub <- 1 # Public upper bound for weights
#'
#' # Grid of possible gamma values for tuning logistic regression model
#' grid.search <- c(100, 1, .0001)
#'
#' # Construct objects for SVM parameter tuning
#' eps <- 1 # Privacy budget should be the same for all models
#' svmdp1 <- svmDP$new("l2", eps, grid.search[1], perturbation.method='output')
#' svmdp2 <- svmDP$new("l2", eps, grid.search[2], perturbation.method='output')
#' svmdp3 <- svmDP$new("l2", eps, grid.search[3], perturbation.method='output')
#' models <- c(svmdp1, svmdp2, svmdp3)
#'
#' # Tune using data and bounds for X based on its construction
#' upper.bounds <- c( 1, 1)
#' lower.bounds <- c(-1,-1)
#' tuned.model <- tune_classification_model(models, X, y, upper.bounds,
#'                                          lower.bounds, weights=weights,
#'                                          weights.upper.bound=wub)
#' tuned.model$gamma # Gives resulting selected hyperparameter
#'
#' # tuned.model result can be used the same as a trained LogisticRegressionDP model
#' # Predict new data points
#' predicted.y <- tuned.model$predict(Xtest)
#' n.errors <- sum(predicted.y!=ytest)
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
tune_classification_model<- function(models, X, y, upper.bounds, lower.bounds,
                      add.bias=FALSE, weights=NULL, weights.upper.bound=NULL){
  # Make sure values are correct
  m <- length(models)
  n <- length(y)
  # Split data into m+1 groups
  validateX <- X[seq(m+1,n,m+1),,drop=FALSE]
  validatey <- y[seq(m+1,n,m+1)]
  z <- numeric(m)
  for (i in 1:m){
    subX <- X[seq(i,n,m+1),,drop=FALSE]
    suby <- y[seq(i,n,m+1)]
    if (!is.null(weights)) subWeights <- weights[seq(i,n,m+1)]
    if (is.null(weights)){
      models[[i]]$fit(subX,suby,upper.bounds,lower.bounds,add.bias)
    } else{
      models[[i]]$fit(subX,suby,upper.bounds,lower.bounds,add.bias,
                      weights=subWeights, weights.upper.bound=weights.upper.bound)
    }
    validatey.hat <- models[[i]]$predict(validateX,add.bias,raw.value=FALSE)
    z[i] <- sum(validatey!=validatey.hat)
  }
  res <- ExponentialMechanism(-z, models[[1]]$eps, 1, candidates=models)
  res[[1]]
}

#' Privacy-preserving Empirical Risk Minimization for Binary Classification
#'
#' @description This class implements differentially private empirical risk
#'   minimization \insertCite{chaudhuri2011}{DPpack}. Either the output or the
#'   objective perturbation method can be used. It is intended to be a framework
#'   for building more specific models via inheritance. See
#'   \code{\link{LogisticRegressionDP}} for an example of this type of
#'   structure.
#'
#' @details To use this class for empirical risk minimization, first use the
#'   \code{new} method to construct an object of this class with the desired
#'   function values and hyperparameters. After constructing the object, the
#'   \code{fit} method can be applied with a provided dataset and data bounds to
#'   fit the model.  In fitting, the model stores a vector of coefficients
#'   \code{coeff} which satisfy differential privacy. These can be released
#'   directly, or used in conjunction with the \code{predict} method to
#'   privately predict the outcomes of new datapoints.
#'
#'   Note that in order to guarantee differential privacy for empirical risk
#'   minimization, certain constraints must be satisfied for the values used to
#'   construct the object, as well as for the data used to fit. These conditions
#'   depend on the chosen perturbation method. Specifically, the provided loss
#'   function must be convex and differentiable with respect to \code{y.hat},
#'   and the absolute value of the first derivative of the loss function must be
#'   at most 1. If objective perturbation is chosen, the loss function must also
#'   be doubly differentiable and the absolute value of the second derivative of
#'   the loss function must be bounded above by a constant c for all possible
#'   values of \code{y.hat} and \code{y}, where \code{y.hat} is the predicted
#'   label and \code{y} is the true label. The regularizer must be 1-strongly
#'   convex and differentiable. It also must be doubly differentiable if
#'   objective perturbation is chosen. Finally, it is assumed that if x
#'   represents a single row of the dataset X, then the l2-norm of x is at most
#'   1 for all x. Note that because of this, a bias term cannot be included
#'   without appropriate scaling/preprocessing of the dataset. To ensure
#'   privacy, the add.bias argument in the \code{fit} and \code{predict} methods
#'   should only be utilized in subclasses within this package where appropriate
#'   preprocessing is implemented, not in this class.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @examples
#' # Build train dataset X and y, and test dataset Xtest and ytest
#' N <- 200
#' K <- 2
#' X <- data.frame()
#' y <- data.frame()
#' for (j in (1:K)){
#'   t <- seq(-.25, .25, length.out = N)
#'   if (j==1) m <- stats::rnorm(N,-.2, .1)
#'   if (j==2) m <- stats::rnorm(N, .2, .1)
#'   Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
#'   ytemp <- data.frame(matrix(j-1, N, 1))
#'   X <- rbind(X, Xtemp)
#'   y <- rbind(y, ytemp)
#' }
#' Xtest <- X[seq(1,(N*K),10),]
#' ytest <- y[seq(1,(N*K),10),,drop=FALSE]
#' X <- X[-seq(1,(N*K),10),]
#' y <- y[-seq(1,(N*K),10),,drop=FALSE]
#'
#' # Construct object for logistic regression
#' mapXy <- function(X, coeff) e1071::sigmoid(X%*%coeff)
#' # Cross entropy loss
#' loss <- function(y.hat,y) -(y*log(y.hat) + (1-y)*log(1-y.hat))
#' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
#' eps <- 1
#' gamma <- 1
#' perturbation.method <- 'objective'
#' c <- 1/4 # Required value for logistic regression
#' mapXy.gr <- function(X, coeff) as.numeric(e1071::dsigmoid(X%*%coeff))*t(X)
#' loss.gr <- function(y.hat, y) -y/y.hat + (1-y)/(1-y.hat)
#' regularizer.gr <- function(coeff) coeff
#' ermdp <- EmpiricalRiskMinimizationDP.CMS$new(mapXy, loss, regularizer, eps,
#'                                              gamma, perturbation.method, c,
#'                                              mapXy.gr, loss.gr,
#'                                              regularizer.gr)
#'
#' # Fit with data
#' # Bounds for X based on construction
#' upper.bounds <- c( 1, 1)
#' lower.bounds <- c(-1,-1)
#' ermdp$fit(X, y, upper.bounds, lower.bounds) # No bias term
#' ermdp$coeff # Gets private coefficients
#'
#' # Predict new data points
#' predicted.y <- ermdp$predict(Xtest)
#' n.errors <- sum(round(predicted.y)!=ytest)
#'
#' @keywords internal
#'
#' @export
EmpiricalRiskMinimizationDP.CMS <- R6::R6Class("EmpiricalRiskMinimizationDP.CMS",
  public=list(
  #' @field mapXy Map function of the form \code{mapXy(X, coeff)} mapping input
  #'   data matrix \code{X} and coefficient vector or matrix \code{coeff} to
  #'   output labels \code{y}.
  mapXy = NULL,
  #' @field mapXy.gr Function representing the gradient of the map function with
  #'   respect to the values in \code{coeff} and of the form \code{mapXy.gr(X,
  #'   coeff)}, where \code{X} is a matrix and \code{coeff} is a matrix or
  #'   numeric vector.
  mapXy.gr = NULL,
  #' @field loss Loss function of the form \code{loss(y.hat, y)}, where
  #'   \code{y.hat} and \code{y} are matrices.
  loss = NULL,
  #' @field loss.gr Function representing the gradient of the loss function with
  #'   respect to \code{y.hat} and of the form \code{loss.gr(y.hat, y)}, where
  #'   \code{y.hat} and \code{y} are matrices.
  loss.gr = NULL,
  #' @field regularizer Regularization function of the form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix.
  regularizer = NULL,
  #' @field regularizer.gr Function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}.
  regularizer.gr = NULL,
  #' @field gamma Nonnegative real number representing the regularization
  #'   constant.
  gamma = NULL,
  #' @field eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  eps = NULL,
  #' @field perturbation.method String indicating whether to use the 'output' or
  #'   the 'objective' perturbation methods \insertCite{chaudhuri2011}{DPpack}.
  perturbation.method = NULL,
  #' @field c Positive real number denoting the upper bound on the absolute
  #'   value of the second derivative of the loss function, as required to
  #'   ensure differential privacy for the objective perturbation method.
  c = NULL,
  #' @field coeff Numeric vector of coefficients for the model.
  coeff = NULL,
  #' @field kernel Value only used in child class \code{\link{svmDP}}. String
  #'   indicating which kernel to use for SVM. Must be one of \{'linear',
  #'   'Gaussian'\}. If 'linear' (default), linear SVM is used. If
  #'   'Gaussian', uses the sampling function corresponding to the Gaussian
  #'   (radial) kernel approximation.
  kernel = NULL,
  #' @field D Value only used in child class \code{\link{svmDP}}. Nonnegative
  #'   integer indicating the dimensionality of the transform space
  #'   approximating the kernel. Higher values of \code{D} provide better kernel
  #'   approximations at a cost of computational efficiency.
  D = NULL,
  #' @field sampling Value only used in child class \code{\link{svmDP}}.
  #'   Sampling function of the form \code{sampling(d)}, where \code{d} is the
  #'   input dimension, returning a (\code{d}+1)-dimensional vector of samples
  #'   corresponding to the Fourier transform of the kernel to be approximated.
  sampling=NULL,
  #' @field phi Value only used in child class \code{\link{svmDP}}. Function of
  #'   the form \code{phi(x, theta)}, where \code{x} is an individual row of the
  #'   original dataset, and theta is a (\code{d}+1)-dimensional vector sampled
  #'   from the Fourier transform of the kernel to be approximated, where
  #'   \code{d} is the dimension of \code{x}. The function returns a numeric
  #'   scalar corresponding to the pre-filtered value at the given row with the
  #'   given sampled vector.
  phi=NULL,
  #' @field kernel.param Value only used in child class \code{\link{svmDP}}.
  #'   Positive real number corresponding to the Gaussian kernel parameter.
  kernel.param=NULL,
  #' @field prefilter Value only used in child class \code{\link{svmDP}}. Matrix
  #'   of pre-filter values used in converting data into transform space.
  prefilter=NULL,
  #' @description Create a new \code{EmpiricalRiskMinimizationDP.CMS} object.
  #' @param mapXy Map function of the form \code{mapXy(X, coeff)} mapping input
  #'   data matrix \code{X} and coefficient vector or matrix \code{coeff} to
  #'   output labels \code{y}. Should return a column matrix of predicted labels
  #'   for each row of \code{X}. See \code{\link{mapXy.sigmoid}} for an example.
  #' @param loss Loss function of the form \code{loss(y.hat, y)}, where
  #'   \code{y.hat} and \code{y} are matrices. Should be defined such that it
  #'   returns a matrix of loss values for each element of \code{y.hat} and
  #'   \code{y}. See \code{\link{loss.cross.entropy}} for an example. It must be
  #'   convex and differentiable, and the absolute value of the first derivative
  #'   of the loss function must be at most 1. Additionally, if the objective
  #'   perturbation method is chosen, it must be doubly differentiable and the
  #'   absolute value of the second derivative of the loss function must be
  #'   bounded above by a constant c for all possible values of \code{y.hat} and
  #'   \code{y}.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix, and
  #'   return the value of the regularizer at \code{coeff}. See
  #'   \code{\link{regularizer.l2}} for an example. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   differentiable. If the objective perturbation method is chosen, it must
  #'   also be doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param perturbation.method String indicating whether to use the 'output' or
  #'   the 'objective' perturbation methods \insertCite{chaudhuri2011}{DPpack}.
  #'   Defaults to 'objective'.
  #' @param c Positive real number denoting the upper bound on the absolute
  #'   value of the second derivative of the loss function, as required to
  #'   ensure differential privacy for the objective perturbation method. This
  #'   input is unnecessary if perturbation.method is 'output', but is required
  #'   if perturbation.method is 'objective'. Defaults to NULL.
  #' @param mapXy.gr Optional function representing the gradient of the map
  #'   function with respect to the values in \code{coeff}. If given, must be of
  #'   the form \code{mapXy.gr(X, coeff)}, where \code{X} is a matrix and
  #'   \code{coeff} is a matrix or numeric vector. Should be defined such that
  #'   the ith row of the output represents the gradient with respect to the ith
  #'   coefficient. See \code{\link{mapXy.gr.sigmoid}} for an example. If not
  #'   given, non-gradient based optimization methods are used to compute the
  #'   coefficient values in fitting the model.
  #' @param loss.gr Optional function representing the gradient of the loss
  #'   function with respect to \code{y.hat} and of the form
  #'   \code{loss.gr(y.hat, y)}, where \code{y.hat} and \code{y} are matrices.
  #'   Should be defined such that the ith row of the output represents the
  #'   gradient of the loss function at the ith set of input values. See
  #'   \code{\link{loss.gr.cross.entropy}} for an example. If not given,
  #'   non-gradient based optimization methods are used to compute the
  #'   coefficient values in fitting the model.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}. Should return a vector. See
  #'   \code{\link{regularizer.gr.l2}} for an example. If \code{regularizer} is
  #'   given as a string, this value is ignored. If not given and
  #'   \code{regularizer} is a function, non-gradient based optimization methods
  #'   are used to compute the coefficient values in fitting the model.
  #'
  #' @return A new \code{EmpiricalRiskMinimizationDP.CMS} object.
  initialize = function(mapXy, loss, regularizer, eps, gamma,
                        perturbation.method = 'objective', c = NULL,
                        mapXy.gr = NULL, loss.gr = NULL, regularizer.gr = NULL){
    self$mapXy <- mapXy
    self$mapXy.gr <- mapXy.gr
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
    self$gamma <- gamma
    if (perturbation.method != 'output' & perturbation.method != 'objective'){
      stop("perturbation.method must be one of 'output' or 'objective'.")
    }
    self$perturbation.method <- perturbation.method
    if (perturbation.method == 'objective' & is.null(c)){
      stop("c must be given for 'objective' perturbation.method.")
    }
    self$c <- c
  },
  #' @description Fit the differentially private empirical risk minimization
  #'   model. This method runs either the output perturbation or the objective
  #'   perturbation algorithm \insertCite{chaudhuri2011}{DPpack}, depending on
  #'   the value of perturbation.method used to construct the object, to
  #'   generate an objective function. A numerical optimization method is then
  #'   run to find optimal coefficients for fitting the model given the training
  #'   data and hyperparameters. The built-in \code{\link{optim}} function using
  #'   the "BFGS" optimization method is used. If \code{mapXy.gr},
  #'   \code{loss.gr}, and \code{regularizer.gr} are all given in the
  #'   construction of the object, the gradient of the objective function is
  #'   utilized by \code{optim} as well. Otherwise, non-gradient based
  #'   optimization methods are used. The resulting privacy-preserving
  #'   coefficients are stored in \code{coeff}.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true labels for each row of \code{X}.
  #' @param upper.bounds Numeric vector of length \code{ncol(X)} giving upper
  #'   bounds on the values in each column of X. The \code{ncol(X)} values are
  #'   assumed to be in the same order as the corresponding columns of \code{X}.
  #'   Any value in the columns of \code{X} larger than the corresponding upper
  #'   bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length \code{ncol(X)} giving lower
  #'   bounds on the values in each column of \code{X}. The \code{ncol(X)}
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of \code{X}. Any value in the columns of \code{X} larger than the
  #'   corresponding upper bound is clipped at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE.
  fit = function(X, y, upper.bounds, lower.bounds, add.bias=FALSE){
    preprocess <- private$preprocess_data(X,y,upper.bounds,lower.bounds,
                                            add.bias)
    X <- preprocess$X
    y <- preprocess$y

    n <- length(y)
    d <- ncol(X)
    if (!is.infinite(self$eps) & self$perturbation.method=='objective'){
      # eps.prime <- self$eps - log(1 + 2*self$c/(n*self$gamma) +
      #                               self$c^2/(n^2*self$gamma^2))
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      eps.prime <- self$eps - log(1 + 2*self$c/self$gamma +
                                    self$c^2/self$gamma^2)
      if (eps.prime > 0) {
        Delta <- 0
      }
      else {
        # Delta <- self$c/(n*(exp(self$eps/4) - 1)) - self$gamma
        # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
        Delta <- self$c/(n*(exp(self$eps/4) - 1)) - self$gamma/n
        eps.prime <- self$eps/2
      }
      beta <- eps.prime/2
      norm.b <- rgamma(1, d, rate=beta)
      direction.b <- stats::rnorm(d);
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b * direction.b
    } else {
      Delta <- 0
      b <- numeric(d)
    }

    tmp.coeff <- private$optimize_coeff(X, y, Delta, b)
    if (self$perturbation.method=='output'){
      # beta <- n*self$gamma*self$eps/2
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      beta <- self$gamma*self$eps/2
      norm.b <- rgamma(1, d, rate=beta)
      direction.b <- stats::rnorm(d)
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b * direction.b
      tmp.coeff <- tmp.coeff + b
    }
    self$coeff <- private$postprocess_coeff(tmp.coeff, preprocess)
    invisible(self)
  },
  #' @description Predict label(s) for given \code{X} using the fitted
  #'   coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as \code{X} used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #'
  #' @return Matrix of predicted labels corresponding to each row of \code{X}.
  predict = function(X, add.bias=FALSE){
    if (add.bias){
      X <- dplyr::mutate(X, bias=1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
    }
    self$mapXy(as.matrix(X), self$coeff)
  }
), private = list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X) giving upper bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param lower.bounds Numeric vector of length ncol(X) giving lower bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # param weights Numeric vector of observation weights of the same length as
  #   \code{y}.
  # param weights.upper.bound Numeric value representing the global or public
  #   upper bound on the weights.
  #
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving empirical risk
  #   minimization algorithm.
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias,
                             weights=NULL, weights.upper.bound=NULL){
    # Make sure values are correct
    if (length(upper.bounds)!=ncol(X)) {
      stop("Length of upper.bounds must be equal to the number of columns of X.");
    }
    if (length(lower.bounds)!=ncol(X)) {
      stop("Length of lower.bounds must be equal to the number of columns of X.");
    }

    # Add bias if needed (highly recommended to not do this due to unwanted
    #       regularization of bias term)
    if (add.bias){
      X <- dplyr::mutate(X, bias=1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
      upper.bounds <- c(1,upper.bounds)
      lower.bounds <- c(1,lower.bounds)
    }

    # Make matrices for multiplication purposes
    X <- as.matrix(X,nrow=length(y))
    y <- as.matrix(y,ncol=1)

    # Clip based on provided bounds
    for (i in 1:length(upper.bounds)){
      X[X[,i]>upper.bounds[i],i] <- upper.bounds[i]
      X[X[,i]<lower.bounds[i],i] <- lower.bounds[i]
    }

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
  #   built-in optim function using the "BFGS" optimization method. If mapXy.gr,
  #   loss.gr, and regularizer.gr are all given in the construction of the
  #   object, the gradient of the objective function is utilized by optim as
  #   well. Otherwise, non-gradient based methods are used.
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
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      as.numeric(sum(self$loss(self$mapXy(X,par),y))/n +
                   self$gamma*self$regularizer(par)/n + t(b)%*%par/n +
                   Delta*par%*%par/2)
    }

    # Get gradient function
    if (!is.null(self$mapXy.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      objective.gr <- function(par, X, y, Delta, b){
        # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
        as.numeric(self$mapXy.gr(X,par)%*%self$loss.gr(self$mapXy(X,par),y)/n +
                     self$gamma*self$regularizer.gr(par)/n + b/n + Delta*par)
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

#' Privacy-preserving Weighted Empirical Risk Minimization
#'
#' @description This class implements differentially private empirical risk
#'   minimization in the case where weighted observation-level losses are
#'   desired (such as weighted SVM \insertCite{Yang2005}{DPpack}). Currently,
#'   only the output perturbation method is implemented.
#'
#' @details To use this class for weighted empirical risk minimization, first
#'   use the \code{new} method to construct an object of this class with the
#'   desired function values and hyperparameters. After constructing the object,
#'   the \code{fit} method can be applied with a provided dataset, data bounds,
#'   weights, and weight bounds to fit the model. In fitting, the model stores a
#'   vector of coefficients \code{coeff} which satisfy differential privacy.
#'   These can be released directly, or used in conjunction with the
#'   \code{predict} method to privately predict the outcomes of new datapoints.
#'
#'   Note that in order to guarantee differential privacy for weighted empirical
#'   risk minimization, certain constraints must be satisfied for the values
#'   used to construct the object, as well as for the data used to fit. These
#'   conditions depend on the chosen perturbation method, though currently only
#'   output perturbation is implemented. Specifically, the provided loss
#'   function must be convex and differentiable with respect to \code{y.hat},
#'   and the absolute value of the first derivative of the loss function must be
#'   at most 1. If objective perturbation is chosen (not currently implemented),
#'   the loss function must also be doubly differentiable and the absolute value
#'   of the second derivative of the loss function must be bounded above by a
#'   constant c for all possible values of \code{y.hat} and \code{y}, where
#'   \code{y.hat} is the predicted label and \code{y} is the true label. The
#'   regularizer must be 1-strongly convex and differentiable. It also must be
#'   doubly differentiable if objective perturbation is chosen. For the data x,
#'   it is assumed that if x represents a single row of the dataset X, then the
#'   l2-norm of x is at most 1 for all x. Note that because of this, a bias term
#'   cannot be included without appropriate scaling/preprocessing of the
#'   dataset. To ensure privacy, the add.bias argument in the \code{fit} and
#'   \code{predict} methods should only be utilized in subclasses within this
#'   package where appropriate preprocessing is implemented, not in this class.
#'   Finally, if weights are provided, they should be nonnegative, of the same
#'   length as y, and be upper bounded by a global or public bound which must
#'   also be provided.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#'   \insertRef{Yang2005}{DPpack}
#'
#' @examples
#' # Build train dataset X and y, and test dataset Xtest and ytest
#' N <- 200
#' K <- 2
#' X <- data.frame()
#' y <- data.frame()
#' for (j in (1:K)){
#'   t <- seq(-.25, .25, length.out = N)
#'   if (j==1) m <- stats::rnorm(N,-.2, .1)
#'   if (j==2) m <- stats::rnorm(N, .2, .1)
#'   Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
#'   ytemp <- data.frame(matrix(j-1, N, 1))
#'   X <- rbind(X, Xtemp)
#'   y <- rbind(y, ytemp)
#' }
#' Xtest <- X[seq(1,(N*K),10),]
#' ytest <- y[seq(1,(N*K),10),,drop=FALSE]
#' X <- X[-seq(1,(N*K),10),]
#' y <- y[-seq(1,(N*K),10),,drop=FALSE]
#'
#' # Construct object for weighted linear SVM
#' mapXy <- function(X, coeff) X%*%coeff
#' # Huber loss from DPpack
#' huber.h <- 0.5
#' loss <- generate.loss.huber(huber.h)
#' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
#' eps <- 1
#' gamma <- 1
#' perturbation.method <- 'output'
#' c <- 1/(2*huber.h) # Required value for SVM
#' mapXy.gr <- function(X, coeff) t(X)
#' loss.gr <- generate.loss.gr.huber(huber.h)
#' regularizer.gr <- function(coeff) coeff
#' wermdp <- WeightedERMDP.CMS$new(mapXy, loss, regularizer, eps,
#'                                 gamma, perturbation.method, c,
#'                                 mapXy.gr, loss.gr,
#'                                 regularizer.gr)
#'
#' # Fit with data
#' # Bounds for X based on construction
#' upper.bounds <- c( 1, 1)
#' lower.bounds <- c(-1,-1)
#' weights <- rep(1, nrow(y)) # Uniform weighting
#' weights[nrow(y)] <- 0.5 # half weight for last observation
#' wub <- 1 # Public upper bound for weights
#' wermdp$fit(X, y, upper.bounds, lower.bounds, weights=weights,
#'            weights.upper.bound=wub)
#' wermdp$coeff # Gets private coefficients
#'
#' # Predict new data points
#' predicted.y <- wermdp$predict(Xtest)
#' n.errors <- sum(round(predicted.y)!=ytest)
#'
#' @keywords internal
#'
#' @export
WeightedERMDP.CMS <- R6::R6Class("WeightedERMDP.CMS",
  inherit=EmpiricalRiskMinimizationDP.CMS,
  public=list(
  #' @description Create a new \code{WeightedERMDP.CMS} object.
  #' @param mapXy Map function of the form \code{mapXy(X, coeff)} mapping input
  #'   data matrix \code{X} and coefficient vector or matrix \code{coeff} to
  #'   output labels \code{y}. Should return a column matrix of predicted labels
  #'   for each row of \code{X}. See \code{\link{mapXy.sigmoid}} for an example.
  #' @param loss Loss function of the form \code{loss(y.hat, y, w)}, where
  #'   \code{y.hat} and \code{y} are matrices and \code{w} is a matrix or vector
  #'   of weights of the same length as \code{y}. Should be defined such that it
  #'   returns a matrix of weighted loss values for each element of \code{y.hat}
  #'   and \code{y}. If \code{w} is not given, the function should operate as if
  #'   uniform weights were given. See \code{\link{generate.loss.huber}} for an
  #'   example. It must be convex and differentiable, and the absolute value of
  #'   the first derivative of the loss function must be at most 1.
  #'   Additionally, if the objective perturbation method is chosen, it must be
  #'   doubly differentiable and the absolute value of the second derivative of
  #'   the loss function must be bounded above by a constant c for all possible
  #'   values of \code{y.hat} and \code{y}.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix, and
  #'   return the value of the regularizer at \code{coeff}. See
  #'   \code{\link{regularizer.l2}} for an example. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   differentiable. If the objective perturbation method is chosen, it must
  #'   also be doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param perturbation.method String indicating whether to use the 'output' or
  #'   the 'objective' perturbation methods \insertCite{chaudhuri2011}{DPpack}.
  #'   Defaults to 'objective'. Currently, only the output perturbation method
  #'   is supported.
  #' @param c Positive real number denoting the upper bound on the absolute
  #'   value of the second derivative of the loss function, as required to
  #'   ensure differential privacy for the objective perturbation method. This
  #'   input is unnecessary if perturbation.method is 'output', but is required
  #'   if perturbation.method is 'objective'. Defaults to NULL.
  #' @param mapXy.gr Optional function representing the gradient of the map
  #'   function with respect to the values in \code{coeff}. If given, must be of
  #'   the form \code{mapXy.gr(X, coeff)}, where \code{X} is a matrix and
  #'   \code{coeff} is a matrix or numeric vector. Should be defined such that
  #'   the ith row of the output represents the gradient with respect to the ith
  #'   coefficient. See \code{\link{mapXy.gr.sigmoid}} for an example. If not
  #'   given, non-gradient based optimization methods are used to compute the
  #'   coefficient values in fitting the model.
  #' @param loss.gr Optional function representing the gradient of the loss
  #'   function with respect to \code{y.hat} and of the form
  #'   \code{loss.gr(y.hat, y, w)}, where \code{y.hat} and \code{y} are matrices
  #'   and \code{w} is a matrix or vector of weights. Should be defined such
  #'   that the ith row of the output represents the gradient of the (possibly
  #'   weighted) loss function at the ith set of input values. See
  #'   \code{\link{generate.loss.gr.huber}} for an example. If not given,
  #'   non-gradient based optimization methods are used to compute the
  #'   coefficient values in fitting the model.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}. Should return a vector. See
  #'   \code{\link{regularizer.gr.l2}} for an example. If \code{regularizer} is
  #'   given as a string, this value is ignored. If not given and
  #'   \code{regularizer} is a function, non-gradient based optimization methods
  #'   are used to compute the coefficient values in fitting the model.
  #'
  #' @return A new \code{WeightedERMDP.CMS} object.
  initialize = function(mapXy, loss, regularizer, eps, gamma,
                        perturbation.method = 'objective', c = NULL,
                        mapXy.gr = NULL, loss.gr = NULL, regularizer.gr = NULL){
    super$initialize(mapXy, loss, regularizer, eps, gamma, perturbation.method,
                     c, mapXy.gr, loss.gr, regularizer.gr)
  },
  #' @description Fit the differentially private weighted empirical risk
  #'   minimization model. This method runs either the output perturbation or
  #'   the objective perturbation algorithm \insertCite{chaudhuri2011}{DPpack}
  #'   (only output is currently implemented), depending on the value of
  #'   perturbation.method used to construct the object, to generate an
  #'   objective function. A numerical optimization method is then run to find
  #'   optimal coefficients for fitting the model given the training data,
  #'   weights, and hyperparameters. The built-in \code{\link{optim}} function
  #'   using the "BFGS" optimization method is used. If \code{mapXy.gr},
  #'   \code{loss.gr}, and \code{regularizer.gr} are all given in the
  #'   construction of the object, the gradient of the objective function is
  #'   utilized by \code{optim} as well. Otherwise, non-gradient based
  #'   optimization methods are used. The resulting privacy-preserving
  #'   coefficients are stored in \code{coeff}.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true labels for each row of \code{X}.
  #' @param upper.bounds Numeric vector of length \code{ncol(X)} giving upper
  #'   bounds on the values in each column of X. The \code{ncol(X)} values are
  #'   assumed to be in the same order as the corresponding columns of \code{X}.
  #'   Any value in the columns of \code{X} larger than the corresponding upper
  #'   bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length \code{ncol(X)} giving lower
  #'   bounds on the values in each column of \code{X}. The \code{ncol(X)}
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of \code{X}. Any value in the columns of \code{X} larger than the
  #'   corresponding upper bound is clipped at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE.
  #' @param weights Numeric vector of observation weights of the same length as
  #'   \code{y}.
  #' @param weights.upper.bound Numeric value representing the global or public
  #'   upper bound on the weights.
  fit = function(X, y, upper.bounds, lower.bounds, add.bias=FALSE, weights=NULL,
                 weights.upper.bound=NULL){
    if (!is.null(weights) & self$perturbation.method!="output"){
      stop(paste("Currently, weighting ERM samples only tested for SVM with output ",
                 "perturbation and should not be used with objective perturbation."))
    }
    preprocess <- private$preprocess_data(X,y,upper.bounds,lower.bounds,
                                          add.bias,weights,weights.upper.bound)
    X <- preprocess$X
    y <- preprocess$y
    weights <- preprocess$weights
    if (is.null(weights)) weights.upper.bound <- 1

    n <- length(y)
    d <- ncol(X)
    if (!is.infinite(self$eps) & self$perturbation.method=='objective'){
      # eps.prime <- self$eps - log(1 + 2*self$c/(n*self$gamma) +
      #                               self$c^2/(n^2*self$gamma^2))
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      eps.prime <- self$eps - log(1 + 2*self$c/self$gamma +
                                    self$c^2/self$gamma^2)
      if (eps.prime > 0) {
        Delta <- 0
      }
      else {
        # Delta <- self$c/(n*(exp(self$eps/4) - 1)) - self$gamma
        # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
        Delta <- self$c/(n*(exp(self$eps/4) - 1)) - self$gamma/n
        eps.prime <- self$eps/2
      }
      beta <- eps.prime/2
      norm.b <- rgamma(1, d, rate=beta)
      direction.b <- stats::rnorm(d);
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b * direction.b
    } else {
      Delta <- 0
      b <- numeric(d)
    }

    tmp.coeff <- private$optimize_coeff(X, y, Delta, b, weights)
    if (self$perturbation.method=='output'){
      # beta <- n*self$gamma*self$eps/(2*weights.upper.bound)
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      beta <- self$gamma*self$eps/(2*weights.upper.bound)
      norm.b <- rgamma(1, d, rate=beta)
      direction.b <- stats::rnorm(d)
      direction.b <- direction.b/sqrt(sum(direction.b^2))
      b <- norm.b * direction.b
      tmp.coeff <- tmp.coeff + b
    }
    self$coeff <- private$postprocess_coeff(tmp.coeff, preprocess)
    invisible(self)
  },
  #' @description Predict label(s) for given \code{X} using the fitted
  #'   coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as \code{X} used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #'
  #' @return Matrix of predicted labels corresponding to each row of \code{X}.
  predict = function(X, add.bias=FALSE){
    super$predict(X, add.bias)
  }
), private=list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X) giving upper bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param lower.bounds Numeric vector of length ncol(X) giving lower bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # param weights Numeric vector of observation weights of the same length as
  #   \code{y}.
  # param weights.upper.bound Numeric value representing the global or public
  #   upper bound on the weights.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving empirical risk
  #   minimization algorithm.
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias,
                             weights=NULL, weights.upper.bound=NULL){
    res <- super$preprocess_data(X, y, upper.bounds, lower.bounds, add.bias)

    # Handle weights
    if (!is.null(weights)){
      if (is.null(weights.upper.bound)){
        stop("Upper bound on weights must be given if weights vector given.")
      }

      if (any(weights<0)) stop("Weights must be nonnegative.")

      weights[weights>weights.upper.bound] <- weights.upper.bound
    }
    list(X=res$X, y=res$y, upper.bounds=res$upper.bounds,
         lower.bounds=res$lower.bounds, weights=weights)
  },
  # description Postprocess coefficients obtained by fitting the model using
  #   differential privacy to ensure they match original inputs X and y.
  #   Effectively undoes the processing done in preprocess_data.
  # param coeff Vector of coefficients obtained by fitting the model.
  # param preprocess List of values returned by preprocess_data and used to
  #   process coeff.
  # return The processed coefficient vector.
  postprocess_coeff = function(coeff, preprocess){
    super$postprocess_coeff(coeff, preprocess)
  },
  # description Run numerical optimization method to find optimal coefficients
  #   for fitting model given training data, observation weights, and
  #   hyperparameters. This function builds the objective function based on the
  #   training data and observation weights, and runs the built-in optim function
  #   using the "BFGS" optimization method. If mapXy.gr, loss.gr, and
  #   regularizer.gr are all given in the construction of the object, the gradient
  #   of the objective function is utilized by optim as well. Otherwise,
  #   non-gradient based methods are used.
  # param X Matrix of data to be fit.
  # param y Vector or matrix of true labels for each row of X.
  # param Delta Nonnegative real number set by the differentially private
  #   algorithm and required for constructing the objective function.
  # param b Numeric perturbation vector randomly drawn by the differentially
  #   private algorithm and required for constructing the objective function.
  # param weights Numeric vector of observation weights of the same length as
  #   \code{y}.
  # return Vector of fitted coefficients.
  optimize_coeff=function(X, y, Delta, b, weights){
    n <- length(y)
    d <- ncol(X)

    # Get objective function
    objective <- function(par, X, y, Delta, b, weights){
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      as.numeric(sum(self$loss(self$mapXy(X,par),y, weights))/n +
                   self$gamma*self$regularizer(par)/n + t(b)%*%par/n +
                   Delta*par%*%par/2)
    }

    # Get gradient function
    if (!is.null(self$mapXy.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      objective.gr <- function(par, X, y, Delta, b, weights){
        as.numeric(self$mapXy.gr(X,par)%*%
                     self$loss.gr(self$mapXy(X,par), y, weights)/n +
                     self$gamma*self$regularizer.gr(par)/n + b/n + Delta*par)
      }
    }
    else objective.gr <- NULL

    # Run optimization
    coeff0 <- numeric(ncol(X))
    opt.res <- optim(coeff0, fn=objective, gr=objective.gr, method="BFGS",
                     X=X, y=y, Delta=Delta, b=b, weights=weights)
    opt.res$par
  }
))

#' Privacy-preserving Logistic Regression
#'
#' @description This class implements differentially private logistic regression
#'   \insertCite{chaudhuri2011}{DPpack}. Either the output or the objective
#'   perturbation method can be used.
#'
#' @details To use this class for logistic regression, first use the \code{new}
#'   method to construct an object of this class with the desired function
#'   values and hyperparameters. After constructing the object, the \code{fit}
#'   method can be applied with a provided dataset and data bounds to fit the
#'   model. In fitting, the model stores a vector of coefficients \code{coeff}
#'   which satisfy differential privacy. These can be released directly, or used
#'   in conjunction with the \code{predict} method to privately predict the
#'   outcomes of new datapoints.
#'
#'   Note that in order to guarantee differential privacy for logistic
#'   regression, certain constraints must be satisfied for the values used to
#'   construct the object, as well as for the data used to fit. These conditions
#'   depend on the chosen perturbation method. The regularizer must be
#'   1-strongly convex and differentiable. It also must be doubly differentiable
#'   if objective perturbation is chosen. Additionally, it is assumed that if x
#'   represents a single row of the dataset X, then the l2-norm of x is at most
#'   1 for all x. In order to ensure this constraint is satisfied, the dataset
#'   is preprocessed and scaled, and the resulting coefficients are
#'   postprocessed and un-scaled so that the stored coefficients correspond to
#'   the original data. Due to this constraint on x, it is best to avoid using a
#'   bias term in the model whenever possible. If a bias term must be used, the
#'   issue can be partially circumvented by adding a constant column to X before
#'   fitting the model, which will be scaled along with the rest of X. The
#'   \code{fit} method contains functionality to add a column of constant 1s to
#'   X before scaling, if desired.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#'   \insertRef{Chaudhuri2009}{DPpack}
#'
#' @examples
#' # Build train dataset X and y, and test dataset Xtest and ytest
#' N <- 200
#' K <- 2
#' X <- data.frame()
#' y <- data.frame()
#' for (j in (1:K)){
#'   t <- seq(-.25, .25, length.out = N)
#'   if (j==1) m <- stats::rnorm(N,-.2, .1)
#'   if (j==2) m <- stats::rnorm(N, .2, .1)
#'   Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
#'   ytemp <- data.frame(matrix(j-1, N, 1))
#'   X <- rbind(X, Xtemp)
#'   y <- rbind(y, ytemp)
#' }
#' Xtest <- X[seq(1,(N*K),10),]
#' ytest <- y[seq(1,(N*K),10),,drop=FALSE]
#' X <- X[-seq(1,(N*K),10),]
#' y <- y[-seq(1,(N*K),10),,drop=FALSE]
#'
#' # Construct object for logistic regression
#' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
#' eps <- 1
#' gamma <- 1
#' lrdp <- LogisticRegressionDP$new(regularizer, eps, gamma)
#'
#' # Fit with data
#' # Bounds for X based on construction
#' upper.bounds <- c( 1, 1)
#' lower.bounds <- c(-1,-1)
#' lrdp$fit(X, y, upper.bounds, lower.bounds) # No bias term
#' lrdp$coeff # Gets private coefficients
#'
#' # Predict new data points
#' predicted.y <- lrdp$predict(Xtest)
#' n.errors <- sum(predicted.y!=ytest)
#'
#' @export
LogisticRegressionDP <- R6::R6Class("LogisticRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.CMS,
  public=list(
  #' @description Create a new \code{LogisticRegressionDP} object.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix, and
  #'   return the value of the regularizer at \code{coeff}. See
  #'   \code{\link{regularizer.l2}} for an example. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param perturbation.method String indicating whether to use the 'output' or
  #'   the 'objective' perturbation methods \insertCite{chaudhuri2011}{DPpack}.
  #'   Defaults to 'objective'.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}. Should return a vector. See
  #'   \code{\link{regularizer.gr.l2}} for an example. If \code{regularizer} is
  #'   given as a string, this value is ignored. If not given and
  #'   \code{regularizer} is a function, non-gradient based optimization methods
  #'   are used to compute the coefficient values in fitting the model.
  #'
  #' @return A new \code{LogisticRegressionDP} object.
  initialize = function(regularizer, eps, gamma,
                        perturbation.method = 'objective',
                        regularizer.gr = NULL){
    super$initialize(mapXy.sigmoid, loss.cross.entropy, regularizer, eps,
                    gamma, perturbation.method, 1/4, mapXy.gr.sigmoid,
                    loss.gr.cross.entropy, regularizer.gr)
  },
  #' @description Fit the differentially private logistic regression model. This
  #'   method runs either the output perturbation or the objective perturbation
  #'   algorithm \insertCite{chaudhuri2011}{DPpack}, depending on the value of
  #'   perturbation.method used to construct the object, to generate an
  #'   objective function. A numerical optimization method is then run to find
  #'   optimal coefficients for fitting the model given the training data and
  #'   hyperparameters. The built-in \code{\link{optim}} function using the
  #'   "BFGS" optimization method is used. If \code{regularizer} is given as
  #'   'l2' or if \code{regularizer.gr} is given in the construction of the
  #'   object, the gradient of the objective function is utilized by
  #'   \code{optim} as well. Otherwise, non-gradient based optimization methods
  #'   are used. The resulting privacy-preserving coefficients are stored in
  #'   \code{coeff}.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true labels for each row of \code{X}.
  #' @param upper.bounds Numeric vector of length \code{ncol(X)} giving upper
  #'   bounds on the values in each column of X. The \code{ncol(X)} values are
  #'   assumed to be in the same order as the corresponding columns of \code{X}.
  #'   Any value in the columns of \code{X} larger than the corresponding upper
  #'   bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length \code{ncol(X)} giving lower
  #'   bounds on the values in each column of \code{X}. The \code{ncol(X)}
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of \code{X}. Any value in the columns of \code{X} larger than the
  #'   corresponding upper bound is clipped at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE.
  fit = function(X, y, upper.bounds, lower.bounds, add.bias=FALSE){
    super$fit(X,y,upper.bounds,lower.bounds,add.bias)
  },
  #' @description Predict label(s) for given \code{X} using the fitted
  #'   coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as \code{X} used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #' @param raw.value Boolean indicating whether to return the raw predicted
  #'   value or the rounded class label. If FALSE (default), outputs the
  #'   predicted labels 0 or 1. If TRUE, returns the raw score from the logistic
  #'   regression.
  #'
  #' @return Matrix of predicted labels or scores corresponding to each row of
  #'   \code{X}.
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
  # param upper.bounds Numeric vector of length ncol(X) giving upper bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param lower.bounds Numeric vector of length ncol(X) giving lower bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving logistic regression
  #   algorithm.
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias){
    res <- super$preprocess_data(X, y, upper.bounds, lower.bounds, add.bias)
    X <- res$X
    y <- res$y
    upper.bounds <- res$upper.bounds
    lower.bounds <- res$lower.bounds

    # Process X
    p <- length(lower.bounds)
    tmp <- matrix(c(lower.bounds,upper.bounds),ncol=2)
    scale1 <- apply(abs(tmp),1,max)
    scale2 <- sqrt(p)
    X.norm <- t(t(X)/scale1)/scale2

    # Process y.
    # Labels must be 0 and 1
    if (any((y!=1) & (y!=0))) stop("y must be a numeric vector of binary labels
                                    0 and 1.")

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
  #   built-in optim function using the "BFGS" optimization method. If mapXy.gr,
  #   loss.gr, and regularizer.gr are all given in the construction of the
  #   object, the gradient of the objective function is utilized by optim as
  #   well. Otherwise, non-gradient based optimization methods are used.
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
    objective <- function(par, X, y, gamma, Delta, b){
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      as.numeric(sum(self$loss(self$mapXy(X,par),y))/n +
                   gamma*self$regularizer(par)/n + t(b)%*%par/n +
                   Delta*par%*%par/2)
    }

    # Get gradient function
    if (!is.null(self$mapXy.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      # NOTE: Gamma from Chaudhuri's paper is gamma/n in this implementation
      objective.gr <- function(par, X, y, gamma, Delta, b){
        as.numeric(t(X)%*%(self$mapXy(X, par)-y)/n +
                     gamma*self$regularizer.gr(par)/n + b/n + Delta*par)
      }
    }
    else objective.gr <- NULL

    # Run optimization
    coeff0 <- numeric(ncol(X))
    opt.res <- optim(coeff0, fn=objective, gr=objective.gr, method="BFGS",
                     X=X, y=y, gamma=gamma, Delta=Delta, b=b)
    opt.res$par
  }
))

#' Generator for Sampling Distribution Function for Gaussian Kernel
#'
#' This function generates and returns a sampling function corresponding to the
#' Fourier transform of a Gaussian kernel with given parameter
#' \insertCite{chaudhuri2011}{DPpack} of form needed for \code{\link{svmDP}}.
#'
#' @param kernel.param Positive real number for the Gaussian (radial) kernel
#'   parameter.
#' @return Sampling function for the Gaussian kernel of form required by
#'   \code{\link{svmDP}}.
#' @examples
#'   kernel.param <- 0.1
#'   sample <- generate.sampling(kernel.param)
#'   d <- 5
#'   sample(d)
#'
#' @keywords internal
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#' @export
generate.sampling <- function(kernel.param){
  function(d){
    omega <- stats::rnorm(d,sd=sqrt(2*kernel.param))
    phi <- stats::runif(1,-pi,pi)
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

#' Privacy-preserving Support Vector Machine
#'
#' @description This class implements differentially private support vector
#'   machine (SVM) \insertCite{chaudhuri2011}{DPpack}. It can be either weighted
#'   \insertCite{Yang2005}{DPpack} or unweighted. Either the output or the
#'   objective perturbation method can be used for unweighted SVM, though only
#'   the output perturbation method is currently supported for weighted SVM.
#'
#' @details To use this class for SVM, first use the \code{new} method to
#'   construct an object of this class with the desired function values and
#'   hyperparameters, including a choice of the desired kernel. After
#'   constructing the object, the \code{fit} method can be applied to fit the
#'   model with a provided dataset, data bounds, and optional observation
#'   weights and weight upper bound. In fitting, the model stores a vector of
#'   coefficients \code{coeff} which satisfy differential privacy. Additionally,
#'   if a nonlinear kernel is chosen, the models stores a mapping function from
#'   the input data X to a higher dimensional embedding V in the form of a
#'   method \code{XtoV} as required \insertCite{chaudhuri2011}{DPpack}. These
#'   can be released directly, or used in conjunction with the \code{predict}
#'   method to privately predict the label of new datapoints. Note that the
#'   mapping function \code{XtoV} is based on an approximation method via
#'   Fourier transforms \insertCite{Rahimi2007,Rahimi2008}{DPpack}.
#'
#'   Note that in order to guarantee differential privacy for the SVM model,
#'   certain constraints must be satisfied for the values used to construct the
#'   object, as well as for the data used to fit. These conditions depend on the
#'   chosen perturbation method. First, the loss function is assumed to be
#'   differentiable (and doubly differentiable if the objective perturbation
#'   method is used). The hinge loss, which is typically used for SVM, is not
#'   differentiable at 1. Thus, to satisfy this constraint, this class utilizes
#'   the Huber loss, a smooth approximation to the hinge loss
#'   \insertCite{Chapelle2007}{DPpack}. The level of approximation to the hinge
#'   loss is determined by a user-specified constant, h, which defaults to 0.5,
#'   a typical value. Additionally, the regularizer must be 1-strongly convex
#'   and differentiable. It also must be doubly differentiable if objective
#'   perturbation is chosen. If weighted SVM is desired, the provided weights
#'   must be nonnegative and bounded above by a global or public value, which
#'   must also be provided.
#'
#'   Finally, it is assumed that if x represents a single row of the dataset X,
#'   then the l2-norm of x is at most 1 for all x. In order to ensure this
#'   constraint is satisfied, the dataset is preprocessed and scaled, and the
#'   resulting coefficients are postprocessed and un-scaled so that the stored
#'   coefficients correspond to the original data. Due to this constraint on x,
#'   it is best to avoid using a bias term in the model whenever possible. If a
#'   bias term must be used, the issue can be partially circumvented by adding a
#'   constant column to X before fitting the model, which will be scaled along
#'   with the rest of X. The \code{fit} method contains functionality to add a
#'   column of constant 1s to X before scaling, if desired.
#'
#' @references \insertRef{chaudhuri2011}{DPpack}
#'
#'   \insertRef{Yang2005}{DPpack}
#'
#'   \insertRef{Chapelle2007}{DPpack}
#'
#'   \insertRef{Rahimi2007}{DPpack}
#'
#'   \insertRef{Rahimi2008}{DPpack}
#'
#' @examples
#' # Build train dataset X and y, and test dataset Xtest and ytest
#' N <- 400
#' X <- data.frame()
#' y <- data.frame()
#' for (i in (1:N)){
#'   Xtemp <- data.frame(x1 = stats::rnorm(1,sd=.28) , x2 = stats::rnorm(1,sd=.28))
#'   if (sum(Xtemp^2)<.15) ytemp <- data.frame(y=0)
#'   else ytemp <- data.frame(y=1)
#'   X <- rbind(X, Xtemp)
#'   y <- rbind(y, ytemp)
#' }
#' Xtest <- X[seq(1,N,10),]
#' ytest <- y[seq(1,N,10),,drop=FALSE]
#' X <- X[-seq(1,N,10),]
#' y <- y[-seq(1,N,10),,drop=FALSE]
#'
#' # Construct object for SVM
#' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
#' eps <- 1
#' gamma <- 1
#' perturbation.method <- 'output'
#' kernel <- 'Gaussian'
#' D <- 20
#' svmdp <- svmDP$new(regularizer, eps, gamma, perturbation.method,
#'                    kernel=kernel, D=D)
#'
#' # Fit with data
#' # Bounds for X based on construction
#' upper.bounds <- c( 1, 1)
#' lower.bounds <- c(-1,-1)
#' weights <- rep(1, nrow(y)) # Uniform weighting
#' weights[nrow(y)] <- 0.5 # Half weight for last observation
#' wub <- 1 # Public upper bound for weights
#' svmdp$fit(X, y, upper.bounds, lower.bounds, weights=weights,
#'           weights.upper.bound=wub) # No bias term
#'
#' # Predict new data points
#' predicted.y <- svmdp$predict(Xtest)
#' n.errors <- sum(predicted.y!=ytest)
#'
#' @export
svmDP <- R6::R6Class("svmDP",
  inherit=WeightedERMDP.CMS,
  public=list(
  #' @description Create a new \code{svmDP} object.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix, and
  #'   return the value of the regularizer at \code{coeff}. See
  #'   \code{\link{regularizer.l2}} for an example. Additionally, in order to
  #'   ensure differential privacy, the function must be 1-strongly convex and
  #'   doubly differentiable.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param perturbation.method String indicating whether to use the 'output' or
  #'   the 'objective' perturbation methods \insertCite{chaudhuri2011}{DPpack}.
  #'   Defaults to 'objective'.
  #' @param kernel String indicating which kernel to use for SVM. Must be one of
  #'   \{'linear', 'Gaussian'\}. If 'linear' (default), linear SVM is used. If
  #'   'Gaussian,' uses the sampling function corresponding to the Gaussian
  #'   (radial) kernel approximation.
  #' @param D Nonnegative integer indicating the dimensionality of the transform
  #'   space approximating the kernel if a nonlinear kernel is used. Higher
  #'   values of D provide better kernel approximations at a cost of
  #'   computational efficiency. This value must be specified if a nonlinear
  #'   kernel is used.
  #' @param kernel.param Positive real number corresponding to the Gaussian
  #'   kernel parameter. Defaults to 1/p, where p is the number of predictors.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}. Should return a vector. See
  #'   \code{\link{regularizer.gr.l2}} for an example. If \code{regularizer} is
  #'   given as a string, this value is ignored. If not given and
  #'   \code{regularizer} is a function, non-gradient based optimization methods
  #'   are used to compute the coefficient values in fitting the model.
  #' @param huber.h Positive real number indicating the degree to which the
  #'   Huber loss approximates the hinge loss. Defaults to 0.5
  #'   \insertCite{Chapelle2007}{DPpack}.
  #'
  #' @return A new svmDP object.
  initialize = function(regularizer, eps, gamma,
                        perturbation.method = 'objective', kernel='linear',
                        D=NULL, kernel.param=NULL, regularizer.gr=NULL,
                        huber.h=0.5){
    super$initialize(mapXy.linear, generate.loss.huber(huber.h), regularizer, eps,
                     gamma, perturbation.method, 1/(2*huber.h), mapXy.gr.linear,
                     generate.loss.gr.huber(huber.h), regularizer.gr)
    self$kernel <- kernel
    if (kernel=="Gaussian"){
      self$phi <- phi.gaussian
      self$kernel.param <- kernel.param
      if (is.null(D)) stop("D must be specified for nonlinear kernel.")
      self$D <- D
    } else if (kernel!="linear") stop("kernel must be one of {'linear',
                                      Gaussian'}")
  },
  #' @description Fit the differentially private SVM model. This method runs
  #'   either the output perturbation or the objective perturbation algorithm
  #'   \insertCite{chaudhuri2011}{DPpack}, depending on the value of
  #'   perturbation.method used to construct the object, to generate an
  #'   objective function. A numerical optimization method is then run to find
  #'   optimal coefficients for fitting the model given the training data,
  #'   weights, and hyperparameters. The built-in \code{\link{optim}} function
  #'   using the "BFGS" optimization method is used. If \code{regularizer} is
  #'   given as 'l2' or if \code{regularizer.gr} is given in the construction of
  #'   the object, the gradient of the objective function is utilized by
  #'   \code{optim} as well. Otherwise, non-gradient based optimization methods
  #'   are used. The resulting privacy-preserving coefficients are stored in
  #'   \code{coeff}.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true labels for each row of \code{X}.
  #' @param upper.bounds Numeric vector of length \code{ncol(X)} giving upper
  #'   bounds on the values in each column of X. The \code{ncol(X)} values are
  #'   assumed to be in the same order as the corresponding columns of \code{X}.
  #'   Any value in the columns of \code{X} larger than the corresponding upper
  #'   bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length \code{ncol(X)} giving lower
  #'   bounds on the values in each column of \code{X}. The \code{ncol(X)}
  #'   values are assumed to be in the same order as the corresponding columns
  #'   of \code{X}. Any value in the columns of \code{X} larger than the
  #'   corresponding upper bound is clipped at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE.
  #' @param weights Numeric vector of observation weights of the same length as
  #'   \code{y}. If not given, no observation weighting is performed.
  #' @param weights.upper.bound Numeric value representing the global or public
  #'   upper bound on the weights. Required if weights are given.
  fit = function(X, y, upper.bounds, lower.bounds, add.bias=FALSE, weights=NULL,
                 weights.upper.bound=NULL){
    if (!is.null(weights) & self$perturbation.method!="output"){
      stop("Weighted SVM is only implemented for output perturbation.")
    }
    if (self$kernel=='Gaussian'){
      if (is.null(self$kernel.param)) self$kernel.param <- 1/ncol(X)
      self$sampling <- generate.sampling(self$kernel.param)
    }
    super$fit(X, y, upper.bounds, lower.bounds, add.bias, weights,
              weights.upper.bound)
  },
  #' @description Convert input data X into transformed data V. Uses sampled
  #'   pre-filter values and a mapping function based on the chosen kernel to
  #'   produce D-dimensional data V on which to train the model or predict
  #'   future values. This method is only used if the kernel is nonlinear. See
  #'   \insertCite{chaudhuri2011;textual}{DPpack} for more details.
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
  #' @description Predict label(s) for given \code{X} using the fitted
  #'   coefficients.
  #' @param X Dataframe of data on which to make predictions. Must be of same
  #'   form as \code{X} used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #' @param raw.value Boolean indicating whether to return the raw predicted
  #'   value or the rounded class label. If FALSE (default), outputs the
  #'   predicted labels 0 or 1. If TRUE, returns the raw score from the SVM
  #'   model.
  #'
  #' @return Matrix of predicted labels or scores corresponding to each row of
  #'   \code{X}.
  predict = function(X, add.bias=FALSE, raw.value=FALSE){
    if (self$kernel!="linear"){
      if (add.bias){
        X <- dplyr::mutate(X, bias =1)
        X <- X[, c(ncol(X), 1:(ncol(X)-1))]
      }
      V <- self$XtoV(as.matrix(X))
      if (raw.value) super$predict(V, FALSE)
      else {
        vals <- sign(super$predict(V, FALSE))
        vals[vals==-1] <- 0
        vals
      }
    } else{
      if (raw.value) super$predict(X, add.bias)
      else {
        vals <- sign(super$predict(X, add.bias))
        vals[vals==-1] <- 0
        vals
      }
    }
  }
), private=list(
  # description Preprocess input data and bounds to ensure they meet the
  #   assumptions necessary for fitting the model. If desired, a bias term can
  #   also be added.
  # param X Dataframe of data to be fit. Will be converted to a matrix.
  # param y Vector or matrix of true labels for each row of X. Will be
  #   converted to a matrix.
  # param upper.bounds Numeric vector of length ncol(X) giving upper bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param lower.bounds Numeric vector of length ncol(X) giving lower bounds on
  #   the values in each column of X. The ncol(X) values are assumed to be in
  #   the same order as the corresponding columns of X. Any value in the
  #   columns of X larger than the corresponding upper bound is clipped at the
  #   bound.
  # param add.bias Boolean indicating whether to add a bias term to X.
  # param weights Numeric vector of observation weights of the same length as
  #   \code{y}.
  # param weights.upper.bound Numeric value representing the global or public
  #   upper bound on the weights.
  # return A list of preprocessed values for X, y, upper.bounds, and
  #   lower.bounds for use in the privacy-preserving SVM algorithm.
  preprocess_data = function(X, y, upper.bounds, lower.bounds, add.bias,
                             weights=NULL, weights.upper.bound=NULL){
    if (self$kernel!="linear"){
      # Add bias if needed (highly recommended to not do this due to unwanted
      #       regularization of bias term)
      if (add.bias){
        X <- dplyr::mutate(X, bias=1)
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

      lbs <- c(numeric(ncol(V)) - sqrt(1/self$D))
      ubs <- c(numeric(ncol(V)) + sqrt(1/self$D))
      res <- super$preprocess_data(V, y, ubs, lbs, add.bias=FALSE,
                                   weights=weights,
                                   weights.upper.bound=weights.upper.bound)
    } else res <- super$preprocess_data(X, y, upper.bounds, lower.bounds,
                                        add.bias, weights, weights.upper.bound)

    X <- res$X
    y <- res$y
    upper.bounds <- res$upper.bounds
    lower.bounds <- res$lower.bounds

    # Process X
    p <- length(lower.bounds)
    tmp <- matrix(c(lower.bounds,upper.bounds),ncol=2)
    scale1 <- apply(abs(tmp),1,max)
    scale2 <- sqrt(p)
    X.norm <- t(t(X)/scale1)/scale2

    # Process y.
    # Labels must be 0 and 1
    if (any((y!=1) & (y!=0))) stop("y must be a numeric vector of binary labels
                                    0 and 1.")

    # Convert to 1 and -1 for SVM algorithm
    y[y==0] <- -1

    list(X=X.norm, y=y, scale1=scale1, scale2=scale2, weights=res$weights)
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

#'Privacy-preserving Empirical Risk Minimization for Regression
#'
#'@description This class implements differentially private empirical risk
#'  minimization using the objective perturbation technique
#'  \insertCite{Kifer2012}{DPpack}. It is intended to be a framework for
#'  building more specific models via inheritance. See
#'  \code{\link{LinearRegressionDP}} for an example of this type of structure.
#'
#'@details To use this class for empirical risk minimization, first use the
#'  \code{new} method to construct an object of this class with the desired
#'  function values and hyperparameters. After constructing the object, the
#'  \code{fit} method can be applied with a provided dataset and data bounds to
#'  fit the model. In fitting, the model stores a vector of coefficients
#'  \code{coeff} which satisfy differential privacy. These can be released
#'  directly, or used in conjunction with the \code{predict} method to privately
#'  predict the outcomes of new datapoints.
#'
#'  Note that in order to guarantee differential privacy for the empirical risk
#'  minimization model, certain constraints must be satisfied for the values
#'  used to construct the object, as well as for the data used to fit.
#'  Specifically, the following constraints must be met. Let \eqn{l} represent
#'  the loss function for an individual dataset row x and output value y and
#'  \eqn{L} represent the average loss over all rows and output values. First,
#'  \eqn{L} must be convex with a continuous Hessian. Second, the l2-norm of the
#'  gradient of \eqn{l} must be bounded above by some constant zeta for all
#'  possible input values in the domain. Third, for all possible inputs to
#'  \eqn{l}, the Hessian of \eqn{l} must be of rank at most one and its
#'  Eigenvalues must be bounded above by some constant lambda. Fourth, the
#'  regularizer must be convex. Finally, the provided domain of \eqn{l} must be
#'  a closed convex subset of the set of all real-valued vectors of dimension p,
#'  where p is the number of columns of X. Note that because of this, a bias
#'  term cannot be included without appropriate scaling/preprocessing of the
#'  dataset. To ensure privacy, the add.bias argument in the \code{fit} and
#'  \code{predict} methods should only be utilized in subclasses within this
#'  package where appropriate preprocessing is implemented, not in this class.
#'
#'@keywords internal
#'
#'@references \insertRef{Kifer2012}{DPpack}
#'
#' @examples
#' # Build example dataset
#' n <- 500
#' X <- data.frame(X=seq(-1,1,length.out = n))
#' true.theta <- c(-.3,.5) # First element is bias term
#' p <- length(true.theta)
#' y <- true.theta[1] + as.matrix(X)%*%true.theta[2:p] + stats::rnorm(n=n,sd=.1)
#'
#' # Construct object for linear regression
#' mapXy <- function(X, coeff) X%*%coeff
#' loss <- function(y.hat, y) (y.hat-y)^2/2
#' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
#' eps <- 1
#' delta <- 1
#' domain <- list("constraints"=function(coeff) coeff%*%coeff-length(coeff),
#'   "jacobian"=function(coeff) 2*coeff)
#' # Set p to be the number of predictors desired including intercept term (length of coeff)
#' zeta <- 2*p^(3/2) # Proper bound for linear regression
#' lambda <- p # Proper bound for linear regression
#' gamma <- 1
#' mapXy.gr <- function(X, coeff) t(X)
#' loss.gr <- function(y.hat, y) y.hat-y
#' regularizer.gr <- function(coeff) coeff
#'
#' ermdp <- EmpiricalRiskMinimizationDP.KST$new(mapXy, loss, 'l2', eps, delta,
#'                                              domain, zeta, lambda,
#'                                              gamma, mapXy.gr, loss.gr,
#'                                              regularizer.gr)
#'
#' # Fit with data
#' # We must assume y is a matrix with values between -p and p (-2 and 2
#' #   for this example)
#' upper.bounds <- c(1, 2) # Bounds for X and y
#' lower.bounds <- c(-1,-2) # Bounds for X and y
#' ermdp$fit(X, y, upper.bounds, lower.bounds, add.bias=TRUE)
#' ermdp$coeff # Gets private coefficients
#'
#' # Predict new data points
#' # Build a test dataset
#' Xtest <- data.frame(X=c(-.5, -.25, .1, .4))
#' predicted.y <- ermdp$predict(Xtest, add.bias=TRUE)
#'
#'@export
EmpiricalRiskMinimizationDP.KST <- R6::R6Class("EmpiricalRiskMinimizationDP.KST",
  public=list(
  #' @field mapXy Map function of the form \code{mapXy(X, coeff)} mapping input
  #'   data matrix \code{X} and coefficient vector or matrix \code{coeff} to
  #'   output labels \code{y}.
  mapXy = NULL,
  #' @field mapXy.gr Function representing the gradient of the map function with
  #'   respect to the values in \code{coeff} and of the form \code{mapXy.gr(X,
  #'   coeff)}, where \code{X} is a matrix and \code{coeff} is a matrix or
  #'   numeric vector.
  mapXy.gr = NULL,
  #' @field loss Loss function of the form \code{loss(y.hat, y)}, where
  #'   \code{y.hat} and \code{y} are matrices.
  loss = NULL,
  #' @field loss.gr Function representing the gradient of the loss function with
  #'   respect to \code{y.hat} and of the form \code{loss.gr(y.hat, y)}, where
  #'   \code{y.hat} and \code{y} are matrices.
  loss.gr = NULL,
  #' @field regularizer Regularization function of the form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix.
  regularizer = NULL,
  #' @field regularizer.gr Function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}.
  regularizer.gr = NULL,
  #' @field gamma Nonnegative real number representing the regularization
  #'   constant.
  gamma = NULL,
  #' @field eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  eps = NULL,
  #' @field delta Nonnegative real number defining the delta privacy parameter.
  #'   If 0, reduces to pure eps-DP.
  delta = NULL,
  #' @field domain List of constraint and jacobian functions representing the
  #'   constraints on the search space for the objective perturbation algorithm
  #'   used in \insertCite{Kifer2012;textual}{DPpack}.
  domain = NULL,
  #' @field zeta Positive real number denoting the upper bound on the l2-norm
  #'   value of the gradient of the loss function, as required to ensure
  #'   differential privacy.
  zeta = NULL,
  #' @field lambda Positive real number corresponding to the upper bound of the
  #'   Eigenvalues of the Hessian of the loss function for all possible inputs.
  lambda = NULL,
  #' @field coeff Numeric vector of coefficients for the model.
  coeff = NULL,
  #' @description Create a new \code{EmpiricalRiskMinimizationDP.KST} object.
  #' @param mapXy Map function of the form \code{mapXy(X, coeff)} mapping input
  #'   data matrix \code{X} and coefficient vector or matrix \code{coeff} to
  #'   output labels \code{y}. Should return a column matrix of predicted labels
  #'   for each row of \code{X}. See \code{\link{mapXy.linear}} for an example.
  #' @param loss Loss function of the form \code{loss(y.hat, y)}, where
  #'   \code{y.hat} and \code{y} are matrices. Should be defined such that it
  #'   returns a matrix of loss values for each element of \code{y.hat} and
  #'   \code{y}. See \code{\link{loss.squared.error}} for an example. This
  #'   function must be convex and the l2-norm of its gradient must be bounded
  #'   above by zeta for some constant zeta for all possible inputs within the
  #'   given domain. Additionally, for all possible inputs within the given
  #'   domain, the Hessian of the loss function must be of rank at most one and
  #'   its Eigenvalues must be bounded above by some constant lambda.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix, and
  #'   return the value of the regularizer at \code{coeff}. See
  #'   \code{\link{regularizer.l2}} for an example. Additionally, in order to
  #'   ensure differential privacy, the function must be convex.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param delta Nonnegative real number defining the delta privacy parameter.
  #'   If 0, reduces to pure eps-DP.
  #' @param domain List of functions representing the constraints on the search
  #'   space for the objective perturbation algorithm. Must contain two
  #'   functions, labeled "constraints" and "jacobian", respectively. The
  #'   "constraints" function accepts a vector of coefficients from the search
  #'   space and returns a value such that the value is nonpositive if and only
  #'   if the input coefficient vector is within the constrained search space.
  #'   The "jacobian" function also accepts a vector of coefficients and returns
  #'   the Jacobian of the constraint function. For example, in linear
  #'   regression, the square of the l2-norm of the coefficient vector
  #'   \eqn{\theta} is assumed to be bounded above by p, where p is the length
  #'   of \eqn{\theta} \insertCite{Kifer2012}{DPpack}. So, domain could be
  #'   defined as `domain <- list("constraints"=function(coeff)
  #'   coeff%*%coeff-length(coeff), "jacobian"=function(coeff) 2*coeff)`.
  #' @param zeta Positive real number denoting the upper bound on the l2-norm
  #'   value of the gradient of the loss function, as required to ensure
  #'   differential privacy.
  #' @param lambda Positive real number corresponding to the upper bound of the
  #'   Eigenvalues of the Hessian of the loss function for all possible inputs.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param mapXy.gr Optional function representing the gradient of the map
  #'   function with respect to the values in \code{coeff}. If given, must be of
  #'   the form \code{mapXy.gr(X, coeff)}, where \code{X} is a matrix and
  #'   \code{coeff} is a matrix or numeric vector. Should be defined such that
  #'   the ith row of the output represents the gradient with respect to the ith
  #'   coefficient. See \code{\link{mapXy.gr.linear}} for an example. If not
  #'   given, non-gradient based optimization methods are used to compute the
  #'   coefficient values in fitting the model.
  #' @param loss.gr Optional function representing the gradient of the loss
  #'   function with respect to \code{y.hat} and of the form
  #'   \code{loss.gr(y.hat, y)}, where \code{y.hat} and \code{y} are matrices.
  #'   Should be defined such that the ith row of the output represents the
  #'   gradient of the loss function at the ith set of input values. See
  #'   \code{\link{loss.gr.squared.error}} for an example. If not given,
  #'   non-gradient based optimization methods are used to compute the
  #'   coefficient values in fitting the model.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}. Should return a vector. See
  #'   \code{\link{regularizer.gr.l2}} for an example. If \code{regularizer} is
  #'   given as a string, this value is ignored. If not given and
  #'   \code{regularizer} is a function, non-gradient based optimization methods
  #'   are used to compute the coefficient values in fitting the model.
  #'
  #' @return A new EmpiricalRiskMinimizationDP.KST object.
  initialize = function(mapXy, loss, regularizer, eps, delta, domain, zeta, lambda,
                        gamma, mapXy.gr=NULL, loss.gr=NULL, regularizer.gr=NULL){
    self$mapXy <- mapXy # Must be of form mapXy(X, coeff)
    self$mapXy.gr <- mapXy.gr
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
  #'   numerical optimization method is then run to find optimal coefficients
  #'   for fitting the model given the training data and hyperparameters. The
  #'   \code{\link[nloptr]{nloptr}} function is used. If mapXy.gr, loss.gr, and
  #'   regularizer.gr are all given in the construction of the object, the
  #'   gradient of the objective function and the Jacobian of the constraint
  #'   function are utilized for the algorithm, and the NLOPT_LD_MMA method is
  #'   used. If one or more of these gradient functions are not given, the
  #'   NLOPT_LN_COBYLA method is used. The resulting privacy-preserving
  #'   coefficients are stored in coeff.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true values for each row of \code{X}.
  #' @param upper.bounds Numeric vector of length \code{ncol(X)+1} giving upper
  #'   bounds on the values in each column of \code{X} and the values of
  #'   \code{y}. The last value in the vector is assumed to be the upper bound
  #'   on \code{y}, while the first \code{ncol(X)} values are assumed to be in
  #'   the same order as the corresponding columns of \code{X}. Any value in the
  #'   columns of \code{X} and in \code{y} larger than the corresponding upper
  #'   bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length \code{ncol(X)+1} giving lower
  #'   bounds on the values in each column of \code{X} and the values of
  #'   \code{y}. The last value in the vector is assumed to be the lower bound
  #'   on \code{y}, while the first \code{ncol(X)} values are assumed to be in
  #'   the same order as the corresponding columns of \code{X}. Any value in the
  #'   columns of \code{X} and in \code{y} larger than the corresponding lower
  #'   bound is clipped at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE.
  fit = function(X, y, upper.bounds, lower.bounds, add.bias=FALSE){
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
    if (self$delta==0){
      norm.b <- rgamma(1, p, rate=self$eps/(2*self$zeta))
      direction.b <- stats::rnorm(p)
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
  #'   form as \code{X} used to fit coefficients.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE. If add.bias was set to TRUE when fitting the
  #'   coefficients, add.bias should be set to TRUE for predictions.
  #'
  #' @return Matrix of predicted y values corresponding to each row of X.
  predict = function(X, add.bias=FALSE){
    if (add.bias){
      X <- dplyr::mutate(X, bias =1)
      X <- X[, c(ncol(X), 1:(ncol(X)-1))]
    }
    self$mapXy(as.matrix(X), self$coeff)
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
    if (length(upper.bounds)!=(ncol(X)+1)){
      stop("Length of upper.bounds must be equal to the number of columns of X plus 1 (for y).");
    }
    if (length(lower.bounds)!=(ncol(X)+1)) {
      stop("Length of lower.bounds must be equal to the number of columns of X plus 1 (for y).");
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
  #   nloptr from the nloptr package. If mapXy.gr, loss.gr, and regularizer.gr are
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
      as.numeric(sum(self$loss(self$mapXy(X,par),y))/n +
                   self$gamma*self$regularizer(par)/n + Delta*par%*%par/(2*n) +
                   t(b)%*%par/n)
    }

    g <- function(par, X, y, Delta, b) self$domain$constraints(par) # For new opt method

    # Get gradient function
    if (!is.null(self$mapXy.gr) && !is.null(self$loss.gr) &&
        !is.null(self$regularizer.gr)) {
      objective.gr <- function(par, X, y, Delta, b){
        as.numeric(self$mapXy.gr(X,par)%*%self$loss.gr(self$mapXy(X,par),y)/n +
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
#' @details To use this class for linear regression, first use the \code{new}
#'   method to construct an object of this class with the desired function
#'   values and hyperparameters. After constructing the object, the \code{fit}
#'   method can be applied with a provided dataset and data bounds to fit the
#'   model. In fitting, the model stores a vector of coefficients \code{coeff}
#'   which satisfy differential privacy. These can be released directly, or used
#'   in conjunction with the \code{predict} method to privately predict the
#'   outcomes of new datapoints.
#'
#'   Note that in order to guarantee differential privacy for linear regression,
#'   certain constraints must be satisfied for the values used to construct the
#'   object, as well as for the data used to fit. The regularizer must be
#'   convex. Additionally, it is assumed that if x represents a single row of
#'   the dataset X, then the l2-norm of x is at most p for all x, where p is the
#'   number of predictors (including any possible intercept term). In order to
#'   ensure this constraint is satisfied, the dataset is preprocessed and
#'   scaled, and the resulting coefficients are postprocessed and un-scaled so
#'   that the stored coefficients correspond to the original data. Due to this
#'   constraint on x, it is best to avoid using an intercept term in the model
#'   whenever possible. If an intercept term must be used, the issue can be
#'   partially circumvented by adding a constant column to X before fitting the
#'   model, which will be scaled along with the rest of X. The \code{fit} method
#'   contains functionality to add a column of constant 1s to X before scaling,
#'   if desired.
#'
#' @references \insertRef{Kifer2012}{DPpack}
#'
#' @examples
#' # Build example dataset
#' n <- 500
#' X <- data.frame(X=seq(-1,1,length.out = n))
#' true.theta <- c(-.3,.5) # First element is bias term
#' p <- length(true.theta)
#' y <- true.theta[1] + as.matrix(X)%*%true.theta[2:p] + stats::rnorm(n=n,sd=.1)
#'
#' # Construct object for linear regression
#' regularizer <- 'l2' # Alternatively, function(coeff) coeff%*%coeff/2
#' eps <- 1
#' delta <- 0 # Indicates to use pure eps-DP
#' gamma <- 1
#' regularizer.gr <- function(coeff) coeff
#'
#' lrdp <- LinearRegressionDP$new('l2', eps, delta, gamma, regularizer.gr)
#'
#' # Fit with data
#' # We must assume y is a matrix with values between -p and p (-2 and 2
#' #   for this example)
#' upper.bounds <- c(1, 2) # Bounds for X and y
#' lower.bounds <- c(-1,-2) # Bounds for X and y
#' lrdp$fit(X, y, upper.bounds, lower.bounds, add.bias=TRUE)
#' lrdp$coeff # Gets private coefficients
#'
#' # Predict new data points
#' # Build a test dataset
#' Xtest <- data.frame(X=c(-.5, -.25, .1, .4))
#' predicted.y <- lrdp$predict(Xtest, add.bias=TRUE)
#'
#' @export
LinearRegressionDP <- R6::R6Class("LinearRegressionDP",
  inherit=EmpiricalRiskMinimizationDP.KST,
  public=list(
  #' @description Create a new LinearRegressionDP object.
  #' @param regularizer String or regularization function. If a string, must be
  #'   'l2', indicating to use l2 regularization. If a function, must have form
  #'   \code{regularizer(coeff)}, where \code{coeff} is a vector or matrix, and
  #'   return the value of the regularizer at \code{coeff}. See
  #'   \code{\link{regularizer.l2}} for an example. Additionally, in order to
  #'   ensure differential privacy, the function must be convex.
  #' @param eps Positive real number defining the epsilon privacy budget. If set
  #'   to Inf, runs algorithm without differential privacy.
  #' @param delta Nonnegative real number defining the delta privacy parameter.
  #'   If 0, reduces to pure eps-DP.
  #' @param gamma Nonnegative real number representing the regularization
  #'   constant.
  #' @param regularizer.gr Optional function representing the gradient of the
  #'   regularization function with respect to \code{coeff} and of the form
  #'   \code{regularizer.gr(coeff)}. Should return a vector. See
  #'   \code{\link{regularizer.gr.l2}} for an example. If \code{regularizer} is
  #'   given as a string, this value is ignored. If not given and
  #'   \code{regularizer} is a function, non-gradient based optimization methods
  #'   are used to compute the coefficient values in fitting the model.
  #'
  #' @return A new LinearRegressionDP object.
  initialize = function(regularizer, eps, delta, gamma, regularizer.gr=NULL){
      domain.linear <- list("constraints"=function(coeff) coeff%*%coeff-length(coeff),
                            "jacobian"=function(coeff) 2*coeff)
      super$initialize(mapXy.linear, loss.squared.error, regularizer, eps, delta,
                       domain.linear, NULL, NULL, gamma,
                       mapXy.gr=mapXy.gr.linear, loss.gr=loss.gr.squared.error,
                       regularizer.gr=regularizer.gr)
    },
  #' @description Fit the differentially private linear regression model. The
  #'   function runs the objective perturbation algorithm
  #'   \insertCite{Kifer2012}{DPpack} to generate an objective function. A
  #'   numerical optimization method is then run to find optimal coefficients
  #'   for fitting the model given the training data and hyperparameters. The
  #'   \code{\link[nloptr]{nloptr}} function is used. If \code{regularizer} is given as
  #'   'l2' or if \code{regularizer.gr} is given in the construction of the
  #'   object, the gradient of the objective function and the Jacobian of the
  #'   constraint function are utilized for the algorithm, and the NLOPT_LD_MMA
  #'   method is used. If this is not the case, the NLOPT_LN_COBYLA method is
  #'   used. The resulting privacy-preserving coefficients are stored in coeff.
  #' @param X Dataframe of data to be fit.
  #' @param y Vector or matrix of true values for each row of \code{X}.
  #' @param upper.bounds Numeric vector of length \code{ncol(X)+1} giving upper
  #'   bounds on the values in each column of \code{X} and the values of
  #'   \code{y}. The last value in the vector is assumed to be the upper bound
  #'   on \code{y}, while the first \code{ncol(X)} values are assumed to be in
  #'   the same order as the corresponding columns of \code{X}. Any value in the
  #'   columns of \code{X} and in \code{y} larger than the corresponding upper
  #'   bound is clipped at the bound.
  #' @param lower.bounds Numeric vector of length \code{ncol(X)+1} giving lower
  #'   bounds on the values in each column of \code{X} and the values of
  #'   \code{y}. The last value in the vector is assumed to be the lower bound
  #'   on \code{y}, while the first \code{ncol(X)} values are assumed to be in
  #'   the same order as the corresponding columns of \code{X}. Any value in the
  #'   columns of \code{X} and in \code{y} larger than the corresponding lower
  #'   bound is clipped at the bound.
  #' @param add.bias Boolean indicating whether to add a bias term to \code{X}.
  #'   Defaults to FALSE.
  fit = function(X, y, upper.bounds, lower.bounds, add.bias=FALSE){
    super$fit(X,y,upper.bounds,lower.bounds,add.bias)
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
    scale.X <- apply(abs(tmp),1,max)
    # Need each row at most norm sqrt(p), already achieved by first scaling
    X.norm <- t(t(X)/scale.X)

    # Process y
    lb.y <- lower.bounds[length(lower.bounds)]
    ub.y <- upper.bounds[length(upper.bounds)]
    if (add.bias) {
      shift.y <- lb.y + (ub.y - lb.y)/2 # Subtracting this centers at 0
      scale.y <- (ub.y-lb.y)/(2*p) # Dividing by this scales to max size p in abs()
    }
    else {
      # Cannot have shift term if no bias term to re-absorb it in post-processing
      shift.y <- 0
      # Dividing by this scales to max size p in abs()
      scale.y <- max(abs(c(lb.y, ub.y)))/p
    }

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

#' Privacy-preserving Hyperparameter Tuning for Linear Regression Models
#'
#' This function implements the privacy-preserving hyperparameter tuning
#' function for linear regression \insertCite{Kifer2012}{DPpack} using the
#' exponential mechanism. It accepts a list of models with various chosen
#' hyperparameters, a dataset X with corresponding values y, upper and lower
#' bounds on the columns of X and the values of y, and a boolean indicating
#' whether to add bias in the construction of each of the models. The data are
#' split into m+1 equal groups, where m is the number of models being compared.
#' One group is set aside as the validation group, and each of the other m
#' groups are used to train each of the given m models. The negative of the sum
#' of the squared error for each model on the validation set is used as the
#' utility values in the exponential mechanism
#' (\code{\link{ExponentialMechanism}}) to select a tuned model in a
#' privacy-preserving way.
#'
#' @param models Vector of linear regression model objects, each initialized
#'   with a different combination of hyperparameter values from the search space
#'   for tuning. Each model should be initialized with the same epsilon privacy
#'   parameter value eps. The tuned model satisfies eps-level differential
#'   privacy.
#' @param X Dataframe of data to be used in tuning the model. Note it is assumed
#'   the data rows and corresponding labels are randomly shuffled.
#' @param y Vector or matrix of true values for each row of X.
#' @param upper.bounds Numeric vector giving upper bounds on the values in each
#'   column of X and the values in y. Should be length ncol(X)+1. The first
#'   ncol(X) values are assumed to be in the same order as the corresponding
#'   columns of X, while the last value in the vector is assumed to be the upper
#'   bound on y. Any value in the columns of X and y larger than the
#'   corresponding upper bound is clipped at the bound.
#' @param lower.bounds Numeric vector giving lower bounds on the values in each
#'   column of X and the values in y. Should be length ncol(X)+1. The first
#'   ncol(X) values are assumed to be in the same order as the corresponding
#'   columns of X, while the last value in the vector is assumed to be the lower
#'   bound on y. Any value in the columns of X and y smaller than the
#'   corresponding lower bound is clipped at the bound.
#' @param add.bias Boolean indicating whether to add a bias term to X. Defaults
#'   to FALSE.
#' @return Single model object selected from the input list models with tuned
#'   parameters.
#' @examples
#' # Build example dataset
#' n <- 500
#' X <- data.frame(X=seq(-1,1,length.out = n))
#' true.theta <- c(-.3,.5) # First element is bias term
#' p <- length(true.theta)
#' y <- true.theta[1] + as.matrix(X)%*%true.theta[2:p] + stats::rnorm(n=n,sd=.1)
#'
#' # Grid of possible gamma values for tuning linear regression model
#' grid.search <- c(100, 1, .0001)
#'
#' # Construct objects for logistic regression parameter tuning
#' # Privacy budget should be the same for all models
#' eps <- 1
#' delta <- 0.01
#' linrdp1 <- LinearRegressionDP$new("l2", eps, delta, grid.search[1])
#' linrdp2 <- LinearRegressionDP$new("l2", eps, delta, grid.search[2])
#' linrdp3 <- LinearRegressionDP$new("l2", eps, delta, grid.search[3])
#' models <- c(linrdp1, linrdp2, linrdp3)
#'
#' # Tune using data and bounds for X and y based on their construction
#' upper.bounds <- c( 1, 2) # Bounds for X and y
#' lower.bounds <- c(-1,-2) # Bounds for X and y
#' tuned.model <- tune_linear_regression_model(models, X, y, upper.bounds,
#'                                             lower.bounds, add.bias=TRUE)
#' tuned.model$gamma # Gives resulting selected hyperparameter
#'
#' # tuned.model result can be used the same as a trained LogisticRegressionDP model
#' tuned.model$coeff # Gives coefficients for tuned model
#'
#' # Build a test dataset for prediction
#' Xtest <- data.frame(X=c(-.5, -.25, .1, .4))
#' predicted.y <- tuned.model$predict(Xtest, add.bias=TRUE)
#'
#' @references \insertRef{Kifer2012}{DPpack}
#'
#' @export
tune_linear_regression_model<- function(models, X, y, upper.bounds, lower.bounds,
                                     add.bias=FALSE){
  # Make sure values are correct
  m <- length(models)
  n <- length(y)
  # Split data into m+1 groups
  validateX <- X[seq(m+1,n,m+1),,drop=FALSE]
  validatey <- y[seq(m+1,n,m+1)]
  z <- numeric(m)
  for (i in 1:m){
    subX <- X[seq(i,n,m+1),,drop=FALSE]
    suby <- y[seq(i,n,m+1)]
    models[[i]]$fit(subX,suby,upper.bounds,lower.bounds,add.bias)
    validatey.hat <- models[[i]]$predict(validateX,add.bias)
    z[i] <- sum((validatey - validatey.hat)^2)
  }
  # Bounds for y give global sensitivity for exponential mechanism in linear
  #    regression case. Sensitivity can be computed as square of difference
  #    between upper and lower bounds.
  ub.y <- upper.bounds[length(upper.bounds)]
  lb.y <- lower.bounds[length(lower.bounds)]
  res <- ExponentialMechanism(-z, models[[1]]$eps, (ub.y-lb.y)^2,
                              candidates=models)
  res[[1]]
}
