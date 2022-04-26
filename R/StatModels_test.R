###############################################
### TESTING ERM.CMS VIA LOGISTIC REGRESSION ###
###############################################
### Build dataset
N <- 200
D <- 2
K <- 2
X <- data.frame()
y <- data.frame()

# May be helpful to set seed for testing
for (j in (1:K)){
  t <- seq(-.25, .25, length.out = N)
  if (j==1) m <- stats::rnorm(N,-.2, .1) # Soft margin
  if (j==2) m <- stats::rnorm(N, .2, .1)
  # if (j==1) m <- stats::rnorm(N,-.2, .05) # Hard margin
  # if (j==2) m <- stats::rnorm(N, .2, .05)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(1, 1)
lower.bounds <- c(-1,-1)

data <- cbind(X, y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

### To verify in unit circle
# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#   r = diameter / 2
#   tt <- seq(0,2*pi,length.out = npoints)
#   xx <- center[1] + r * cos(tt)
#   yy <- center[2] + r * sin(tt)
#   return(data.frame(x = xx, y = yy))
# }
# cir <- circleFun(c(0,0),2,npoints = 100)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                          size = 2, show.legend = F) +
  ggplot2::scale_colour_discrete(name = "Label") + ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::xlim(-.8, .8) + ggplot2::ylim(-.7, .7)

### No DP, no regularization
# gamma <- 0
# eps <- Inf

### No DP, regularization
# gamma <- .01
# eps <- Inf

### DP, no regularization
# gamma <- 0
# eps <- 1

### DP, regularization
gamma <- .01
eps <- 1

pm <- 'objective'
c <- 1/4 # Necessary constant for logistic regression if objective used
# pm <- 'output'
# c <- NULL
ermdp <- EmpiricalRiskMinimizationDP.CMS$new(mapXy.sigmoid, loss.cross.entropy,
                                             'l2', eps, gamma, pm, c,
                                             mapXy.gr.sigmoid,
                                             loss.gr.cross.entropy)

ermdp$fit(X, y, upper.bounds, lower.bounds)

theta <- ermdp$coeff

# Decision boundary grid
grid <- expand.grid(seq(-.8, .8, length.out = 100),
                    seq(-.7, .7, length.out = 100))
Z <- round(ermdp$predict(grid))
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# Decision boundary visualization
ggplot2::ggplot() + ggplot2::geom_point(data = data, ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                      size = 2, show.legend = F) +
  ggplot2::geom_tile(data = gridPred, ggplot2::aes(x = grid[, 1],y = grid[, 2],
                                 fill=as.character(Z)), alpha = 0.3,
            show.legend = F) + ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12)+ ggplot2::xlim(-.8,.8) + ggplot2::ylim(-.7,.7)
### END TESTING LOGISTIC REGRESSION USING ERM ###


#########################################
### TESTING LOGISTIC REGRESSION CLASS ###
#########################################
### Build dataset (no bias, satisfies constraints)
N <- 200
D <- 2
K <- 2
X <- data.frame()
y <- data.frame()

# May be helpful to set seed for testing
for (j in (1:K)){
  t <- seq(-.25, .25, length.out = N)
  if (j==1) m <- stats::rnorm(N,-.2, .1) # Soft margin
  if (j==2) m <- stats::rnorm(N, .2, .1)
  # if (j==1) m <- stats::rnorm(N,-.2, .05) # Hard margin
  # if (j==2) m <- stats::rnorm(N, .2, .05)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(1, 1)
lower.bounds <- c(-1,-1)

data <- cbind(X, y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

### To verify in unit circle
# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#   r = diameter / 2
#   tt <- seq(0,2*pi,length.out = npoints)
#   xx <- center[1] + r * cos(tt)
#   yy <- center[2] + r * sin(tt)
#   return(data.frame(x = xx, y = yy))
# }
# cir <- circleFun(c(0,0),2,npoints = 100)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                          size = 2, show.legend = F) +
  ggplot2::scale_colour_discrete(name = "Label") + ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::xlim(-.8, .8) + ggplot2::ylim(-.7, .7)
### End build dataset

### Build dataset (with bias, needs scaling)
# N <- 200 # number of points per class
# D <- 2 # dimensionality, we use 2D data for easy visulization
# K <- 2 # number of classes, binary for logistic regression
# X <- data.frame() # data matrix (each row = single example, can view as xy coordinates)
# y <- data.frame() # class labels
#
# # May be helpful to set seed for testing
#
# for (j in (1:K)){
#   # t, m are parameters of parametric equations x1, x2
#   t <- seq(0,1,length.out = N)
#   # add randomness
#   m <- stats::rnorm(N, j+0.5, 0.25)
#   Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
#   ytemp <- data.frame(matrix(j-1, N, 1))
#   X <- rbind(X, Xtemp)
#   y <- rbind(y, ytemp)
# }
#
# upper.bounds <- c(3,3)
# lower.bounds <- c(0,0)
#
# data <- cbind(X,y)
# y <- as.matrix(data[,3])
# colnames(data) <- c(colnames(X), 'label')
#
# ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)), size = 2) +
#   ggplot2::scale_colour_discrete(name  ="Label") +
#   ggplot2::ylim(0, 3) + ggplot2::coord_fixed(ratio = 1) +
#   ggplot2::ggtitle('Data to be classified') +
#   ggplot2::theme_bw(base_size = 12) +
#   ggplot2::theme(legend.position=c(0.85, 0.87))
### End build dataset

### No DP, no regularization
# gamma <- 0
# eps <- Inf

### No DP, regularization
# gamma <- .01
# eps <- Inf

### DP, no regularization
# gamma <- 0
# eps <- 1

### DP, regularization
gamma <- .01
eps <- 1

pm <- 'objective'
# pm <- 'output'
lrdp <- LogisticRegressionDP$new('l2', eps, gamma, pm)

### No bias, satisfies constraints
lrdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds)

### With bias, needs scaling
# lrdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds,add.bias=TRUE)

# lrdp$predict(X, add.bias=TRUE)

theta <- lrdp$coeff

# Decision boundary grid
### No bias, satisfies constraints
grid <- expand.grid(seq(-.8, .8, length.out = 100),
                    seq(-.7, .7, length.out = 100))
Z <- lrdp$predict(grid)
### With bias, needs scaling
# grid <- expand.grid(seq(0, 3, length.out = 100), seq(0, 3, length.out = 100))
# Z <- lrdp$predict(grid, add.bias=TRUE)

gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

### No bias, satisfies constraints
ggplot2::ggplot() + ggplot2::geom_point(data = data, ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                      size = 2, show.legend = F) +
  ggplot2::geom_tile(data = gridPred, ggplot2::aes(x = grid[, 1],y = grid[, 2],
                                 fill=as.character(Z)), alpha = 0.3,
            show.legend = F) + ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12)+ ggplot2::xlim(-.8,.8) + ggplot2::ylim(-.7,.7)

### With bias, needs scaling
# ggplot2::ggplot() +   ggplot2::geom_point(data = data,
#                         ggplot2::aes(x=x1, y=x2, color = as.character(label)),
#                         size = 2, show.legend = F) +
#   ggplot2::geom_tile(data = gridPred,
#             ggplot2::aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
#             alpha = 0.3, show.legend = F)+
#   ggplot2::ggtitle('Decision Boundary for Logistic Regression') +
#   ggplot2::coord_fixed(ratio = 1) +
#   ggplot2::theme_bw(base_size = 12)
### END TESTING LOGISTIC REGRESSION CLASS ###
#
#########################################
### TESTING SVM CLASS (LINEAR CASE) ###
#########################################
### Build dataset (no bias, satisfies constraints)
N <- 200
D <- 2
K <- 2
X <- data.frame()
y <- data.frame()

# May be helpful to set seed for testing
for (j in (1:K)){
  t <- seq(-.25, .25, length.out = N)
  if (j==1) m <- stats::rnorm(N,-.2, .1) # Soft margin
  if (j==2) m <- stats::rnorm(N, .2, .1)
  # if (j==1) m <- stats::rnorm(N,-.2, .05) # Hard margin
  # if (j==2) m <- stats::rnorm(N, .2, .05)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(1, 1)
lower.bounds <- c(-1,-1)

data <- cbind(X, y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

### To verify in unit circle
# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#   r = diameter / 2
#   tt <- seq(0,2*pi,length.out = npoints)
#   xx <- center[1] + r * cos(tt)
#   yy <- center[2] + r * sin(tt)
#   return(data.frame(x = xx, y = yy))
# }
# cir <- circleFun(c(0,0),2,npoints = 100)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                          size = 2, show.legend = F) +
  ggplot2::scale_colour_discrete(name = "Label") + ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::xlim(-.8, .8) + ggplot2::ylim(-.7, .7)
### End build dataset

### Build dataset (with bias, needs scaling)
# N <- 200 # number of points per class
# D <- 2 # dimensionality, we use 2D data for easy visulization
# K <- 2 # number of classes, binary for logistic regression
# X <- data.frame() # data matrix (each row = single example, can view as xy coordinates)
# y <- data.frame() # class labels
#
# # May be helpful to set seed for testing
#
# for (j in (1:K)){
#   # t, m are parameters of parametric equations x1, x2
#   t <- seq(0,1,length.out = N)
#   # add randomness
#   m <- stats::rnorm(N, j+0.5, 0.25)
#   Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
#   ytemp <- data.frame(matrix(j-1, N, 1))
#   X <- rbind(X, Xtemp)
#   y <- rbind(y, ytemp)
# }
#
# upper.bounds <- c(3,3)
# lower.bounds <- c(0,0)
#
# data <- cbind(X,y)
# y <- as.matrix(data[,3])
# colnames(data) <- c(colnames(X), 'label')
#
# ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)), size = 2) +
#   ggplot2::scale_colour_discrete(name  ="Label") +
#   ggplot2::ylim(0, 3) + ggplot2::coord_fixed(ratio = 1) +
#   ggplot2::ggtitle('Data to be classified') +
#   ggplot2::theme_bw(base_size = 12) +
#   ggplot2::theme(legend.position=c(0.85, 0.87))
### End build dataset

### No DP, no regularization
# gamma <- 0
# eps <- Inf

### No DP, regularization
# gamma <- .01
# eps <- Inf

### DP, no regularization
# gamma <- 0
# eps <- 1

### DP, regularization
gamma <- .01
eps <- 1

pm <- 'objective'
# pm <- 'output'
svmdp <- svmDP$new('l2', eps, gamma, pm, kernel='linear')

### No bias, satisfies constraints
svmdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds)

### With bias, needs scaling
# svmdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds,add.bias=TRUE)

# svmdp$predict(X, add.bias=TRUE)

theta <- svmdp$coeff

# Decision boundary grid
### No bias, satisfies constraints
grid <- expand.grid(seq(-.8, .8, length.out = 100),
                    seq(-.7, .7, length.out = 100))
Z <- svmdp$predict(grid)
### With bias, needs scaling
# grid <- expand.grid(seq(0, 3, length.out = 100), seq(0, 3, length.out = 100))
# Z <- svmdp$predict(grid, add.bias=TRUE)

gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

### No bias, satisfies constraints
ggplot2::ggplot() + ggplot2::geom_point(data = data, ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                      size = 2, show.legend = F) +
  ggplot2::geom_tile(data = gridPred, ggplot2::aes(x = grid[, 1],y = grid[, 2],
                                 fill=as.character(Z)), alpha = 0.3,
            show.legend = F) + ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12)+ ggplot2::xlim(-.8,.8) + ggplot2::ylim(-.7,.7)

### With bias, needs scaling
# ggplot2::ggplot() +   ggplot2::geom_point(data = data,
#                         ggplot2::aes(x=x1, y=x2, color = as.character(label)),
#                         size = 2, show.legend = F) +
#   ggplot2::geom_tile(data = gridPred,
#             ggplot2::aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
#             alpha = 0.3, show.legend = F)+
#   ggplot2::ggtitle('Decision Boundary for Logistic Regression') +
#   ggplot2::coord_fixed(ratio = 1) +
#   ggplot2::theme_bw(base_size = 12)
### END TESTING SVM CLASS ###

#######################################################
### TESTING GAUSSIAN KERNEL SVM AND COMPARING TO REGULAR SVM ### (Linear)
#######################################################
# Build dataset
N <- 200
D <- 2
K <- 2
X <- data.frame()
y <- data.frame()

# May be helpful to set seed for testing

for (j in (1:K)){
  t <- seq(-.25,.25,length.out = N)
  # add randomness
  if (j==1) m <- stats::rnorm(N,-.2,.1) # Soft margin
  if (j==2) m <- stats::rnorm(N, .2,.1)
  # if (j==1) m <- stats::rnorm(N,-.2, .05) # Hard margin
  # if (j==2) m <- stats::rnorm(N, .2, .05)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(1,1)
lower.bounds <- c(-1,-1)

data <- cbind(X,y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

# To verify in unit circle
# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#   r = diameter / 2
#   tt <- seq(0,2*pi,length.out = npoints)
#   xx <- center[1] + r * cos(tt)
#   yy <- center[2] + r * sin(tt)
#   return(data.frame(x = xx, y = yy))
# }
# cir <- circleFun(c(0,0),2,npoints = 100)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                          size = 2, show.legend = F) +
  ggplot2::scale_colour_discrete(name  ="Label") + #ggplot2::ylim(-3, 3) +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::xlim(-.8,.8) + ggplot2::ylim(-.7,.7)

# Hard-margin
# eps <- Inf
# # gamma <- 1
# gamma <- 0
# D <- 20

# Soft-margin
# eps <- Inf
# gamma <- .01
# D <- 20

# With DP
eps <- 1
gamma <- .01
D <- 5

pm <- 'objective'
# pm <- 'output'
ksvmdp <- svmDP$new('l2', eps, gamma, pm, kernel='Gaussian', D)

ksvmdp$fit(X,y) # Note bounds unnecessary for Gaussian kernel

theta <- ksvmdp$coeff

# Grid
grid <- expand.grid(seq(-.8, .8, length.out = 100), seq(-.7,.7, length.out = 100))
Z <- ksvmdp$predict(grid)
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# decision boundary visualization
ggplot2::ggplot() +   ggplot2::geom_point(data = data,
                        ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                        size = 2, show.legend = F) +
  ggplot2::geom_tile(data = gridPred,
            ggplot2::aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
            alpha = 0.3, show.legend = F)+
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::xlim(-.8,.8) + ggplot2::ylim(-.7,.7)

# Standard SVM
stdSVM <- e1071::svm(as.matrix(X),y,kernel="radial")

# Grid
grid <- expand.grid(seq(-.8, .8, length.out = 100), seq(-.7,.7, length.out = 100))
Z <- sign(predict(stdSVM, newdata=grid,decision.values=TRUE))
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# decision boundary visualization
ggplot2::ggplot() +   ggplot2::geom_point(data = data,
                        ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                        size = 2, show.legend = F) +
  ggplot2::geom_tile(data = gridPred,
            ggplot2::aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
            alpha = 0.3, show.legend = F)+
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) + ggplot2::xlim(-.8,.8) + ggplot2::ylim(-.7,.7)

### END TESTING AND COMPARING KERNEL SVM ###

#######################################################
### TESTING GAUSSIAN KERNEL SVM ON NONLINEAR DATASET ### (Nonlinear)
#######################################################
# Build dataset
N <- 400 # number of points
D <- 2 # dimensionality, we use 2D data for easy visulization
K <- 2 # number of classes, binary for logistic regression
X <- data.frame() # data matrix (each row = single example, can view as xy coordinates)
y <- data.frame() # class labels

# May be helpful to set seed for testing

for (i in (1:N)){
  Xtemp <- data.frame(x1 = stats::rnorm(1,sd=.28) , x2 = stats::rnorm(1,sd=.28))
  if (sum(Xtemp^2)<.15) ytemp <- data.frame(y=0)
  else ytemp <- data.frame(y=1)
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(1,1)
lower.bounds <- c(-1,-1)

data <- cbind(X,y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

# To verify in unit circle
# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#   r = diameter / 2
#   tt <- seq(0,2*pi,length.out = npoints)
#   xx <- center[1] + r * cos(tt)
#   yy <- center[2] + r * sin(tt)
#   return(data.frame(x = xx, y = yy))
# }
# cir <- circleFun(c(0,0),2,npoints = 100)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)), size = 2,
                          show.legend = F) +
  ggplot2::ylim(-1, 1) + ggplot2::xlim(-1,1)+
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12)

eps <- 5
gamma <- 0.01
D <- 20

# pm <- 'objective'
pm <- 'output'
ksvmdp <- svmDP$new('l2', eps, gamma, pm, kernel='Gaussian', D)

# Note bounds unnecessary with Gaussian kernel
ksvmdp$fit(X,y)

theta <- ksvmdp$coeff

# Grid
grid <- expand.grid(seq(-1, 1, length.out = 100), seq(-1, 1, length.out = 100))
Z <- ksvmdp$predict(grid)
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# decision boundary visualization
ggplot2::ggplot() +   ggplot2::geom_point(data = data,
                        ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                        size = 2, show.legend = F) + ggplot2::ylim(-1, 1) + ggplot2::xlim(-1,1)+
  ggplot2::geom_tile(data = gridPred,
            ggplot2::aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
            alpha = 0.3, show.legend = F)+
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12)
### END TESTING KERNEL SVM ON NONLINEAR DATASET ###


############################################################
### TESTING PARAMETER TUNING FUNCTION FOR CLASSIFICATION ###
############################################################
# Build dataset
N <- 200 # number of points per class
D <- 2 # dimensionality, we use 2D data for easy visulization
K <- 2 # number of classes, binary for logistic regression
X <- data.frame() # data matrix (each row = single example, can view as xy coordinates)
y <- data.frame() # class labels

# May be helpful to set seed for testing

for (j in (1:K)){
  t <- seq(-.25,.25,length.out = N)
  if (j==1) m <- stats::rnorm(N,-.2,.1)
  if (j==2) m <- stats::rnorm(N, .2,.1)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(1,1)
lower.bounds <- c(-1,-1)

data <- cbind(X,y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

# To verify in unit circle
# circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
#   r = diameter / 2
#   tt <- seq(0,2*pi,length.out = npoints)
#   xx <- center[1] + r * cos(tt)
#   yy <- center[2] + r * sin(tt)
#   return(data.frame(x = xx, y = yy))
# }
# cir <- circleFun(c(0,0),2,npoints = 100)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=x1, y=x2, color = as.character(label)), size = 2) +
  ggplot2::scale_colour_discrete(name  ="Label") +
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position=c(0.85, 0.87))

eps <- 1

lrdp1 <- LogisticRegressionDP$new("l2", eps, 100)
lrdp2 <- LogisticRegressionDP$new("l2", eps, 1)
lrdp3 <- LogisticRegressionDP$new("l2", eps, .0001)

models <- c(lrdp1, lrdp2, lrdp3)
model <- tune_classification_model(models, X, y, upper.bounds, lower.bounds)
theta <- model$coeff

# Grid
grid <- expand.grid(seq(-.8, .8, length.out = 100), seq(-.7,.7, length.out = 100))
Z <- model$predict(grid)
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# decision boundary
ggplot2::ggplot() +   ggplot2::geom_point(data = data,
                        ggplot2::aes(x=x1, y=x2, color = as.character(label)),
                        size = 2, show.legend = F) +
  ggplot2::geom_tile(data = gridPred,
            ggplot2::aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
            alpha = 0.3, show.legend = F)+
  ggplot2::coord_fixed(ratio = 1) +
  ggplot2::theme_bw(base_size = 12)
model$gamma

### END TESTING PARAMETER TUNING FUNCTION ###


###############################################
### TESTING ERMDP.KST VIA LINEAR REGRESSION ###
###############################################
# Build dataset
n <- 500
p <- 2 # for easy visualization
X <- data.frame(X=seq(-1,1,length.out = n))
true.theta <- c(-.3,.5) # Normal case (within search space bounds - shows variation)
# true.theta <- c(1,2) # Too large coeff case (outside search space bounds - no variation)
y <- true.theta[1] + as.matrix(X)%*%true.theta[2:length(true.theta)] + stats::rnorm(n=n,sd=.1)
y[y< -p] <- -p
y[y>p] <- p

data <- cbind(X,y)
y <- as.matrix(data[,2])
colnames(data) <- c(colnames(X), 'y')
# Includes everything
ub <- c(1,2)
lb <- c(-1,-2)
# Limits some values
# ub <- c(1, .75)
# lb <- c(-1, -.75)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=X, y=y), size = 2) +
  ggplot2::ylim(min(y), max(y)) + #ggplot2::coord_fixed(ratio = 1) +
  ggplot2::ggtitle('Data to be regressed') +
  ggplot2::theme_bw(base_size = 12)

eps <- 1
delta <- .01
domain <- list("constraints"=function(coeff) coeff%*%coeff - length(coeff),
               "jacobian"=function(coeff) 2*coeff)
zeta <- 2*p^(3/2)
lambda <- p
gamma <- 1

ermdp.kst <- EmpiricalRiskMinimizationDP.KST$new(mapXy.linear, loss.squared.error,
                                                 'l2', eps, delta, domain, zeta,
                                                 lambda, gamma, mapXy.gr.linear,
                                                 loss.gr.squared.error)
ermdp.kst$fit(X,y,ub,lb,add.bias=TRUE)

theta <- ermdp.kst$coeff

grid <- seq(-1,1,length.out=100)

gridData <- data.frame(X=grid,y=theta[1]+theta[2]*grid)

# Regression line visualization
ggplot2::ggplot() + ggplot2::geom_point(data = data,
                      ggplot2::aes(x=X, y=y),
                      size = 2, show.legend = F) +
  ggplot2::geom_line(data=gridData, ggplot2::aes(x=X,y=y), color='red') + ggplot2::ylim(min(y), max(y)) +
  ggplot2::ggtitle('Regression Line for Linear Regression') +
  ggplot2::theme_bw(base_size = 12)
theta

### END TESTING ERMDP.KST WITH LINEAR REGRESSION ###


#################################
### TESTING LINEAR REGRESSION ###
#################################
# Build dataset
n <- 500
X <- data.frame(X=seq(-1,1,length.out = n))
true.theta <- c(-.3,.5) # Normal case
# true.theta <- c(1,2) # Too large coeff case
y <- true.theta[1] + as.matrix(X)%*%true.theta[2:length(true.theta)] + stats::rnorm(n=n,sd=.1)
p <- length(true.theta)
y[y< -p] <- -p
y[y> p] <- p

data <- cbind(X,y)
y <- as.matrix(data[,2])
colnames(data) <- c(colnames(X), 'y')
ub <- c(1, p)
lb <- c(-1, -p)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=X, y=y), size = 2) +
  ggplot2::ylim(min(y), max(y)) + #ggplot2::coord_fixed(ratio = 1) +
  ggplot2::ggtitle('Data to be regressed') +
  ggplot2::theme_bw(base_size = 12)

regularizer <- function(coeff) coeff%*%coeff/2
regularizer.gr <- function(coeff) coeff
eps <- 1
delta <- 0.01
gamma <- 1

linrdp <- LinearRegressionDP$new(regularizer, eps, delta, gamma, regularizer.gr)
linrdp$fit(X,y,ub,lb,add.bias=TRUE)

theta <- linrdp$coeff

grid <- seq(-1,1,length.out=100)

gridData <- data.frame(X=grid,y=theta[1]+theta[2]*grid)

# Regression line visualization
ggplot2::ggplot() + ggplot2::geom_point(data = data,
                      ggplot2::aes(x=X, y=y),
                      size = 2, show.legend = F) +
  ggplot2::geom_line(data=gridData, ggplot2::aes(x=X,y=y), color='red') + ggplot2::ylim(min(y), max(y)) +
  ggplot2::ggtitle('Regression Line for Linear Regression') +
  ggplot2::theme_bw(base_size = 12)
theta

### END TESTING LINEAR REGRESSION ###

###############################################################
### TESTING PARAMETER TUNING FUNCTION FOR LINEAR REGRESSION ###
###############################################################
# Build dataset
n <- 500
X <- data.frame(X=seq(-1,1,length.out = n))
true.theta <- c(-.3,.5) # Normal case
# true.theta <- c(1,2) # Too large coeff case
y <- true.theta[1] + as.matrix(X)%*%true.theta[2:length(true.theta)] + stats::rnorm(n=n,sd=.1)
p <- length(true.theta)
y[y< -p] <- -p
y[y> p] <- p

data <- cbind(X,y)
y <- as.matrix(data[,2])
colnames(data) <- c(colnames(X), 'y')
ub <- c(1, p)
lb <- c(-1, -p)

ggplot2::ggplot(data) + ggplot2::geom_point(ggplot2::aes(x=X, y=y), size = 2) +
  ggplot2::ylim(min(y), max(y)) + #ggplot2::coord_fixed(ratio = 1) +
  ggplot2::ggtitle('Data to be regressed') +
  ggplot2::theme_bw(base_size = 12)

eps <- 1
delta <- 0.01

linrdp1 <- LinearRegressionDP$new("l2", eps, delta, 100)
linrdp2 <- LinearRegressionDP$new("l2", eps, delta, 1)
linrdp3 <- LinearRegressionDP$new("l2", eps, delta, .0001)

models <- c(linrdp1, linrdp2, linrdp3)
model <- tune_linear_regression_model(models, X, y, ub, lb, add.bias=TRUE)
theta <- model$coeff

grid <- seq(-1,1,length.out=100)

gridData <- data.frame(X=grid,y=theta[1]+theta[2]*grid)

# Regression line visualization
ggplot2::ggplot() + ggplot2::geom_point(data = data,
                      ggplot2::aes(x=X, y=y),
                      size = 2, show.legend = F) +
  ggplot2::geom_line(data=gridData, ggplot2::aes(x=X,y=y), color='red') + ggplot2::ylim(min(y), max(y)) +
  ggplot2::ggtitle('Regression Line for Linear Regression') +
  ggplot2::theme_bw(base_size = 12)
model$gamma

### END TESTING PARAMETER TUNING FUNCTION FOR LINEAR REGRESSION ###
