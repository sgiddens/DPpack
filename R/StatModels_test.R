library(ggplot2)

#############################################
### TESTING LOGISTIC REGRESSION USING ERM ###
#############################################
sigmoid <- function(z){1/(1+exp(-z))}

h <- function(X, coeff){
  sigmoid(X%*%coeff)
}

h.gr <- function(X, coeff){
  as.numeric(sigmoid(X%*%coeff)*(1-sigmoid(X%*%coeff)))*t(X)
}

loss <- function(y.hat,y){
  -(y*log(y.hat) + (1-y)*log(1-y.hat))
}

loss.gr <- function(y.hat,y){
  -y/y.hat + (1-y)/(1-y.hat)
}

regularizer <- function(coeff){
  coeff%*%coeff/2
}

regularizer.gr <- function(coeff){
  coeff
}

lambda <- 0
# eps <- 2
eps <- Inf
c <- 1/4 # From paper

N <- 200 # number of points per class
D <- 2 # dimensionality, we use 2D data for easy visulization
K <- 2 # number of classes, binary for logistic regression
X <- data.frame() # data matrix (each row = single example, can view as xy coordinates)
y <- data.frame() # class labels

set.seed(56)

for (j in (1:K)){
  # t, m are parameters of parametric equations x1, x2
  t <- seq(0,1,length.out = N)
  # add randomness
  m <- rnorm(N, j+0.5, 0.25)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(3,3,1)
lower.bounds <- c(0,0,0)

X <- dplyr::mutate(X, bias =1)
X <- X[, c(ncol(X), 1:(ncol(X)-1))]
upper.bounds <- c(1,upper.bounds)
lower.bounds <- c(1,lower.bounds)
lb <- lower.bounds[1:(length(lower.bounds)-1)]
ub <- upper.bounds[1:(length(upper.bounds)-1)]
tmp <- matrix(c(lb,ub),ncol=2)
scale1 <- apply(abs(tmp),1,max)
scale2 <- sqrt(sum(scale1^2))
X.norm <- t(t(X)/scale1)/scale2

data <- cbind(X,y)
y <- as.matrix(data[,4])
colnames(data) <- c(colnames(X), 'label')

ggplot(data) + geom_point(aes(x=x1, y=x2, color = as.character(label)), size = 2) +
  scale_colour_discrete(name  ="Label") +
  ylim(0, 3) + coord_fixed(ratio = 1) +
  ggtitle('Data to be classified') +
  theme_bw(base_size = 12) +
  theme(legend.position=c(0.85, 0.87))

EmpiricalRiskMinimizationDP.CMS$undebug("fit")
EmpiricalRiskMinimizationDP.CMS$debug("optimize_coeff")

set.seed(round(runif(1,0,100)))

ermdp <- EmpiricalRiskMinimizationDP.CMS$new(h, loss, 'l2', eps, lambda, c, h.gr,
                                loss.gr)

ermdp$fit(X.norm,y)#,upper.bounds=upper.bounds,lower.bounds=lower.bounds)

ermdp$coeff <- ermdp$coeff/(scale1*scale2)
# round(ermdp$predict(X,add.bias=TRUE))

theta <- ermdp$coeff

# generate a grid for decision boundary, this is the test set
grid <- expand.grid(seq(0, 3, length.out = 100), seq(0, 3, length.out = 100))
# grid <- mutate(grid, bias =1)
# grid <- as.matrix(grid[, c(ncol(grid), 1:(ncol(grid)-1))])
# predict the probability
# probZ <- logisticProb(theta, grid)
# predict the label
# Z <- logisticPred(probZ)
Z <- round(ermdp$predict(grid,add.bias=TRUE))
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# decision boundary visualization
ggplot() +   geom_point(data = data,
                        aes(x=x1, y=x2, color = as.character(label)),
                        size = 2, show.legend = F) +
  geom_tile(data = gridPred,
            aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
            alpha = 0.3, show.legend = F)+
  ggtitle('Decision Boundary for Logistic Regression') +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 12) #+geom_line(data=my.data,aes(x=my.x,y=my.y))+
# geom_line(data=my.data2,aes(x=my.x2,y=my.y2)) +
# geom_point(data=data.frame(t(colMeans(data[,1:2]))), aes(x=x1,y=x2))

error <- sum(abs(round(ermdp$predict(X,add.bias=TRUE))-y))/length(y)
### END TESTING LOGISTIC REGRESSION USING ERM ###


#########################################
### TESTING LOGISTIC REGRESSION CLASS ###
#########################################
regularizer <- function(coeff){
  coeff%*%coeff/2
}

regularizer.gr <- function(coeff){
  coeff
}

lambda <- 0.1
# lambda <- 0
eps <- 2
# eps <- Inf

N <- 200 # number of points per class
D <- 2 # dimensionality, we use 2D data for easy visulization
K <- 2 # number of classes, binary for logistic regression
X <- data.frame() # data matrix (each row = single example, can view as xy coordinates)
y <- data.frame() # class labels

set.seed(56)

for (j in (1:K)){
  # t, m are parameters of parametric equations x1, x2
  t <- seq(0,1,length.out = N)
  # add randomness
  m <- rnorm(N, j+0.5, 0.25)
  Xtemp <- data.frame(x1 = 3*t , x2 = m - t)
  ytemp <- data.frame(matrix(j-1, N, 1))
  X <- rbind(X, Xtemp)
  y <- rbind(y, ytemp)
}

upper.bounds <- c(3,3,1)
lower.bounds <- c(0,0,0)

data <- cbind(X,y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'label')

ggplot(data) + geom_point(aes(x=x1, y=x2, color = as.character(label)), size = 2) +
  scale_colour_discrete(name  ="Label") +
  ylim(0, 3) + coord_fixed(ratio = 1) +
  ggtitle('Data to be classified') +
  theme_bw(base_size = 12) +
  theme(legend.position=c(0.85, 0.87))

EmpiricalRiskMinimizationDP.CMS$undebug("optimize_coeff")

set.seed(round(runif(1,0,100)))

lrdp <- LogisticRegressionDP$new('l2', eps, lambda)

lrdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds,add.bias=TRUE)

# lrdp$predict(X, add.bias=TRUE)

theta <- lrdp$coeff

# generate a grid for decision boundary, this is the test set
grid <- expand.grid(seq(0, 3, length.out = 100), seq(0, 3, length.out = 100))
# grid <- mutate(grid, bias =1)
# grid <- as.matrix(grid[, c(ncol(grid), 1:(ncol(grid)-1))])
# predict the probability
# probZ <- logisticProb(theta, grid)
# predict the label
# Z <- logisticPred(probZ)
Z <- lrdp$predict(grid, add.bias=TRUE)
gridPred = cbind(grid, Z)
colnames(gridPred)[3] <- 'label'
gridPred <- data.frame(gridPred)

# decision boundary visualization
ggplot() +   geom_point(data = data,
                        aes(x=x1, y=x2, color = as.character(label)),
                        size = 2, show.legend = F) +
  geom_tile(data = gridPred,
            aes(x = grid[, 1],y = grid[, 2], fill=as.character(Z)),
            alpha = 0.3, show.legend = F)+
  ggtitle('Decision Boundary for Logistic Regression') +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 12) #+geom_line(data=my.data,aes(x=my.x,y=my.y))+
# geom_line(data=my.data2,aes(x=my.x2,y=my.y2)) +
# geom_point(data=data.frame(t(colMeans(data[,1:2]))), aes(x=x1,y=x2))

error <- sum(abs(round(lrdp$predict(X,add.bias=TRUE))-y))/length(y)
### END TESTING LOGISTIC REGRESSION CLASS ###

#######################################
### TESTING LINEAR REGRESSION CLASS ###
#######################################

## Test with visualization
regularizer <- function(coeff){
  coeff%*%coeff/2
}

lambda <- 0.0005
eps <- 1
# eps <- Inf
c <- 2 # From paper

N <- 400 # number of points
D <- 1 # dimensionality, we use 1D data for easy visulization
set.seed(56)
X <- data.frame(x1=runif(N,0,3)) # data matrix (each row = single example)
true.beta <- c(5,2)
y <- true.beta[1]+true.beta[2]*X + rnorm(N,sd=.5)
colnames(y) <- c("y")

data <- cbind(X,y)
y <- as.matrix(data[,2])
colnames(data) <- c(colnames(X), 'y')
upper.bounds <- c(3)
lower.bounds <- c(0)

ggplot(data) + geom_point(aes(x=x1, y=y), size = 2) +
  ylim(2, 13) + #coord_fixed(ratio = 1) +
  ggtitle('Data to be regressed') +
  theme_bw(base_size = 12)

set.seed(round(runif(1,0,100)))

linrdp <- LinearRegressionDP$new(regularizer, lambda, eps, c)

linrdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds)

theta <- linrdp$coeff

grid <- seq(0,3,length.out=100)

gridData <- data.frame(x1=grid,y=theta[1]+theta[2]*grid)

# Regression line visualization
ggplot() + geom_point(data = data,
                        aes(x=x1, y=y),
                        size = 2, show.legend = F) +
  geom_line(data=gridData, aes(x=x1,y=y), color='red') +
  ggtitle('Regression Line for Linear Regression') +
  theme_bw(base_size = 12)

## Test with more features
regularizer <- function(coeff){
  coeff%*%coeff/2
}

lambda <- 0.0005
# lambda <- 0
eps <- 8
# eps <- Inf
c <- 2 # TODO: Verify this

N <- 400 # number of points
D <- 5 # dimensionality
set.seed(56)
X <- data.frame(tmp=numeric(N))
upper.bounds <- numeric()
lower.bounds <- numeric()
for (i in 1:D){
  upper.bounds[i] <- 3*i-4
  lower.bounds[i] <- 3*i-6
  X <- cbind(X, runif(N,lower.bounds[i],upper.bounds[i]))
}
X<-X[,2:(D+1)]
colnames(X) <- paste("x",as.character(1:D),sep="")
true.beta <- c(5,1,-10,3,2,8)

y <- true.beta[1]+as.matrix(X)%*%true.beta[2:(D+1)] + rnorm(N,sd=.5)
colnames(y) <- c("y")

data <- cbind(X,y)
y <- as.matrix(data[,(D+1)])
colnames(data) <- c(colnames(X), 'y')

set.seed(round(runif(1,0,100)))

linrdp <- LinearRegressionDP$new(regularizer, lambda, eps, c)

linrdp$fit(X,y,upper.bounds=upper.bounds,lower.bounds=lower.bounds)

theta <- linrdp$coeff
abs(theta-true.beta)

### END TESTING LINEAR REGRESSION CLASS ###

################################################
### TESTING LINEAR REGRESSION REGULARIZATION ###
################################################
data <- data.frame(x1=seq(30,50,length.out=50),x2=seq(1,5,length.out=50))
data$y <- data$x1+data$x2 + rnorm(50)
X <- as.matrix(data[,1:2])
y <- as.matrix(data[3])
n <- 500
beta.hat.seq <- matrix(nrow=n,ncol=2)
pts <- seq(0,5000,length.out=n)
for (i in 1:n){
  beta.hat.seq[i,] <- solve((t(X)%*%X + pts[i]*diag(2)),t(X)%*%y)
}
beta.hat.seq <- data.frame(beta.hat.seq)
colnames(beta.hat.seq) <- c("beta1", "beta2")

loss <- function(X, y, beta, lambda) sum((y-X%*%beta)^2)/(2*length(y)) + lambda*beta%*%beta

n.circle <- 51
circle <- matrix(nrow=n.circle, ncol=3)
circle.pts <- seq(-pi/8,pi/8,length.out=n.circle)
for (t in 1:n.circle){
  circle[t,] <- c(cos(circle.pts[t]),sin(circle.pts[t]),
                  loss(X,y,c(cos(circle.pts[t]),sin(circle.pts[t])),500))
}
circle <- data.frame(circle)
colnames(circle) <- c("beta1","beta2","loss")

ggplot() + geom_point(aes(beta1,beta2), data=beta.hat.seq) + xlim(.5,1.5) +
  ylim(-.25,.5) + geom_point(aes(beta1,beta2,color=loss),data=circle) +
  geom_point(aes(beta1,beta2),data=circle[which.min(circle$loss),],col='green')

### END TESTING LINEAR REGRESSION REGULARIZATION ###


################################################
### TESTING ERMDP.KST WITH LINEAR REGRESSION ###
################################################
# Build dataset
n <- 500
p <- 2 # for easy visualization
X <- data.frame(const=numeric(n)+1, X=seq(-1,1,length.out = n))
# true.theta <- c(-.3,.5) # Normal case (shows variation)
true.theta <- c(1,2) # Too large coeff case (no variation)
y <- as.matrix(X)%*%true.theta + rnorm(n=n,sd=.1)
y[y< -p] <- -p
y[y>p] <- p

data <- cbind(X,y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'y')
# Includes everything
ub <- c(1, 1, 2)
lb <- c(1, -1, -2)
# Limits some values
# ub <- c(1, 1, .75)
# lb <- c(1, -1, -.75)

ggplot(data) + geom_point(aes(x=X, y=y), size = 2) +
  ylim(min(y), max(y)) + #coord_fixed(ratio = 1) +
  ggtitle('Data to be regressed') +
  theme_bw(base_size = 12)

h <- function(X, coeff) X%*%coeff
h.gr <- function(X, coeff) t(X)
loss <- function(y.hat, y) (y.hat-y)^2/2
loss.gr <- function(y.hat, y) y.hat-y
regularizer <- function(coeff) coeff%*%coeff/2
regularizer.gr <- function(coeff) coeff
eps <- 1
delt <- 1
domain <- list("constraints"=function(coeff) coeff%*%coeff - length(coeff),
               "jacobian"=function(coeff) 2*coeff)
zeta <- 2*p^(3/2)
lambda <- p
gamma <- 1

ermdp.kst <- EmpiricalRiskMinimizationDP.KST$new(h, loss, regularizer, eps,
                                                 delt, domain, zeta, lambda,
                                                 gamma, h.gr, loss.gr,
                                                 regularizer.gr)
EmpiricalRiskMinimizationDP.KST$undebug("fit")
ermdp.kst$fit(X,y,ub,lb)

theta <- ermdp.kst$coeff

grid <- seq(-1,1,length.out=100)

gridData <- data.frame(X=grid,y=theta[1]+theta[2]*grid)

# Regression line visualization
ggplot() + geom_point(data = data,
                      aes(x=X, y=y),
                      size = 2, show.legend = F) +
  geom_line(data=gridData, aes(x=X,y=y), color='red') + ylim(min(y), max(y)) +
  ggtitle('Regression Line for Linear Regression') +
  theme_bw(base_size = 12)
theta

### END TESTING ERMDP.KST WITH LINEAR REGRESSION ###


#################################
### TESTING LINEAR REGRESSION ###
#################################
# Build dataset
n <- 500
p <- 2 # for easy visualization
X <- data.frame(const=numeric(n)+1, X=seq(-1,.5,length.out = n))
# true.theta <- c(-.3,.5) # Normal case (shows variation)
true.theta <- c(1,2) # Too large coeff case (no variation)
y <- as.matrix(X)%*%true.theta + rnorm(n=n,sd=.1)
y[y< -p] <- -p
y[y> p] <- p

data <- cbind(X,y)
y <- as.matrix(data[,3])
colnames(data) <- c(colnames(X), 'y')
ub <- c(1, 1, 2)
lb <- c(1, -1, -2)

ggplot(data) + geom_point(aes(x=X, y=y), size = 2) +
  ylim(min(y), max(y)) + #coord_fixed(ratio = 1) +
  ggtitle('Data to be regressed') +
  theme_bw(base_size = 12)

regularizer <- function(coeff) coeff%*%coeff/2
regularizer.gr <- function(coeff) coeff
eps <- 1
delt <- 1
gamma <- 1

linrdp <- LinearRegressionDP$new(regularizer, eps, delt, gamma,
                                                 regularizer.gr)
LinearRegressionDP$undebug("fit")
LinearRegressionDP$undebug("postprocess_coeff")
linrdp$fit(X,y,ub,lb)

theta <- linrdp$coeff

grid <- seq(-1,.5,length.out=100)

gridData <- data.frame(X=grid,y=theta[1]+theta[2]*grid)

# Regression line visualization
ggplot() + geom_point(data = data,
                      aes(x=X, y=y),
                      size = 2, show.legend = F) +
  geom_line(data=gridData, aes(x=X,y=y), color='red') + ylim(min(y), max(y)) +
  ggtitle('Regression Line for Linear Regression') +
  theme_bw(base_size = 12)
theta

### END TESTING LINEAR REGRESSION ###
