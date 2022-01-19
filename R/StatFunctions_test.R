cur.dir <- '~/Desktop/DP Research/DP R Package/DPpack/R/';
library(MASS) # For CARS93 dataset

test_meanDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50);
  data2d <- matrix(data=rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(meanDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(meanDP(data1d,1,lower.bounds = c(1,2),upper.bounds=3),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(meanDP(data1d,1,lower.bounds=-3, upper.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds2d:")
  a = tryCatch(meanDP(data2d,1,lower.bounds = -.5,upper.bounds=c(3,3)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds 2d:")
  a = tryCatch(meanDP(data2d,1,lower.bounds=c(-3,-3), upper.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Missing bounds:")
  a = tryCatch(meanDP(data2d,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(meanDP(data2d,1,lower.bounds=c(-3,-3), upper.bounds=c(3,3),
                      mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(meanDP(data2d,1,lower.bounds=c(-3,-3),upper.bounds=c(3,3),
                      mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Simple bounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'unbounded'
  lb = -3
  ub = 3
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Unbounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Simple unbounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -.5
  ub = 3
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Laplace with lower.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram right of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = .5
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Laplace with upper.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Gaussian pDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_varDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50,sd=2);
  data2d <- matrix(data=rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(varDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(varDP(data1d,1,lower.bounds = c(1,2), upper.bounds=6),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(varDP(data1d,1,lower.bounds=-6,upper.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds2d:")
  a = tryCatch(varDP(data2d,1,lower.bounds = -.5,upper.bounds=c(6,6)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds 2d:")
  a = tryCatch(varDP(data2d,1,lower.bounds=c(-6,-6),upper.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Missing bounds:")
  a = tryCatch(varDP(data2d,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(varDP(data2d,1,lower.bounds=c(-6,-6),upper.bounds=c(6,6),
                     mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(varDP(data2d,1,lower.bounds=c(-6,-6),upper.bounds=c(6,6),
                     mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = var(x)
  hist(data,freq=FALSE,main="Simple bounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'unbounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Unbounded
  }

  tv = var(x)
  hist(data,freq=FALSE,main="Simple unbounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -.5
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = var(x)
  hist(data,freq=FALSE,main="Laplace with lower.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = .5
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = var(x)
  hist(data,freq=FALSE,main="Laplace with upper.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = var(x)
  hist(data,freq=FALSE,main="Gaussian pDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = var(x)
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_sdDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50,sd=2);
  data2d <- matrix(data=rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(sdDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(sdDP(data1d,1,lower.bounds = c(1,2), upper.bounds = 6),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(sdDP(data1d,1,lower.bounds=-6,upper.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds2d:")
  a = tryCatch(sdDP(data2d,1,lower.bounds = -.5,upper.bounds=c(6,6)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds 2d:")
  a = tryCatch(sdDP(data2d,1,lower.bounds=c(-6,-6),upper.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Missing bounds:")
  a = tryCatch(varDP(data2d,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(sdDP(data2d,1,lower.bounds=c(-6,-6),upper.bounds=c(6,6),
                    mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(sdDP(data2d,1,lower.bounds=c(-6,-6),upper.bounds=c(6,6),
                    mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = sd(x)
  hist(data,freq=FALSE,main="Simple bounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'unbounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Unbounded
  }

  tv = sd(x)
  hist(data,freq=FALSE,main="Simple unbounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -.5
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = sd(x)
  hist(data,freq=FALSE,main="Laplace with lower.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = .5
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = sd(x)
  hist(data,freq=FALSE,main="Laplace with upper.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = sd(x)
  hist(data,freq=FALSE,main="Gaussian pDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp, ap)$Bounded
  }

  tv = sd(x)
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_covDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50,sd=2);
  data2d <- matrix(data=rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(covDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1 = c(1,2), upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1 = c(1,2),
                     lower.bound2=-5,upper.bound2=7),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         No bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -5
  ub2 = 7
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)$Bounded
  }

  tv = cov(x1,x2)
  hist(data,freq=FALSE,main="Simple bounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- 1;
  ws = 'unbounded'
  lb1 = -6
  ub1 = 6
  lb2 = -5
  ub2 = 7
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)$Unbounded
  }

  tv = cov(x1,x2)
  hist(data,freq=FALSE,main="Simple unbounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- 1;
  ws = 'bounded'
  lb1 = -.5
  ub1 = 6
  lb2 = -5
  ub2 = 7
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)$Bounded
  }

  tv = cov(x1,x2)
  hist(data,freq=FALSE,main="Laplace with lower.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = -.5
  lb2 = -5
  ub2 = 7
  mech = 'laplace'
  delta = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)$Bounded
  }

  tv = cov(x1,x2)
  hist(data,freq=FALSE,main="Laplace with upper.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -5
  ub2 = 7
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)$Bounded
  }

  tv = cov(x1,x2)
  hist(data,freq=FALSE,main="Gaussian pDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -5
  ub2 = 7
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)$Bounded
  }

  tv = cov(x1,x2)
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_histogramDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50)

  print("         No input:")
  a = tryCatch(histogramDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(histogramDP(data1d,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(histogramDP(data1d,1,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, mech, delta, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Simple bounded Laplace");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'unbounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, mech, delta, tdp, an)$Unbounded

  par(mfrow=c(1,2))
  plot(result,main="Simple unbounded Laplace");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, mech, delta, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Gaussian pDP");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, mech, delta, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Gaussian aDP");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Normalized:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- TRUE;
  ws = 'bounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, mech, delta, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Normalized");
  hist(x,freq=FALSE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Allow Negative:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- TRUE;
  ws = 'bounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = TRUE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, mech, delta, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Allow Negative");
  hist(x,freq=FALSE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_tableDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- Cars93$Type;
  y <- Cars93$Origin;

  print("         No input:")
  a = tryCatch(tableDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(tableDP(x,y,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(tableDP(x,y,1,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delta,tdp,an)$Bounded);
  Sys.sleep(5)

  print("         Simple unbounded Laplace:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'unbounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delta,tdp,an)$Unbounded);
  Sys.sleep(5)

  print("         Gaussian pDP:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delta,tdp,an)$Bounded);
  Sys.sleep(5)

  print("         Gaussian aDP:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delta,tdp,an)$Bounded);
  Sys.sleep(5)

  print("         Allow Negative:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  an = TRUE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delta,tdp,an)$Bounded);
  Sys.sleep(5)

  ### END TEST FUNCTIONALITY ###
}

test_pooledVarDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- rnorm(10,sd=2);
  y <- rnorm(15,sd=2);
  z <- rnorm(20,sd=2);

  print("         No input:")
  a = tryCatch(pooledVarDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound = c(1,2),upper.bound=6),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         No bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound=6,
                           mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound=6,
                           mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Simple bounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'unbounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Unbounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Simple unbounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -2
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Laplace with lower.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 2
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Laplace with upper.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Gaussian pDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Approx n.max:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- rnorm(n1,sd=2);
  y <- rnorm(n2,sd=2);
  z <- rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = TRUE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*var(x) + (n2-1)*var(y) + (n3-1)*var(z))/(n1+n2+n3-3);
  hist(data,freq=FALSE,main="Approx n.max");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)
  ### END TEST FUNCTIONALITY ###
}

test_pooledCovDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- matrix(rnorm(20,sd=2),ncol=2);
  y <- matrix(rnorm(30,sd=2),ncol=2);
  z <- matrix(rnorm(40,sd=2),ncol=2);

  print("         No input:")
  a = tryCatch(pooledCovDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound1:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = c(1,2),upper.bound1=6,
                           lower.bound2=-6,upper.bound2=6),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound1:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1=-6,upper.bound1=c(1,2),
                           lower.bound2=-6,upper.bound2=6),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound2:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 =-6,upper.bound1=6,
                           lower.bound2=c(1,2),upper.bound2=6),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound2:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = -6,upper.bound1=6,
                           lower.bound2=-6,upper.bound2=c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = -6,upper.bound1=6,
                           lower.bound2=-6,upper.bound2=6, mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = -6,upper.bound1=6,
                           lower.bound2=-6,upper.bound2=6,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta,
                           type.DP=tdp, approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Simple bounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'unbounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Unbounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Simple unbounded Laplace");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -2
  ub1 = 6
  lb2 = -2
  ub2 = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Laplace with lower.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 2
  lb2 = -6
  ub2 = 2
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Laplace with upper.bounds");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Gaussian pDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'gaussian'
  delta = .5
  tdp = 'aDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Approximate n.max:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'laplace'
  delta = NULL
  tdp = 'pDP'
  anm = TRUE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)$Bounded
  }

  tv = ((n1-1)*cov(x[,1],x[,2])+(n2-1)*cov(y[,1],y[,2])+(n3-1)*cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  hist(data,freq=FALSE,main="Approximate n.max");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_quantileDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- rnorm(50);

  print("         No input:")
  a = tryCatch(quantileDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad quant:")
  a = tryCatch(quantileDP(x,-1,1,-3,3),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(quantileDP(x,.25,1,lower.bound = c(1,2),upper.bound = 3),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(quantileDP(x,.25,1,lower.bound=-3, upper.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         No bounds:")
  a = tryCatch(quantileDP(x,.25,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(quantileDP(x,.25,1,-3,3,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Exponential:")
  x <- rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Simple bounded Exponential");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Exponential:")
  x <- rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'unbounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Unbounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Simple unbounded Exponential");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Exponential with lower.bound:")
  x <- rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'bounded'
  lb = -2
  ub = 3
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Exponential with lower.bound");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Exponential with upper.bounds:")
  x <- rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 2
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Exponential with upper.bound");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Different quants:")
  x <- rnorm(50,sd=2);
  quant=.1
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Different quants");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Different quants:")
  x <- rnorm(50,sd=2);
  quant=.75
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Different quants");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Different quants:")
  x <- rnorm(50,sd=2);
  quant=.9
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = quantile(x,quant)
  hist(data,freq=FALSE,main="Different quants");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_medianDP <- function(){
  #
  #
  #
  # source(paste(cur.dir,'Mechanisms.R',sep=""))
  # source(paste(cur.dir,'StatFunctions.R',sep=""))
  # source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  # Input tested in quantileDP
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  # Most tests done in quantileDP
  print("         Simple Bounded Exponential:")
  x <- rnorm(50,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'exponential'
  delta = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- medianDP(x, eps, lb, ub, ws, mech, delta)$Bounded
  }

  tv = median(x)
  hist(data,freq=FALSE,main="Simple bounded Exponential");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

