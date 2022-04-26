test_meanDP <- function(){
  ### Test input ###
  data1d <- stats::rnorm(50);
  data2d <- matrix(data=stats::rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(meanDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound:")
  a = tryCatch(meanDP(data1d,1,lower.bound = c(1,2),upper.bound=3,
                      type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound:")
  a = tryCatch(meanDP(data1d,1,lower.bound=-3, upper.bound = c(1,2),
                      type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound 2d:")
  a = tryCatch(meanDP(data2d,1,lower.bound=c(-3,-3), upper.bound = -.5,
                      type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound 2d:")
  a = tryCatch(meanDP(data2d,1,lower.bound=-.5,upper.bound=c(3,3),
                      type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Missing bounds:")
  a = tryCatch(meanDP(data2d,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(meanDP(data2d,1,lower.bound=-3, upper.bound=3,
                      mechanism='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(meanDP(data2d,1,lower.bound=-3, upper.bound=3,
                      which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Warning for 'both':")
  a = tryCatch(meanDP(data2d,1,lower.bound=-3, upper.bound=3,
                      which.sensitivity='both',type.DP='pDP'),
               warning=function(w) print(paste("PASS --",w)));
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Simple bounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'unbounded'
  lb = -3
  ub = 3
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Simple unbounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -.5
  ub = 3
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Laplace with lower.bounds");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram right of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bound:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = .5
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Laplace with upper.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Gaussian pDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- data1d;
  eps <- .9;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'Gaussian'
  delta = .01
  tdp = 'aDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Gaussian aDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Multidimensional x:")
  x <- data2d;
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = mean(x)
  graphics::hist(data,freq=FALSE,main="Multidimensional x");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_varDP <- function(){
  ### Test input ###
  data1d <- stats::rnorm(50,sd=2);
  data2d <- matrix(data=stats::rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(varDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad x:")
  a = tryCatch(varDP(data2d,1,lower.bound = -6, upper.bound=6,
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound:")
  a = tryCatch(varDP(data1d,1,lower.bound = c(1,2), upper.bound=6,
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound:")
  a = tryCatch(varDP(data1d,1,lower.bound=-6,upper.bound = c(1,2),
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Missing bounds:")
  a = tryCatch(varDP(data1d,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(varDP(data1d,1,lower.bound=-6,upper.bound=6,
                     mechanism='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(varDP(data1d,1,lower.bound=-6, upper.bound=6,
                      which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Warning for 'both':")
  a = tryCatch(varDP(data1d,1,lower.bound=-6, upper.bound=6,
                      which.sensitivity='both',type.DP='pDP'),
               warning=function(w) print(paste("PASS --",w)));
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::var(x)
  graphics::hist(data,freq=FALSE,main="Simple bounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'unbounded'
  lb = -6
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::var(x)
  graphics::hist(data,freq=FALSE,main="Simple unbounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bound:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -.5
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::var(x)
  graphics::hist(data,freq=FALSE,main="Laplace with lower.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bound:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = .5
  mech = 'Laplace'
  delta = 0
  tdp = 'pDp'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::var(x)
  graphics::hist(data,freq=FALSE,main="Laplace with upper.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::var(x)
  graphics::hist(data,freq=FALSE,main="Gaussian pDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- data1d;
  eps <- .9;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Gaussian'
  delta = .5
  tdp = 'aDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::var(x)
  graphics::hist(data,freq=FALSE,main="Gaussian aDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_sdDP <- function(){
  ### Test input ###
  data1d <- stats::rnorm(50,sd=2);
  data2d <- matrix(data=stats::rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(sdDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad x:")
  a = tryCatch(sdDP(data2d,1,lower.bound = -6, upper.bound=6,
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound:")
  a = tryCatch(sdDP(data1d,1,lower.bound = c(1,2), upper.bound=6,
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound:")
  a = tryCatch(sdDP(data1d,1,lower.bound=-6,upper.bound = c(1,2),
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Missing bounds:")
  a = tryCatch(sdDP(data1d,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(sdDP(data1d,1,lower.bound=-6,upper.bound=6,
                     mechanism='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(sdDP(data1d,1,lower.bound=-6, upper.bound=6,
                     which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Warning for 'both':")
  a = tryCatch(sdDP(data1d,1,lower.bound=-6, upper.bound=6,
                     which.sensitivity='both',type.DP='pDP'),
               warning=function(w) print(paste("PASS --",w)));
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::sd(x)
  graphics::hist(data,freq=FALSE,main="Simple bounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'unbounded'
  lb = -6
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::sd(x)
  graphics::hist(data,freq=FALSE,main="Simple unbounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bound:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -.5
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::sd(x)
  graphics::hist(data,freq=FALSE,main="Laplace with lower.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Laplace with upper.bound:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = .5
  mech = 'Laplace'
  delta = 0
  tdp = 'pDp'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::sd(x)
  graphics::hist(data,freq=FALSE,main="Laplace with upper.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram left of line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::sd(x)
  graphics::hist(data,freq=FALSE,main="Gaussian pDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x <- data1d;
  eps <- .9;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Gaussian'
  delta = .5
  tdp = 'aDP'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, lb, ub, ws, mech, delta, tdp)
  }

  tv = stats::sd(x)
  graphics::hist(data,freq=FALSE,main="Gaussian aDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_covDP <- function(){
  ### Test input ###
  data1d <- stats::rnorm(50,sd=2);
  data2d <- matrix(data=stats::rnorm(1000),ncol=2);

  print("         No input:")
  a = tryCatch(covDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1 = c(1,2), upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1 = c(1,2),
                     lower.bound2=-5,upper.bound2=7,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         No bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7,mechanism='abc',
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7,which.sensitivity='abc',
                     type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Warning for 'both':")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1=-6,upper.bound1=6,
                     lower.bound2=-5,upper.bound2=7,which.sensitivity='both',
                     type.DP='pDP'),
               warning=function(w) print(paste("PASS --",w)));
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDp'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)
  }

  tv = stats::cov(x1,x2)
  graphics::hist(data,freq=FALSE,main="Simple bounded Laplace");
  graphics::abline(v=tv,col='blue')
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)
  }

  tv = stats::cov(x1,x2)
  graphics::hist(data,freq=FALSE,main="Simple unbounded Laplace");
  graphics::abline(v=tv,col='blue')
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)
  }

  tv = stats::cov(x1,x2)
  graphics::hist(data,freq=FALSE,main="Laplace with lower.bounds");
  graphics::abline(v=tv,col='blue')
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
  mech = 'Laplace'
  delta = 0
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)
  }

  tv = stats::cov(x1,x2)
  graphics::hist(data,freq=FALSE,main="Laplace with upper.bounds");
  graphics::abline(v=tv,col='blue')
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
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)
  }

  tv = stats::cov(x1,x2)
  graphics::hist(data,freq=FALSE,main="Gaussian pDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  x1 <- data1d;
  x2 = data1d+1;
  eps <- .9;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -5
  ub2 = 7
  mech = 'Gaussian'
  delta = .01
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, lb1, ub1, lb2, ub2, ws, mech, delta, tdp)
  }

  tv = stats::cov(x1,x2)
  graphics::hist(data,freq=FALSE,main="Gaussian aDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_histogramDP <- function(){
  ### Test input ###
  data1d <- stats::rnorm(50)

  print("         No input:")
  a = tryCatch(histogramDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(histogramDP(data1d,1,mechanism='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(histogramDP(data1d,1,which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- stats::rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges", normal, ws, mech, delta, tdp, an)

  # graphics::par(mfrow=c(1,2))
  plot(result,main="Simple bounded Laplace");
  graphics::hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)
  # graphics::par(mfrow=c(1,1))

  print("         Simple unbounded Laplace:")
  x <- stats::rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'unbounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges", normal, ws, mech, delta, tdp, an)

  # graphics::par(mfrow=c(1,2))
  plot(result,main="Simple unbounded Laplace");
  graphics::hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)
  # graphics::par(mfrow=c(1,1))

  print("         Gaussian pDP:")
  x <- stats::rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges", normal, ws, mech, delta, tdp, an)

  # graphics::par(mfrow=c(1,2))
  plot(result,main="Gaussian pDP");
  graphics::hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)
  # graphics::par(mfrow=c(1,1))

  print("         Gaussian aDP:")
  x <- stats::rnorm(100,mean=3,sd=2);
  eps <- .9;
  normal <- FALSE;
  ws = 'bounded'
  mech = 'Gaussian'
  delta = .01
  tdp = 'aDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges", normal, ws, mech, delta, tdp, an)

  # graphics::par(mfrow=c(1,2))
  plot(result,main="Gaussian aDP");
  graphics::hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)
  # graphics::par(mfrow=c(1,1))

  print("         Normalized:")
  x <- stats::rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- TRUE;
  ws = 'bounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges", normal, ws, mech, delta, tdp, an)

  # graphics::par(mfrow=c(1,2))
  plot(result,main="Normalized");
  graphics::hist(x,freq=FALSE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)
  # graphics::par(mfrow=c(1,1))

  print("         Allow Negative:")
  x <- stats::rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- TRUE;
  ws = 'bounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = TRUE;

  result <- histogramDP(x, eps, "Sturges", normal, ws, mech, delta, tdp, an)

  # graphics::par(mfrow=c(1,2))
  plot(result,main="Allow Negative");
  graphics::hist(x,freq=FALSE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)
  # graphics::par(mfrow=c(1,1))

  ### END TEST FUNCTIONALITY ###
}

test_tableDP <- function(){
  ### Test input ###
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  z <- MASS::Cars93$AirBags;

  print("         No input:")
  a = tryCatch(tableDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(tableDP(x,y,eps=1,mechanism='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(tableDP(x,y,eps=1,which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps=eps,which.sensitivity=ws,mechanism=mech,delta=delta,
                type.DP=tdp,allow.negative=an));
  Sys.sleep(5)

  print("         Simple unbounded Laplace:")
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  eps <- 1;
  ws = 'unbounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps=eps,which.sensitivity=ws,mechanism=mech,delta=delta,
                type.DP=tdp,allow.negative=an));
  Sys.sleep(5)

  print("         Gaussian pDP:")
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps=eps,which.sensitivity=ws,mechanism=mech,delta=delta,
                type.DP=tdp,allow.negative=an));
  Sys.sleep(5)

  print("         Gaussian aDP:")
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  eps <- .9;
  ws = 'bounded'
  mech = 'Gaussian'
  delta = .01
  tdp = 'aDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps=eps,which.sensitivity=ws,mechanism=mech,delta=delta,
                type.DP=tdp,allow.negative=an));
  Sys.sleep(5)

  print("         Allow Negative:")
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = TRUE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps=eps,which.sensitivity=ws,mechanism=mech,delta=delta,
                type.DP=tdp,allow.negative=an));
  Sys.sleep(5)

  print("         More than 2 inputs:")
  x <- MASS::Cars93$Type;
  y <- MASS::Cars93$Origin;
  z <- MASS::Cars93$AirBags;
  eps <- 1;
  ws = 'bounded'
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y,z))
  print("Sanitized:")
  print(tableDP(x,y,z,eps=eps,which.sensitivity=ws,mechanism=mech,delta=delta,
                type.DP=tdp,allow.negative=an));
  Sys.sleep(5)

  ### END TEST FUNCTIONALITY ###
}

test_pooledVarDP <- function(){
  ### Test input ###
  x <- stats::rnorm(10,sd=2);
  y <- stats::rnorm(15,sd=2);
  z <- stats::rnorm(20,sd=2);

  print("         No input:")
  a = tryCatch(pooledVarDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound = c(1,2),upper.bound=6,
                           type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound = c(1,2),
                           type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         No bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound=6,
                           mechanism='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound=6,
                           which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Warning for 'both':")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound=-6,upper.bound=6,
                           which.sensitivity='both',type.DP='pDP'),
               warning=function(w) print(paste("PASS --",w)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Simple bounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- 1;
  ws = 'unbounded'
  lb = -6
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Simple unbounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -2
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Laplace with lower.bounds");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 2
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Laplace with upper.bounds");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Gaussian pDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- .9;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Gaussian'
  delta = .01
  tdp = 'aDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Gaussian aDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Approx n.max:")
  n1 = 10
  n2 = 15
  n3 = 20
  x <- stats::rnorm(n1,sd=2);
  y <- stats::rnorm(n2,sd=2);
  z <- stats::rnorm(n3,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = TRUE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps,
                           lower.bound=lb, upper.bound=ub, which.sensitivity=ws,
                           mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::var(x) + (n2-1)*stats::var(y) + (n3-1)*stats::var(z))/(n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Approx n.max");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)
  ### END TEST FUNCTIONALITY ###
}

test_pooledCovDP <- function(){
  ### Test input ###
  x <- matrix(stats::rnorm(20,sd=2),ncol=2);
  y <- matrix(stats::rnorm(30,sd=2),ncol=2);
  z <- matrix(stats::rnorm(40,sd=2),ncol=2);

  print("         No input:")
  a = tryCatch(pooledCovDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound1:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = c(1,2),upper.bound1=6,
                           lower.bound2=-6,upper.bound2=6,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound1:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1=-6,upper.bound1=c(1,2),
                           lower.bound2=-6,upper.bound2=6,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound2:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 =-6,upper.bound1=6,
                           lower.bound2=c(1,2),upper.bound2=6,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound2:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = -6,upper.bound1=6,
                           lower.bound2=-6,upper.bound2=c(1,2),type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = -6,upper.bound1=6,
                           lower.bound2=-6,upper.bound2=6, mechanism='abc',
                           type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1=-6,upper.bound1=6,
                           lower.bound2=-6,upper.bound2=6,
                           which.sensitivity='abc',type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta,
                           type.DP=tdp, approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Simple bounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Laplace:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'unbounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Simple unbounded Laplace");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -2
  ub1 = 6
  lb2 = -2
  ub2 = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Laplace with lower.bounds");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 2
  lb2 = -6
  ub2 = 2
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Laplace with upper.bounds");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'Gaussian'
  delta = .01
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Gaussian pDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Gaussian aDP:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- .9;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'Gaussian'
  delta = .01
  tdp = 'aDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Gaussian aDP");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Approximate n.max:")
  n1 = 10;
  n2 = 15;
  n3 = 20;
  x <- matrix(stats::rnorm(2*n1,sd=2),ncol=2);
  y <- matrix(stats::rnorm(2*n2,sd=2),ncol=2);
  z <- matrix(stats::rnorm(2*n3,sd=2),ncol=2);
  eps <- 1;
  ws = 'bounded'
  lb1 = -6
  ub1 = 6
  lb2 = -6
  ub2 = 6
  mech = 'Laplace'
  delta = 0
  tdp = 'pDP'
  anm = TRUE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           which.sensitivity=ws, mechanism=mech, delta=delta, type.DP=tdp,
                           approx.n.max=anm)
  }

  tv = ((n1-1)*stats::cov(x[,1],x[,2])+(n2-1)*stats::cov(y[,1],y[,2])+(n3-1)*stats::cov(z[,1],z[,2]))/
    (n1+n2+n3-3);
  graphics::hist(data,freq=FALSE,main="Approximate n.max");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_quantileDP <- function(){
  ### Test input ###
  x <- stats::rnorm(50);

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

  print("         Bad which.sensitivity:")
  a = tryCatch(quantileDP(x,.25,1,-3,3,which.sensitivity='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Warning for 'both':")
  a = tryCatch(quantileDP(x,.25,1,-3,3,which.sensitivity='both'),
               warning=function(w) print(paste("PASS --",w)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Exponential:")
  x <- stats::rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Simple bounded Exponential");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Simple unbounded Exponential:")
  x <- stats::rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'unbounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Simple unbounded Exponential");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Exponential with lower.bound:")
  x <- stats::rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'bounded'
  lb = -2
  ub = 3
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Exponential with lower.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Exponential with upper.bounds:")
  x <- stats::rnorm(50,sd=2);
  quant=.25
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 2
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Exponential with upper.bound");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Different quants:")
  x <- stats::rnorm(50,sd=2);
  quant=.1
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Different quants");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Different quants:")
  x <- stats::rnorm(50,sd=2);
  quant=.75
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Different quants");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  print("         Different quants:")
  x <- stats::rnorm(50,sd=2);
  quant=.9
  eps <- 1;
  ws = 'bounded'
  lb = -3
  ub = 3
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, lb, ub, ws, mech)
  }

  tv = stats::quantile(x,quant)
  graphics::hist(data,freq=FALSE,main="Different quants");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

test_medianDP <- function(){
  ### Test input ###
  # Input tested in quantileDP
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  # Most tests done in quantileDP
  print("         Simple Bounded Exponential:")
  x <- stats::rnorm(50,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = -6
  ub = 6
  mech = 'exponential'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- medianDP(x, eps, lb, ub, ws, mech)
  }

  tv = stats::median(x)
  graphics::hist(data,freq=FALSE,main="Simple bounded Exponential");
  graphics::abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}
