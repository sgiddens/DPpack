cur.dir <- '~/Desktop/DP Research/DP R Package/DPpack/R/';
library(MASS) # For CARS93 dataset

test_meanDP <- function(){
  #
  #
  #
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50);
  data2d <- matrix(data=rnorm(1000),ncol=100);

  print("         No input:")
  a = tryCatch(meanDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(meanDP(data1d,1,lower.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(meanDP(data1d,1,upper.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds2d:")
  a = tryCatch(meanDP(data2d,1,lower.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds 2d:")
  a = tryCatch(meanDP(data2d,1,upper.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(meanDP(data2d,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
  a = tryCatch(meanDP(data2d,1,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Unbounded
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
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = .5
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- meanDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
  }

  tv = mean(x)
  hist(data,freq=FALSE,main="Gaussian aDP");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###

  ### TEST DEFINITION ###
  # vals = numeric(100)
  # for (j in 1:100){
  #   x = c(-4,-4);
  #   y = c(-4,4);
  #   eps = .5
  #   ws = 'bounded';
  #   lb = -4
  #   ub = 4
  #   mech = 'laplace'
  #   delt = NULL
  #   tdp = 'pDP'
  #   ap = NULL
  #   n = 10000
  #   data.x = numeric(n)
  #   data.y = numeric(n)
  #   for (i in 1:n){
  #     data.x[i] <- meanDP(x,eps,ws,lb,ub,mech,delt,tdp,ap)$Bounded;
  #     data.y[i] <- meanDP(y,eps,ws,lb,ub,mech,delt,tdp,ap)$Bounded;
  #   }
  #   S = c(-.5,0)
  #   p.M.x = sum(data.x>S[1] & data.x<S[2])/n
  #   p.M.y = sum(data.y>S[1] & data.y<S[2])/n
  #   vals[j] <- ((log(p.M.y) - log(p.M.x)) <= eps)
  # }
  # print(sum(vals)/100)
  ### END TEST DEFINITION ###
}

test_varDP <- function(){
  #
  #
  #
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50,sd=2);
  data2d <- matrix(data=rnorm(1000),ncol=100);

  print("         No input:")
  a = tryCatch(varDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(varDP(data1d,1,lower.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(varDP(data1d,1,upper.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds2d:")
  a = tryCatch(varDP(data2d,1,lower.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds 2d:")
  a = tryCatch(varDP(data2d,1,upper.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(varDP(data2d,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
  a = tryCatch(varDP(data2d,1,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Unbounded
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
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = .5
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- varDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50,sd=2);
  data2d <- matrix(data=rnorm(1000),ncol=100);

  print("         No input:")
  a = tryCatch(sdDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(sdDP(data1d,1,lower.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(sdDP(data1d,1,upper.bounds = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds2d:")
  a = tryCatch(sdDP(data2d,1,lower.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds 2d:")
  a = tryCatch(sdDP(data2d,1,upper.bounds = -.5),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(sdDP(data2d,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
  a = tryCatch(sdDP(data2d,1,mechanism='gaussian'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded Laplace:")
  x <- data1d;
  eps <- 1;
  ws = 'bounded'
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Unbounded
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
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = .5
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- sdDP(x, eps, ws, lb, ub, mech, delt, tdp, ap)$Bounded
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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50,sd=2);
  data2d <- matrix(data=rnorm(1000),ncol=100);

  print("         No input:")
  a = tryCatch(covDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,lower.bound1 = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(covDP(data1d,data1d+1,1,upper.bound1 = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(covDP(data1d,data1d+1,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
  a = tryCatch(covDP(data1d,data1d+1,1,mechanism='gaussian'),
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, ws, lb1, ub1, lb2, ub2, mech, delt, tdp)$Bounded
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, ws, lb1, ub1, lb2, ub2, mech, delt, tdp)$Unbounded
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
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, ws, lb1, ub1, lb2, ub2, mech, delt, tdp)$Bounded
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
  lb1 = NULL
  ub1 = -.5
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDp'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, ws, lb1, ub1, lb2, ub2, mech, delt, tdp)$Bounded
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, ws, lb1, ub1, lb2, ub2, mech, delt, tdp)$Bounded
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  ap = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- covDP(x1, x2, eps, ws, lb1, ub1, lb2, ub2, mech, delt, tdp)$Bounded
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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  data1d <- rnorm(50)

  print("         No input:")
  a = tryCatch(histogramDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(histogramDP(data1d,1,lower.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(histogramDP(data1d,1,upper.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(histogramDP(data1d,1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
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
  lb = -10
  ub = 15
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

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
  lb = -10
  ub = 15
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Unbounded

  par(mfrow=c(1,2))
  plot(result,main="Simple unbounded Laplace");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Laplace with lower.bounds:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  lb = 1
  ub = 15
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Laplace with lower.bounds");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Laplace with upper.bounds:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  lb = -9
  ub = 5
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

  par(mfrow=c(1,2))
  plot(result,main="Laplace with upper.bounds");
  hist(x,freq=TRUE,main="Original");
  print("Verify similarity...")
  Sys.sleep(2)

  print("         Gaussian pDP:")
  x <- rnorm(100,mean=3,sd=2);
  eps <- 1;
  normal <- FALSE;
  ws = 'bounded'
  lb = -9
  ub = 15
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

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
  lb = -9
  ub = 15
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

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
  lb = -9
  ub = 15
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

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
  lb = -9
  ub = 15
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = TRUE;

  result <- histogramDP(x, eps, "Sturges",normal, ws, lb, ub, mech, delt, tdp, an)$Bounded

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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

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

  print("         Bad delt:")
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
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delt,tdp,an)$Bounded);
  Sys.sleep(5)

  print("         Simple unbounded Laplace:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'unbounded'
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delt,tdp,an)$Unbounded);
  Sys.sleep(5)

  print("         Gaussian pDP:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delt,tdp,an)$Bounded);
  Sys.sleep(5)

  print("         Gaussian aDP:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  an = FALSE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delt,tdp,an)$Bounded);
  Sys.sleep(5)

  print("         Allow Negative:")
  x <- Cars93$Type;
  y <- Cars93$Origin;
  eps <- 1;
  ws = 'bounded'
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  an = TRUE;

  print("Verify similarity...")
  print("Original:")
  print(table(x,y))
  print("Sanitized:")
  print(tableDP(x,y,eps,ws,mech,delt,tdp,an)$Bounded);
  Sys.sleep(5)

  ### END TEST FUNCTIONALITY ###
}

test_pooledVarDP <- function(){
  #
  #
  #
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- rnorm(10,sd=2);
  y <- rnorm(15,sd=2);
  z <- rnorm(20,sd=2);

  print("         No input:")
  a = tryCatch(pooledVarDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,lower.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,upper.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
  a = tryCatch(pooledVarDP(x,y,z,eps=1,mechanism='gaussian'),
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
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb = NULL
  ub = 2
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb = NULL
  ub = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  anm = FALSE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb = NULL
  ub = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = TRUE
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledVarDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound=lb, upper.bound=ub,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- matrix(rnorm(20,sd=2),ncol=2);
  y <- matrix(rnorm(30,sd=2),ncol=2);
  z <- matrix(rnorm(40,sd=2),ncol=2);

  print("         No input:")
  a = tryCatch(pooledCovDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound1:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound1 = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound1:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,upper.bound1 = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bound2:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,lower.bound2 = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bound2:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,upper.bound2 = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,mechanism='abc'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delt:")
  a = tryCatch(pooledCovDP(x,y,z,eps=1,mechanism='gaussian'),
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
                           approx.n.max=anm)$Bounded
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  ub1 = NULL
  lb2 = -2
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb1 = NULL
  ub1 = 2
  lb2 = NULL
  ub2 = 2
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'pDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'gaussian'
  delt = .5
  tdp = 'aDP'
  anm = FALSE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  lb1 = NULL
  ub1 = NULL
  lb2 = NULL
  ub2 = NULL
  mech = 'laplace'
  delt = NULL
  tdp = 'pDP'
  anm = TRUE;
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- pooledCovDP(x,y,z, eps=eps, which.sensitivity=ws,
                           lower.bound1=lb1, upper.bound1=ub1,
                           lower.bound2=lb2, upper.bound2=ub2,
                           mechanism=mech, delt=delt, type.DP=tdp,
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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  x <- rnorm(50);

  print("         No input:")
  a = tryCatch(quantileDP(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad quant:")
  a = tryCatch(quantileDP(x,-1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad lower.bounds:")
  a = tryCatch(quantileDP(x,.25,1,lower.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad upper.bounds:")
  a = tryCatch(quantileDP(x,.25,1,upper.bound = c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad mechanism:")
  a = tryCatch(quantileDP(x,1,mechanism='abc'),
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
  lb = NULL
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Unbounded
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
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Bounded
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
  lb = NULL
  ub = 2
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Bounded
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
  lb = NULL
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- quantileDP(x, quant, eps, ws, lb, ub, mech, delt)$Bounded
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
  source(paste(cur.dir,'Mechanisms.R',sep=""))
  source(paste(cur.dir,'StatFunctions.R',sep=""))
  source(paste(cur.dir,'DataAccess.R',sep=""))

  ### Test input ###
  # Input tested in quantileDP
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  # Most tests done in quantileDP
  print("         Simple Bounded Exponential:")
  x <- rnorm(50,sd=2);
  eps <- 1;
  ws = 'bounded'
  lb = NULL
  ub = NULL
  mech = 'exponential'
  delt = NULL
  n = 10000;
  data <- numeric(n);
  for (i in 1:n){
    data[i] <- medianDP(x, eps, ws, lb, ub, mech, delt)$Bounded
  }

  tv = median(x)
  hist(data,freq=FALSE,main="Simple bounded Exponential");
  abline(v=tv,col='blue')
  print("Verify histogram centered at line...")
  Sys.sleep(2)

  ### END TEST FUNCTIONALITY ###
}

