test_Laplace <- function(){
  ### Test input ###
  print("         No input:")
  a = tryCatch(LaplaceMechanism(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad true.values:");
  a = tryCatch(LaplaceMechanism('a',1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad eps:")
  a = tryCatch(LaplaceMechanism(0,'a',1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,-1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,0,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad sensitivities:")
  a = tryCatch(LaplaceMechanism(0,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,-1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(c(0,1,2),1,c(0,1)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad alloc.proportions:")
  a = tryCatch(LaplaceMechanism(0,1,1,alloc.proportions='a'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,1,alloc.proportions=-1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(c(0,1),1,1,alloc.proportions=c(1,2,3)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple example:");
  tv = 0;
  eps = 1;
  sens = 1;
  ap = NULL;
  th.s = sens/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Simple example");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,rmutil::dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values:");
  tv = 5;
  eps = 1;
  sens = 1;
  ap = NULL;
  th.s = sens/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing true.values");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,rmutil::dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values/epsilon:");
  tv = -5;
  eps = .1;
  sens = 1;
  ap = NULL;
  th.s = sens/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing true.values/epsilon");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,rmutil::dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing epsilon/sensitivity:");
  tv = -5;
  eps = .1;
  sens = .5;
  ap = NULL;
  th.s = sens/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing epsilon/sensitivity");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,rmutil::dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using multiple:");
  tv = c(-5,0,10);
  eps = 1;
  sens = c(.5,1,2);
  ap = NULL;
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = LaplaceMechanism(tv,eps,sens,ap);
  }
  th.s = sum(sens)/eps;
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Using multiple:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,rmutil::dlaplace(x,m=tv[j],s=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Single sensitivity for multiple values:")
  tv = c(-5,0)
  eps = 1
  sens = 1
  ap = NULL
  n = 10000
  data = matrix(NaN,nrow=n,ncol=length(tv))
  for (i in 1:n){
    data[i,] = LaplaceMechanism(tv,eps,sens,ap)
  }
  th.s = sens/eps
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Single sensitivity:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,rmutil::dlaplace(x,m=tv[j],s=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Using allocation:");
  tv = c(-5,0);
  eps = 1;
  sens = c(.5,1);
  ap = c(1,9);
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = LaplaceMechanism(tv,eps,sens,ap);
  }
  ap = ap/sum(ap);
  th.s = sens/(ap*eps);
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Using allocation:",j));
    x = seq(tv[j]-5*th.s[j], tv[j]+5*th.s[j],.1);
    graphics::lines(x,rmutil::dlaplace(x,m=tv[j],s=th.s[j]));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }
  ### END TEST FUNCTIONALITY ###
}

test_Gaussian <- function(){
  ### Test input ###
  print("         No input:")
  a = tryCatch(GaussianMechanism(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad true.values:");
  a = tryCatch(GaussianMechanism('a',1,1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad eps:")
  a = tryCatch(GaussianMechanism(0,'a',1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,-1,1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,0,1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,3,1,1,type.DP='aDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(GaussianMechanism(0,1,'a',1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,-1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,0,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad sensitivities:")
  a = tryCatch(GaussianMechanism(0,1,1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,c(1,2),type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,-1,type.DP='pDP'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad type.DP:")
  a = tryCatch(GaussianMechanism(0,1,1,1,type.DP='pdp'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad alloc.proportions:")
  a = tryCatch(GaussianMechanism(0,1,1,1,type.DP='pDP',alloc.proportions='a'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,1,type.DP='pDP',alloc.proportions=-1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(c(0,1),1,1,1,type.DP='pDP',alloc.proportions=c(1,2,3)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple example (pDP):");
  tv = 0;
  eps = .5;
  delta = .01;
  sens = 1;
  dp = 'pDP'
  ap = NULL;
  th.s = sens*(sqrt(stats::qnorm(delta/2)^2+2*eps)-stats::qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  graphics::hist(data,freq=FALSE,main="Simple example (pDP)");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple example (aDP):");
  tv = 0;
  eps = .5;
  delta = .01
  sens = 1;
  dp = 'aDP'
  ap = NULL;
  th.s = sens*sqrt(2*log(1.25/delta))/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  graphics::hist(data,freq=FALSE,main="Simple bounded (aDP)");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values:");
  tv = 5;
  eps = .5;
  delta = .01
  sens = 1;
  dp = 'pDP';
  ap = NULL;
  th.s = sens*(sqrt(stats::qnorm(delta/2)^2+2*eps)-stats::qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing true.values");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values/epsilon:");
  tv = -5;
  eps = .1;
  delta = .01
  sens = 1;
  dp = 'aDP';
  ap = NULL;
  th.s = sens*sqrt(2*log(1.25/delta))/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing true.values/epsilon");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing delta/sensitivities:");
  tv = -5;
  eps = .5;
  delta = .1
  sens = 3;
  dp = 'pDP';
  ap = NULL;
  th.s = sens*(sqrt(stats::qnorm(delta/2)^2+2*eps)-stats::qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing delta/sensitivities");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using multiple:");
  tv = c(-5,0,10);
  eps = .5
  delta = .01
  sens = c(.5,1,2);
  dp = 'aDP';
  ap = NULL;
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  th.s = sqrt(sum(sens^2))*sqrt(2*log(1.25/delta))/eps;
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Using multiple:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,stats::dnorm(x,m=tv[j],sd=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Single sensitivity for multiple values:")
  tv = c(-5,0)
  eps = .5;
  delta = .01
  sens = 1
  dp = 'pDP';
  ap = NULL
  n = 10000
  data = matrix(NaN,nrow=n,ncol=length(tv))
  for (i in 1:n){
    data[i,] = GaussianMechanism(tv,eps,delta,sens,dp,ap)
  }
  th.s = sens*(sqrt(stats::qnorm(delta/2)^2+2*eps)-stats::qnorm((delta)/2))/(2*eps)
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Single sensitivity:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,stats::dnorm(x,m=tv[j],sd=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Using allocation:");
  tv = c(-5,0);
  eps = .5;
  delta = .01
  sens = c(.5,1);
  dp = 'pDP';
  ap = c(1,9);
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = GaussianMechanism(tv,eps,delta,sens,dp,ap);
  }
  ap = ap/sum(ap);
  th.s = sens*(sqrt(stats::qnorm((ap*delta)/2)^2+2*(ap*eps))-stats::qnorm((ap*delta)/2))/(2*(ap*eps));
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Using allocation:",j));
    x = seq(tv[j]-5*th.s[j], tv[j]+5*th.s[j],.1);
    graphics::lines(x,stats::dnorm(x,m=tv[j],sd=th.s[j]));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }
  ### END TEST FUNCTIONALITY ###
}

test_analytic_Gaussian <- function(){
  ### Test input ###
  print("         No input:")
  a = tryCatch(AnalyticGaussianMechanism(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad true.values:");
  a = tryCatch(AnalyticGaussianMechanism('a',1,1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad eps:")
  a = tryCatch(AnalyticGaussianMechanism(0,'a',1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,-1,1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,0,1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(AnalyticGaussianMechanism(0,1,'a',1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,1,-1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,1,0,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad sensitivities:")
  a = tryCatch(AnalyticGaussianMechanism(0,1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,1,1,c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,1,1,-1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad alloc.proportions:")
  a = tryCatch(AnalyticGaussianMechanism(0,1,1,1,alloc.proportions='a'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(0,1,1,1,alloc.proportions=-1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(AnalyticGaussianMechanism(c(0,1),1,1,1,alloc.proportions=c(1,2,3)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple example:");
  tv = 0;
  eps = .5;
  delta = .01;
  sens = 1;
  ap = NULL;
  th.s = calibrateAnalyticGaussianMechanism(eps, delta, sens)
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Simple example");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values:");
  tv = 5;
  eps = .5;
  delta = .01
  sens = 1;
  ap = NULL;
  th.s = calibrateAnalyticGaussianMechanism(eps, delta, sens)
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing true.values");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values/epsilon:");
  tv = -5;
  eps = 2;
  delta = .01
  sens = 1;
  ap = NULL;
  th.s = calibrateAnalyticGaussianMechanism(eps, delta, sens)
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing true.values/epsilon");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing delta/sensitivities:");
  tv = -5;
  eps = .5;
  delta = .1
  sens = 3;
  ap = NULL;
  th.s = calibrateAnalyticGaussianMechanism(eps, delta, sens)
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap);
  }
  graphics::hist(data,freq=FALSE,main="Changing delta/sensitivities");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  graphics::lines(x,stats::dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using multiple:");
  tv = c(-5,0,10);
  eps = .5
  delta = .01
  sens = c(.5,1,2);
  ap = NULL;
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap);
  }
  th.s = calibrateAnalyticGaussianMechanism(eps, delta, sqrt(sum(sens^2)))
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Using multiple:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,stats::dnorm(x,m=tv[j],sd=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Single sensitivity for multiple values:")
  tv = c(-5,0)
  eps = .5;
  delta = .01
  sens = 1
  ap = NULL
  n = 10000
  data = matrix(NaN,nrow=n,ncol=length(tv))
  for (i in 1:n){
    data[i,] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap)
  }
  th.s = calibrateAnalyticGaussianMechanism(eps, delta, sens)
  for (j in 1:length(tv)){
    graphics::hist(data[,j],freq=FALSE,main=paste("Single sensitivity:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,stats::dnorm(x,m=tv[j],sd=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Using allocation:");
  tv = c(-5,0);
  eps = .5;
  delta = .01
  sens = c(.5,1);
  ap = c(1,9);
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = AnalyticGaussianMechanism(tv,eps,delta,sens,ap);
  }
  ap = ap/sum(ap);
  for (j in 1:length(tv)){
    th.s = calibrateAnalyticGaussianMechanism(ap[j]*eps, ap[j]*delta, sens[j])
    graphics::hist(data[,j],freq=FALSE,main=paste("Using allocation:",j));
    x = seq(tv[j]-5*th.s, tv[j]+5*th.s,.1);
    graphics::lines(x,stats::dnorm(x,m=tv[j],sd=th.s));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }
  ### END TEST FUNCTIONALITY ###
}

test_Exponential <- function(){
  ### Test input ###
  print("         No input:")
  a = tryCatch(ExponentialMechanism(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad utility:");
  a = tryCatch(ExponentialMechanism('a',1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad eps:")
  a = tryCatch(ExponentialMechanism(0,'a',1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,-1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,0,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad sensitivity:")
  a = tryCatch(ExponentialMechanism(0,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,-1),
               error=function(e) print(paste("PASS --",e)));
  print("")

  print("         Bad measure:")
  a = tryCatch(ExponentialMechanism(c(0,1,2),1,1,measure = c(1,1)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,1,measure = -1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad candidates:")
  a = tryCatch(ExponentialMechanism(c(0,1,2),1,1,candidates = c(-1,1)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple example:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  sens = 1;
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,sens,m,c);
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*sens));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Simple example");
  graphics::lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Changing utility:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = 1;
  sens = 1;
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,sens,m,c);
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*sens));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Different Utility");
  graphics::lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Changing utility/epsilon:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = .1;
  sens = 1;
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,sens,m,c);
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*sens));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Different utility/epsilon");
  graphics::lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Changing epsilon/sensitivity:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = .1;
  sens = 5;
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,sens,m,c);
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*sens));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Different epsilon/sensitivity");
  graphics::lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Using measure:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  sens = 1;
  m = c(.5,1,2,.5,5,1);
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,sens,m,c);
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*sens));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Using Measure");
  graphics::lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Returning Index:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  sens = 1;
  m = NULL;
  c = NULL;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,sens,m,c);
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*sens));
  th.probs = th.probs/sum(th.probs);
  plot(1:length(Z),table(data)/sum(table(data)),main="Returning Index");
  graphics::lines(1:length(Z),th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)
  ### END TEST FUNCTIONALITY ###
}
