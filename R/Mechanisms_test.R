library(ggplot2)
library(rmutil)
cur.dir <- '~/Desktop/DP Research/DP R Package/DPpack/R/';

test_Laplace <- function(){
  #
  #
  #
  source(paste(cur.dir,'Mechanisms.R',sep=""))

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

  print("         Bad bounded.sensitivities:")
  a = tryCatch(LaplaceMechanism(0,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(c(0,1),1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,bounded.sensitivities=-1),
               error=function(e) print(paste("PASS --",e)));
  print("")

  print("         Bad unbounded.sensitivities:")
  a = tryCatch(LaplaceMechanism(0,1,which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,unbounded.sensitivities=c(1,2),
                                which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(c(0,1),1,unbounded.sensitivities=1,
                                which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,unbounded.sensitivities=-1,
                                which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(LaplaceMechanism(0,1,which.sensitivity='unbnd'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(0,1,unbounded.sensitivities=c(1,2),
                                which.sensitivity='bounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(c(0,1),1,unbounded.sensitivities=1,
                                which.sensitivity='both'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(LaplaceMechanism(c(0,1),1,1,2,
                                which.sensitivity='both'),
               error=function(e) print(paste("PASS --",e)));
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

  print("         Only bounded if bounded==unbounded:")
  res = LaplaceMechanism(0,1,1,1,'both');
  if (is.null(res$Unbounded)) {print("PASS")}
  else{print("FAIL")};
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded example:");
  tv = 0;
  eps = 1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  ap = NULL;
  th.s = bs/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Simple bounded");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple unbounded example:");
  tv = 0;
  eps = 1;
  bs = NULL;
  us = 1;
  ws = 'unbounded';
  ap = NULL;
  th.s = us/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Unbounded;
  }
  hist(data,freq=FALSE,main="Simple unbounded");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple both example:");
  tv = 0;
  eps = 1;
  bs = 2;
  us = 1;
  ws = 'both';
  ap = NULL;
  th.s.b = bs/eps;
  th.s.u = us/eps;
  n = 10000;
  data.b = numeric(n);
  data.u = numeric(n);
  for (i in 1:n){
    data.b[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
    data.u[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Unbounded;
  }

  hist(data.b,freq=FALSE,main="Simple both bounded");
  x = seq(tv-5*th.s.b, tv+5*th.s.b,.1);
  lines(x,dlaplace(x,m=tv,s=th.s.b));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  hist(data.u,freq=FALSE,main="Simple both unbounded");
  x = seq(tv-5*th.s.u, tv+5*th.s.u,.1);
  lines(x,dlaplace(x,m=tv,s=th.s.u));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values:");
  tv = 5;
  eps = 1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  ap = NULL;
  th.s = bs/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Changing true.values");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values/epsilon:");
  tv = -5;
  eps = .1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  ap = NULL;
  th.s = bs/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Changing true.values/epsilon");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing epsilon/bounded.sensitivity:");
  tv = -5;
  eps = .1;
  bs = .5;
  us = NULL;
  ws = 'bounded';
  ap = NULL;
  th.s = bs/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Changing epsilon/bounded.sensitivity");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using unbounded:");
  tv = 10;
  eps = .5;
  bs = NULL;
  us = .2;
  ws = 'unbounded';
  ap = NULL;
  th.s = us/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Unbounded;
  }
  hist(data,freq=FALSE,main="Using unbounded");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dlaplace(x,m=tv,s=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using multiple:");
  tv = c(-5,0,10);
  eps = 1;
  bs = c(.5,1,2);
  us = NULL;
  ws = 'bounded';
  ap = NULL;
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
  }
  th.s = length(tv)*bs/eps;
  for (j in 1:length(tv)){
    hist(data[,j],freq=FALSE,main=paste("Using multiple:",j));
    x = seq(tv[j]-5*th.s[j], tv[j]+5*th.s[j],.1);
    lines(x,dlaplace(x,m=tv[j],s=th.s[j]));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Using allocation:");
  tv = c(-5,0);
  eps = 1;
  bs = c(.5,1);
  us = NULL;
  ws = 'bounded';
  ap = c(1,9);
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = LaplaceMechanism(tv,eps,bs,us,ws,ap)$Bounded;
  }
  ap = ap/sum(ap);
  th.s = bs/(ap*eps);
  for (j in 1:length(tv)){
    hist(data[,j],freq=FALSE,main=paste("Using allocation:",j));
    x = seq(tv[j]-5*th.s[j], tv[j]+5*th.s[j],.1);
    lines(x,dlaplace(x,m=tv[j],s=th.s[j]));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }
  ### END TEST FUNCTIONALITY ###
}

test_Gaussian <- function(){
  #
  #
  #
  source(paste(cur.dir,'Mechanisms.R',sep=""))

  ### Test input ###
  print("         No input:")
  a = tryCatch(GaussianMechanism(),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad true.values:");
  a = tryCatch(GaussianMechanism('a',1,1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad eps:")
  a = tryCatch(GaussianMechanism(0,'a',1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,-1,1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,0,1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad delta:")
  a = tryCatch(GaussianMechanism(0,1,'a',1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,-1,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,0,1),error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad bounded.sensitivities:")
  a = tryCatch(GaussianMechanism(0,1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(c(0,1),1,1,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,1,bounded.sensitivities=-1),
               error=function(e) print(paste("PASS --",e)));
  print("")

  print("         Bad unbounded.sensitivities:")
  a = tryCatch(GaussianMechanism(0,1,1,which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,unbounded.sensitivities=c(1,2),
                                which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(c(0,1),1,1,unbounded.sensitivities=1,
                                which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,unbounded.sensitivities=-1,
                                which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(GaussianMechanism(0,1,1,which.sensitivity='unbnd'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,unbounded.sensitivities=c(1,2),
                                which.sensitivity='bounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(c(0,1),1,1,unbounded.sensitivities=1,
                                which.sensitivity='both'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(c(0,1),1,1,1,2,
                                which.sensitivity='both'),
               error=function(e) print(paste("PASS --",e)));
  print("")

  print("         Bad type.DP:")
  a = tryCatch(GaussianMechanism(0,1,1,1,type.DP='pdp'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad alloc.proportions:")
  a = tryCatch(GaussianMechanism(0,1,1,1,alloc.proportions='a'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(0,1,1,1,alloc.proportions=-1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(GaussianMechanism(c(0,1),1,1,1,alloc.proportions=c(1,2,3)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Only bounded if bounded==unbounded:")
  res = GaussianMechanism(0,1,1,1,1,'both');
  if (is.null(res$Unbounded)) {print("PASS")}
  else{print("FAIL")};
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded example (pDP):");
  tv = 0;
  eps = .5;
  delta = .5;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  dp = 'pDP'
  ap = NULL;
  th.s = bs*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Simple bounded (pDP)");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple bounded example (aDP):");
  tv = 0;
  eps = .5;
  delta = .5;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  dp = 'aDP'
  ap = NULL;
  th.s = bs*sqrt(2*log(1.25/delta))/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Simple bounded (aDP)");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple unbounded example (pDP):");
  tv = 0;
  eps = .5;
  delta = .5;
  bs = NULL;
  us = 1;
  ws = 'unbounded';
  dp = 'pDP'
  ap = NULL;
  th.s = us*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Unbounded;
  }
  hist(data,freq=FALSE,main="Simple unbounded (pDP)");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple unbounded example (aDP):");
  tv = 0;
  eps = .5;
  delta = .5;
  bs = NULL;
  us = 1;
  ws = 'unbounded';
  dp = 'aDP'
  ap = NULL;
  th.s = us*sqrt(2*log(1.25/delta))/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Unbounded;
  }
  hist(data,freq=FALSE,main="Simple unbounded (aDP)");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple both example (pDP):");
  tv = 0;
  eps = .5;
  delta = .5;
  bs = 2;
  us = 1;
  ws = 'both';
  dp = 'pDP';
  ap = NULL;
  th.s.b = bs*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  th.s.u = us*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  n = 10000;
  data.b = numeric(n);
  data.u = numeric(n);
  for (i in 1:n){
    data.b[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
    data.u[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Unbounded;
  }

  hist(data.b,freq=FALSE,main="Simple both bounded (pDP)");
  x = seq(tv-5*th.s.b, tv+5*th.s.b,.1);
  lines(x,dnorm(x,m=tv,sd=th.s.b));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  hist(data.u,freq=FALSE,main="Simple both unbounded (pDP)");
  x = seq(tv-5*th.s.u, tv+5*th.s.u,.1);
  lines(x,dnorm(x,m=tv,sd=th.s.u));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Simple both example (aDP):");
  tv = 0;
  eps = .5;
  delta = .5;
  bs = 2;
  us = 1;
  ws = 'both';
  dp = 'aDP';
  ap = NULL;
  th.s.b = bs*sqrt(2*log(1.25/delta))/eps;
  th.s.u = us*sqrt(2*log(1.25/delta))/eps;
  n = 10000;
  data.b = numeric(n);
  data.u = numeric(n);
  for (i in 1:n){
    data.b[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
    data.u[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Unbounded;
  }

  hist(data.b,freq=FALSE,main="Simple both bounded (aDP)");
  x = seq(tv-5*th.s.b, tv+5*th.s.b,.1);
  lines(x,dnorm(x,m=tv,sd=th.s.b));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  hist(data.u,freq=FALSE,main="Simple both unbounded (aDP)");
  x = seq(tv-5*th.s.u, tv+5*th.s.u,.1);
  lines(x,dnorm(x,m=tv,sd=th.s.u));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values:");
  tv = 5;
  eps = .5;
  delta = .5;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  dp = 'pDP';
  ap = NULL;
  th.s = bs*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Changing true.values");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing true.values/epsilon:");
  tv = -5;
  eps = .1;
  delta = .5;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  dp = 'aDP';
  ap = NULL;
  th.s = bs*sqrt(2*log(1.25/delta))/eps;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Changing true.values/epsilon");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Changing delta/bounded.sensitivity:");
  tv = -5;
  eps = .5;
  delta = .1
  bs = 3;
  us = NULL;
  ws = 'bounded';
  dp = 'pDP';
  ap = NULL;
  th.s = bs*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  hist(data,freq=FALSE,main="Changing delta/bounded.sensitivity");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using unbounded:");
  tv = 10;
  eps = .5;
  delta = .5
  bs = NULL;
  us = .2;
  ws = 'unbounded';
  dp = 'pDP';
  ap = NULL;
  th.s = us*(sqrt(qnorm(delta/2)^2+2*eps)-qnorm(delta/2))/(2*eps);
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Unbounded;
  }
  hist(data,freq=FALSE,main="Using unbounded");
  x = seq(tv-5*th.s, tv+5*th.s,.1);
  lines(x,dnorm(x,m=tv,sd=th.s));
  print("Verify line matches histogram...")
  Sys.sleep(2)

  print("         Using multiple:");
  tv = c(-5,0,10);
  eps = .5
  delta = .5;
  bs = c(.5,1,2);
  us = NULL;
  ws = 'bounded';
  dp = 'aDP';
  ap = NULL;
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  ap = c(1,1,1)/3;
  th.s = bs*sqrt(2*log(1.25/(ap*delta)))/(ap*eps);
  for (j in 1:length(tv)){
    hist(data[,j],freq=FALSE,main=paste("Using multiple:",j));
    x = seq(tv[j]-5*th.s[j], tv[j]+5*th.s[j],.1);
    lines(x,dnorm(x,m=tv[j],sd=th.s[j]));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }

  print("         Using allocation:");
  tv = c(-5,0);
  eps = .5;
  delta = .5;
  bs = c(.5,1);
  us = NULL;
  ws = 'bounded';
  dp = 'pDP';
  ap = c(1,9);
  n = 10000;
  data = matrix(NaN,nrow=n,ncol=length(tv));
  for (i in 1:n){
    data[i,] = GaussianMechanism(tv,eps,delta,bs,us,ws,dp,ap)$Bounded;
  }
  ap = ap/sum(ap);
  th.s = bs*(sqrt(qnorm((ap*delta)/2)^2+2*(ap*eps))-qnorm((ap*delta)/2))/(2*(ap*eps));
  for (j in 1:length(tv)){
    hist(data[,j],freq=FALSE,main=paste("Using allocation:",j));
    x = seq(tv[j]-5*th.s[j], tv[j]+5*th.s[j],.1);
    lines(x,dnorm(x,m=tv[j],sd=th.s[j]));
    print("Verify line matches histogram...")
    Sys.sleep(2)
  }
  ### END TEST FUNCTIONALITY ###
}

test_Exponential <- function(){
  #
  #
  #
  source(paste(cur.dir,'Mechanisms.R',sep=""))

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

  print("         Bad bounded.sensitivities:")
  a = tryCatch(ExponentialMechanism(0,1),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,c(1,2)),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,-1),
               error=function(e) print(paste("PASS --",e)));
  print("")

  print("         Bad unbounded.sensitivities:")
  a = tryCatch(ExponentialMechanism(0,1,1,which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,unbounded.sensitivities=c(1,2),
                                 which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,unbounded.sensitivities=-1,
                                 which.sensitivity='unbounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  print("")

  print("         Bad which.sensitivity:")
  a = tryCatch(ExponentialMechanism(0,1,which.sensitivity='unbnd'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,unbounded.sensitivities=1,
                                 which.sensitivity='bounded'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
  a = tryCatch(ExponentialMechanism(0,1,c(0,1),unbounded.sensitivities=1,
                                 which.sensitivity='both'),
               error=function(e) print(paste("PASS --",e)));
  if (!is.character(a)) print("FAIL");
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

  print("         Only bounded if bounded==unbounded:")
  res = ExponentialMechanism(0,1,1,1,'both');
  if (is.null(res$Unbounded)) {print("PASS")}
  else{print("FAIL")};
  print("")
  ### END TEST INPUT ###

  ### TEST FUNCTIONALITY ###
  print("         Simple bounded example:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Simple Bounded");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Simple unbounded example:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  bs = NULL;
  us = 1;
  ws = 'unbounded';
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Unbounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*us));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Simple Unbounded");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Simple both example:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  bs = 1;
  us = 2;
  ws = 'both';
  m = NULL;
  c = Z;
  n = 10000;
  data.b = numeric(n);
  data.u = numeric(n);
  for (i in 1:n){
    data.b[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
    data.u[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Unbounded;
  }

  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data.b)/sum(table(data.b)),main="Simple both bounded");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*us));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data.u)/sum(table(data.u)),main="Simple both unbounded");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Changing utility:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = 1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Different Utility");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Changing utility/epsilon:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = .1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Different utility/epsilon");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Changing epsilon/bounded.sensitivity:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = .1;
  bs = 5;
  us = NULL;
  ws = 'bounded';
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Different epsilon/bounded.sensitivity");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Using unbounded:");
  Z = c(0,1,2,3,4,5);
  u = c(-5,-1,0,3,-6,-2)
  eps = 1;
  bs = NULL;
  us = 1;
  ws = 'unbounded';
  m = NULL;
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Unbounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*us));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Using unbounded");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Using measure:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  m = c(.5,1,2,.5,5,1);
  c = Z;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(Z,table(data)/sum(table(data)),main="Using Measure");
  lines(Z,th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)

  print("         Returning Index:");
  Z = c(0,1,2,3,4,5);
  u = -abs((1:length(Z))-.5-.5*length(Z));
  eps = 1;
  bs = 1;
  us = NULL;
  ws = 'bounded';
  m = NULL;
  c = NULL;
  n = 10000;
  data = numeric(n);
  for (i in 1:n){
    data[i] = ExponentialMechanism(u,eps,bs,us,ws,m,c)$Bounded;
  }
  if (is.null(m)) m = rep(1,length(Z));
  th.probs = m*exp(eps*u/(2*bs));
  th.probs = th.probs/sum(th.probs);
  plot(1:length(Z),table(data)/sum(table(data)),main="Returning Index");
  lines(1:length(Z),th.probs);
  print("Verify line matches dots...")
  Sys.sleep(2)
  ### END TEST FUNCTIONALITY ###
}

