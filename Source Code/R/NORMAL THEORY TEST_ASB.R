set.seed(1729)

normal_test=function(m=100,n=100,p,exponent=2,repl,alpha=0.05,del=0,sig=1)
{
  N=m+n
  test_stat <<- numeric(repl)
  result=0
  for (reps in 1:repl)
  {
    dat=matrix(0,p,(m+n))
    for (i in 1:p)
    {
      dat[i,1:m]=rnorm(m)
      dat[i,(m+1):(m+n)]=rnorm(n,del/sqrt(p),sig)
    }

    dat <<- dat
    datt <<- t(dat)

    if (p==1)
    {
      mu0=mean(dat)
      mu1=mean(dat[1:m])
      mu2=mean(dat[(m+1):(m+n)])
    }
    else
    {
      mu0=apply(dat,1,mean)
      mu1=apply(dat[,1:m],1,mean)
      mu2=apply(dat[,(m+1):(m+n)],1,mean)
    }
    var0=determinant(as.matrix(var(datt))*(N-1)/N)$mod[1]
    var1=determinant(as.matrix(var(datt[1:m,])*(m-1)/m))$mod[1]
    var2=determinant(as.matrix(var(datt[(m+1):(m+n),])*(n-1)/n))$mod[1]

    test_stat[reps] <<- N*var0-m*var1-n*var2

    if (test_stat[reps]>qchisq((1-alpha),(p*(p+3)/2))) {result=result+1}
  }
  cat("\nNo. of rejections using normal theory is",result,"\n")
}

normal_test(p=1,del=0.3,repl=10000)
normal_test(p=2,del=0.5,repl=10000)
normal_test(p=5,del=0.75,repl=10000)
normal_test(p=10,del=1,repl=10000)
normal_test(p=20,del=1.2,repl=10000)

normal_test(p=1,sig=1.3,repl=10000)
normal_test(p=2,sig=1.2,repl=10000)
normal_test(p=5,sig=1.2,repl=10000)
normal_test(p=10,sig=1.1,repl=10000)
normal_test(p=20,sig=1.075,repl=10000)
