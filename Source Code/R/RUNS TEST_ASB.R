#m,n sample sizes of two p dimensional distributions, exponent is used for
#calculating Euclidean distance, repl is the number of replications of the expt.,
#alpha is the prob. of type 1 error, del=location change, sig=scale change,
#num_mst is the number of orthogonal msts to be used

set.seed(1729)

runs_test=function(m=100,n=100,p,exponent=2,repl,num_mst=3,alpha=0.05,del=0,sig=1,distr="normal")
{
  N=m+n
  result <<- numeric(num_mst)
  a0=2*m*n/(N*(N-1))
  a1=2*(m-1)*(n-1)/((N-2)*(N-3))
  for (reps in 1:repl)
  {
    dat=matrix(0,p,(m+n))
    for (i in 1:p)
    {
      dat[i,1:m]=rnorm(m)
      dat[i,(m+1):(m+n)]=rnorm(n,del/sqrt(p),sig)
    }

    if (distr=="lognormal")
    for (i in 1:p)
    {
      dat[i,1:m]=rlnorm(m)
      dat[i,(m+1):(m+n)]=rlnorm(n,meanlog=del)
    }

    dat <<- dat

    mat <<- as.matrix(dist(t(dat),diag=TRUE,upper=TRUE))
    mat[upper.tri(mat,diag=TRUE)] <<- Inf

    finder=function(num0)
    {
      while (parent[num0] !=num0) {num0=parent[num0]}
      return(num0)      
    }
    merger=function(root1,root2) 
    {
      if (sizes[root1]<sizes[root2]) {parent[root1] <<- root2; sizes[root2] <<- sizes[root1]+sizes[root2]}
      else {parent[root2] <<- root1; sizes[root1] <<- sizes[root1]+sizes[root2]}
    }

    msts <<- vector(mode="list",length=num_mst)
    orth_mst=function(num)
    {
      z=sort.list(mat)[1:((N*(N-1)/2)-(num-1)*(N-1))]
      mst=matrix(0,(N-1),2)    #Edges of the minimal spanning tree stored here in rows
      parent <<- c(1:N)
      sizes <<- rep(1,N)
      count=1
      for (i in 1:length(z))
      {
        temp1=1+(z[i]+N-1)%%N  #row
        temp2=ceiling(z[i]/N)  #column
        temp3=finder(temp1)
        temp4=finder(temp2)
        if (temp3 !=temp4)
        {
          merger(temp3,temp4)
          mst[count,]=c(temp1,temp2)
          mat[temp1,temp2] <<- Inf
          count=count+1
        }
      }
      msts[[num]] <<- mst
    }

    for (i in 1:num_mst)
    {
      orth_mst(i)
    }

    tester=function(mstrees)
    {
      for (s in 1:num_mst)
      {
        temp6=c()
        for (s0 in 1:s) {temp6=rbind(temp6,mstrees[[s0]])}
        R0=0
        for (i in 1:nrow(temp6)) 
        {
          if ((temp6[i,1]>m) && (temp6[i,2]<(m+1))) {R0=R0+1}
        }
        c0=0
        for (i in 1:nrow(temp6))
        {
          for (j in i:nrow(temp6))
          {
            if ((temp6[i,1]==temp6[j,1]) || (temp6[i,1]==temp6[j,2]) || (temp6[i,2]==temp6[j,1]) || (temp6[i,2]==temp6[j,2])) {c0=c0+1}
          }
        }
        c0=c0-nrow(temp6)
        e0=s*(N-1)
        b0=a0*(e0+c0+a1*(e0*(e0-1)-2*c0)-a0*e0*e0)
        if (((R0-(a0*e0))/sqrt(b0))<qnorm(alpha)) {result[s] <<- result[s]+1}
      }
    }
    tester(msts)
  }
  cat("\nNo. of rejections using orthogonal mst nos.",c(1:num_mst),"are",result,"respectively.\n")
}

runs_test(p=1,del=0.3,repl=100)
runs_test(p=2,del=0.5,repl=100)
runs_test(p=5,del=0.75,repl=100)
runs_test(p=10,del=1,repl=100)
runs_test(p=20,del=1.2,repl=100)

runs_test(p=1,sig=1.3,repl=100)
runs_test(p=2,sig=1.2,repl=100)
runs_test(p=5,sig=1.2,repl=100)
runs_test(p=10,sig=1.1,repl=100)
runs_test(p=20,sig=1.075,repl=100)

runs_test(p=1,del=0.4,repl=100,distr="lognormal")
runs_test(p=2,del=0.4,repl=100,distr="lognormal")
runs_test(p=5,del=0.3,repl=100,distr="lognormal")
runs_test(p=10,del=0.3,repl=100,distr="lognormal")
runs_test(p=20,del=0.3,repl=100,distr="lognormal")

