set.seed(1729)

deg_test=function(m=100,n=100,p,exponent=2,repl,alpha=0.05,del=0,sig=1,distr="normal")
{
  N=m+n
  cutoff <<- numeric(repl)
  result=0
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

    msts <<- vector(mode="list",length=1)
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

    orth_mst(1)

    temp0=table(msts[[1]])
    leaves=as.numeric(rownames(temp0)[temp0==1])
    d1=length(leaves)
    x1=length(leaves[leaves<=m])
    #if (phyper(x1,d1,(N-d1),m)<=alpha) {result=result+1}
    cutoff[reps] <<- qhyper(alpha,d1,(N-d1),m)
    if (x1<cutoff[reps]) {result=result+1}
    if (x1==(cutoff[reps])) {result=result+(runif(1)<((alpha-phyper((cutoff[reps]-1),d1,(N-d1),m))/dhyper(cutoff[reps],d1,(N-d1),m)))}
  }
  cat("\nNo. of rejections using degree 1 vertices is",result,"\n")
}

deg_test(p=1,del=0.3,repl=100)
deg_test(p=2,del=0.5,repl=100)
deg_test(p=5,del=0.75,repl=100)
deg_test(p=10,del=1,repl=100)
deg_test(p=20,del=1.2,repl=100)

deg_test(p=1,sig=1.3,repl=100)
deg_test(p=2,sig=1.2,repl=100)
deg_test(p=5,sig=1.2,repl=100)
deg_test(p=10,sig=1.1,repl=100)
deg_test(p=20,sig=1.075,repl=100)

deg_test(p=1,del=0.4,repl=100,distr="lognormal")
deg_test(p=2,del=0.4,repl=100,distr="lognormal")
deg_test(p=5,del=0.3,repl=100,distr="lognormal")
deg_test(p=10,del=0.3,repl=100,distr="lognormal")
deg_test(p=20,del=0.3,repl=100,distr="lognormal")

