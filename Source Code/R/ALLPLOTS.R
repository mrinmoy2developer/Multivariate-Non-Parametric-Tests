library(reticulate)
library(rgl)
library(profvis)

set.seed(1729)

runs_test=function(m=100,n=100,p,exponent=2,repl,num_mst=3,alpha=0.05,del=0,sig=1,distr)
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
  #cat("\nNo. of rejections using orthogonal mst nos.",c(1:num_mst),"are",result,"respectively.\n")
  return(result)
}

deg_test=function(m=100,n=100,p,exponent=2,repl,alpha=0.05,del=0,sig=1,distr)
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
  #cat("\nNo. of rejections using degree 1 vertices is",result,"\n")
  return(result)
}

normal_test=function(m=100,n=100,p,exponent=2,repl,alpha=0.05,del=0,sig=1,distr)
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

    if (distr=="lognormal")
    for (i in 1:p)
    {
      dat[i,1:m]=rlnorm(m)
      dat[i,(m+1):(m+n)]=rlnorm(n,meanlog=del)
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
  #cat("\nNo. of rejections using normal theory is",result,"\n")
  return(result)
}

plots=function(dims0,par0,trials0,param,testname,distrib)
{
  a1 <<- seq(dims0[1],dims0[2],by=dims0[3])
  b1 <<- (a1-a1[1])/(max(a1)-min(a1))
  a2 <<- seq(par0[1],par0[2],by=par0[3])
  b2 <<- (a2-a2[1])/(max(a2)-min(a2))
  val=matrix(0,length(a1),length(a2))

  rownames0=numeric(length(a1))
  colnames0=numeric(length(a2))
  for (i in 1:length(a1))
  {
    rownames0[i]=paste("Dim",a1[i],sep=" ")
  }

  for (i in 1:length(a2))
  {
    colnames0[i]=paste(param,"=",a2[i],sep="")
  }  

  rownames(val)=rownames0
  colnames(val)=colnames0

  for (i in 1:nrow(val))
  {
print(i)
    for (j in 1:ncol(val))
    {
      tempname0=paste(testname,"_","test",sep="")
      tempname1=list()
      tempname1[["p"]]=a1[i]
      tempname1[[param]]=a2[j]
      tempname1[["repl"]]=trials0
      tempname1[["distr"]]=distrib
      output=do.call(eval(as.symbol(tempname0)),tempname1)
      val[i,j]=output[1]
    } 
  }
  val=val/trials0
  #surface3d(b1,b2,val,col="green")
  persp3d(b1,b2,val,col="green",xlab="Dimensions",ylab=param,zlab=paste("Rejection Fraction ","(",testname," test",")",sep=""))
  filename1=paste("./",testname,"_",distrib,"_",param,".stl",sep="")
  writeSTL(filename1)
  persp3d(a1,a2,val,col="green",xlab="Dimensions",ylab=param,zlab=paste("Rejection Fraction ","(",testname," test",")",sep=""))
  filename0=paste("./",testname,"_",distrib,"_",param,".csv",sep="")
  write.csv(val,filename0)
  val <<- val
}

a=commandArgs(1)
#Normal theory test

set.seed(1729)
if(a==1)system.time({plots(c(1,20,1),c(0,1.3,0.01),1000,"del","normal","normal")})    #This was done
#plots(c(1,20,1),c(0,1.3,0.05),10,"del","normal","normal")      #Bachha version

set.seed(1729)
if(a==2)system.time({plots(c(1,20,1),c(1,1.4,0.005),1000,"sig","normal","normal")})   #This was done
#plots(c(1,20,1),c(1,1.4,0.02),10,"sig","normal","normal")      #Bachha version

set.seed(1729)
if(a==3)system.time({plots(c(1,20,1),c(0,0.5,0.01),1000,"del","normal","lognormal")}) #Needs to be done
#plots(c(1,10,1),c(0,0.5,0.05),10,"del","normal","lognormal")   #Bachha version


#Degree test

set.seed(1729)
if(a==4)system.time({plots(c(1,20,1),c(0,1.3,0.01),100,"del","deg","normal")})   #Needs to be done
#plots(c(1,20,1),c(0,1.3,0.05),1,"del","deg","normal")          #Bachha version

set.seed(1729)
if(a==5)system.time({plots(c(1,20,1),c(1,1.4,0.1),100,"sig","deg","normal")})      #Needs to be done
#plots(c(1,20,1),c(1,1.4,0.05),1,"sig","deg","normal")          #Bachha version

set.seed(1729)
if(a==6)system.time({plots(c(1,20,1),c(0,0.5,0.01),100,"del","deg","lognormal")})     #Needs to be done
#plots(c(1,20,1),c(0,0.5,0.05),1,"del","deg","lognormal")       #Bachha version


#Runs Test

set.seed(1729)
if(a==7)system.time({plots(c(1,20,1),c(0,1.3,0.01),100,"del","runs","normal")})        #Needs to be done
# plots(c(1,20,1),c(0,1.3,0.05),1,"del","runs","normal")          #Bachha version

set.seed(1729)
if(a==8)system.time({plots(c(1,20,1),c(1,1.4,0.005),100,"sig","runs","normal")})       #Needs to be done
#plots(c(1,20,1),c(1,1.4,0.05),1,"sig","runs","normal")          #Bachha version

set.seed(1729)
if(a==9)system.time({plots(c(1,20,1),c(0,0.5,0.01),100,"del","runs","lognormal")})     #Needs to be done
#plots(c(1,20,1),c(0,0.5,0.05),1,"del","runs","lognormal")       #Bachha version

