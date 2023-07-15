# library(rgl)
library(MASS)
library(parallel)
library(doParallel)
library(igraph)
numCores = detectCores()
# cl = makeCluster(numCores)
# clusterEvalQ(cl, {
#     library(MASS)
#     library(parallel)
#     library(doParallel)
#     library(igraph)
# })
registerDoParallel(numCores)

set.seed(1729)

runs_test_old=function(m=100,n=100,p=2,exponent=2,repl=100,num_mst=3,alpha=0.05,del=0,sig=1,distr="normal",pval=F)
{
    N=m+n
    result <<- numeric(num_mst)
    p_val <<- numeric(num_mst)
    a0=2*m*n/(N*(N-1))
    a1=2*(m-1)*(n-1)/((N-2)*(N-3))
    for (reps in 1:repl)
    {
        if(distr=='normal')dat=matrix(c(rnorm(p*m),rnorm(p*n,del/sqrt(p),sig)),p)
        else if(distr=="lognormal")dat=matrix(c(rlnorm(p*m),rlnorm(p*n,del/sqrt(p),sig)),p)
        
        # dat <<- dat
        
        mat <<- as.matrix(dist(t(dat),diag=TRUE,upper=TRUE))
        mat[upper.tri(mat,diag=TRUE)] <<- Inf
        
        finder=function(num0)
        {
            while (pnt[num0] !=num0) {num0=pnt[num0]}
            return(num0)      
        }
        merger=function(root1,root2) 
        {
            if (sz[root1]<sz[root2]) {pnt[root1] <<- root2; sz[root2] <<- sz[root1]+sz[root2]}
            else {pnt[root2] <<- root1; sz[root1] <<- sz[root1]+sz[root2]}
        }
        
        msts <<- vector(mode="list",length=num_mst)
        orth_mst=function(num)
        {
            z=sort.list(mat)[1:((N*(N-1)/2)-(num-1)*(N-1))]
            mst=matrix(0,(N-1),2)    #Edges of the minimal spanning tree stored here in rows
            pnt <<- c(1:N)
            sz <<- rep(1,N)
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
                stat_val=(R0-(a0*e0))/sqrt(b0)
                p_val[s] <<- p_val[s]+pnorm(stat_val)
                if (stat_val<qnorm(alpha)) {result[s] <<- result[s]+1}
            }
        }
        tester(msts)
    }
    #cat("\nNo. of rejections using orthogonal mst nos.",c(1:num_mst),"are",result,"respectively.\n")
    # print(result)
    if (pval==TRUE) {return(p_val)}
    else {return(result)}
}
runs_test=function(m=100,n=100,p=2,exponent=2,repl=100,num_mst=3,alpha=0.05,del=0,sig=1,distr="normal",pval=F)
{
    N=m+n
    result<-numeric(num_mst)
    p_val<- numeric(num_mst)
    a0=2*m*n/(N*(N-1))
    a1=2*(m-1)*(n-1)/((N-2)*(N-3))
    for (reps in 1:repl){
    # loop=function(reps){
    # foreach(reps=1:repl,.packages = 'igraph')%dopar%{
        if(distr=='normal')dat=matrix(c(rnorm(p*m),rnorm(p*n,del/sqrt(p),sig)),p)
        else if(distr=="lognormal")dat=matrix(c(rlnorm(p*m),rlnorm(p*n,del/sqrt(p),sig)),p)
        mat<- as.matrix(dist(t(dat),diag=TRUE,upper=TRUE))
        # mat[upper.tri(mat,diag=TRUE)] <- Inf
        G<-graph_from_adjacency_matrix(mat,weighted=TRUE,mode="undirected")
        msts <- vector(mode="list",length=num_mst)
        for (i in 1:num_mst){
            MST<-mst(G)
            msts[[i]]<-as_edgelist(MST,F)
            G=G-MST
        }
        # print(msts)
        tester=function(mstrees)
        {
            for (s in 1:num_mst)
            {
                # temp6=rbind(mstrees[1:s])
                temp6=c()
                for (s0 in 1:s) {temp6=rbind(temp6,mstrees[[s0]])}
                R0<-sum(apply(temp6,MARGIN=1,function(row)((row[1]>m&&row[2]<m+1)||(row[2]>m&&row[1]<m+1))))
                # adj=get.adjacency(graph_from_edgelist(temp6,directed=F))
                # adj2=as.matrix(adj%*%adj)
                # c0=sum(adj2[upper.tri(adj2)])
                deg<- degree(graph_from_edgelist(temp6,directed=F))
                c0=sum(deg*(deg-1)/2)
                c0=c0-nrow(temp6)
                e0=s*(N-1)
                b0=a0*(e0+c0+a1*(e0*(e0-1)-2*c0)-a0*e0*e0)
                stat_val=(R0-(a0*e0))/sqrt(b0)
                p_val[s] <<- p_val[s]+pnorm(stat_val)
                if (stat_val<qnorm(alpha)) {result[s]<<-result[s]+1}
            }
        }
        tester(msts)
    }
    # parSapply(cl,1:repl,loop)
    # mclapply(1:repl,loop)
    # lapply(1:repl,loop)
    # stopCluster(cl)
    #cat("\nNo. of rejections using orthogonal mst nos.",c(1:num_mst),"are",result,"respectively.\n")
    print(result)
    if (pval==TRUE) {return(p_val)}
    else {return(result)}
}
runs_test_newest=function(m=100,n=100,p=2,exponent=2,repl=3,num_mst=3,alpha=0.05,del=0,sig=1,distr="normal",pval=F)
{
    N=m+n
    # p_val<- numeric(num_mst)
    a0=2*m*n/(N*(N-1))
    a1=2*(m-1)*(n-1)/((N-2)*(N-3))
    # for (reps in 1:repl){
        # loop=function(reps){
    foreach(reps=1:repl,.combine=rbind,.packages='igraph')%dopar%{
        result<- numeric(num_mst)
        if(distr=='normal')dat=matrix(c(rnorm(p*m),rnorm(p*n,del/sqrt(p),sig)),p)
        else if(distr=="lognormal")dat=matrix(c(rlnorm(p*m),rlnorm(p*n,del/sqrt(p),sig)),p)
        mat<- as.matrix(dist(t(dat),diag=TRUE,upper=TRUE))
        G<-graph_from_adjacency_matrix(mat,weighted=TRUE,mode="undirected")
        msts <- vector(mode="list",length=num_mst)
        for (i in 1:num_mst){
            MST<-mst(G)
            msts[[i]]<-as_edgelist(MST,F)
            G=G-MST
            temp6=c()
            for (s0 in 1:i) {temp6=rbind(temp6,msts[[s0]])}
            R0<-sum(apply(temp6,MARGIN=1,function(row)((row[1]>m&&row[2]<m+1)||(row[2]>m&&row[1]<m+1))))
            # adj=get.adjacency(graph_from_edgelist(temp6,directed=F))
            # adj2=as.matrix(adj%*%adj)
            # c0=sum(adj2[upper.tri(adj2)])
            deg<- degree(graph_from_edgelist(temp6,directed=F))
            c0=sum(deg*(deg-1)/2)
            # c0=c0-nrow(temp6)
            e0=i*(N-1)
            b0=a0*(e0+c0+a1*(e0*(e0-1)-2*c0)-a0*e0*e0)
            stat_val=(R0-(a0*e0))/sqrt(b0)
            # p_val[i] <- p_val[i]+pnorm(stat_val)
            if (stat_val<qnorm(alpha)) {result[i]<-result[i]+1}
        }
        result
    }->out
    result=colSums(out)
    # parSapply(cl,1:repl,loop)
    # mclapply(1:repl,loop)
    # lapply(1:repl,loop)
    # stopCluster(cl)
    #cat("\nNo. of rejections using orthogonal mst nos.",c(1:num_mst),"are",result,"respectively.\n")
    # print(result)
    if (pval==TRUE) {return(p_val)}
    else {return(result)}
}

deg_test=function(m=100,n=100,p,exponent=2,repl,alpha=0.05,del=0,sig=1,distr,pval)
{
    N=m+n
    cutoff <<- numeric(repl)
    p_val <<- 0
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
            while (pnt[num0] !=num0) {num0=pnt[num0]}
            return(num0)      
        }
        merger=function(root1,root2) 
        {
            if (sz[root1]<sz[root2]) {pnt[root1] <<- root2; sz[root2] <<- sz[root1]+sz[root2]}
            else {pnt[root2] <<- root1; sz[root1] <<- sz[root1]+sz[root2]}
        }
        
        msts <<- vector(mode="list",length=1)
        orth_mst=function(num)
        {
            z=sort.list(mat)[1:((N*(N-1)/2)-(num-1)*(N-1))]
            mst=matrix(0,(N-1),2)    #Edges of the minimal spanning tree stored here in rows
            pnt <<- c(1:N)
            sz <<- rep(1,N)
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
        p_val <<- p_val+phyper(x1,d1,(N-d1),m)
        if (x1<cutoff[reps]) {result=result+1}
        if (x1==(cutoff[reps])) {result=result+(runif(1)<((alpha-phyper((cutoff[reps]-1),d1,(N-d1),m))/dhyper(cutoff[reps],d1,(N-d1),m)))}
    }
    #cat("\nNo. of rejections using degree 1 vertices is",result,"\n")
    if (pval==TRUE) {return(p_val)}
    else {return(result)}
}

normal_test=function(m=100,n=100,p,exponent=2,repl,alpha=0.05,del=0,sig=1,distr,pval=F)
{
    N=m+n
    test_stat <<- numeric(repl)
    p_val <<- 0
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
        p_val <<- p_val+1-pchisq(test_stat[reps],(p*(p+3)/2))
        if (test_stat[reps]>qchisq((1-alpha),(p*(p+3)/2))) {result=result+1}
    }
    #cat("\nNo. of rejections using normal theory is",result,"\n")
    if (pval==TRUE) {return(p_val)}
    else {return(result)}
}

plots_old=function(dims0,par0,trials0,param,testname='normal',distrib,pvals=F)
{
    a1 <<- seq(dims0[1],dims0[2],by=dims0[3])
    a10=length(a1)
    b1 <<- (a1-a1[1])/(max(a1)-min(a1))
    a2 <<- seq(par0[1],par0[2],by=par0[3])
    a20=length(a2)
    b2 <<- (a2-a2[1])/(max(a2)-min(a2))
    ind="rejfrac"
    ind0="Rejection Fraction"
    if (pvals==TRUE) 
    {
        ind="pval"
        ind0="Average p-value"
    }
    v=1
    if (testname=="runs") {v=3}
    val=array(0,dim=c(a10,a20,v))
    rownames0=numeric(a10)
    colnames0=numeric(a20)
    for (i in 1:a10)
    {
        rownames0[i]=paste("Dim",a1[i],sep=" ")
    }
    
    for (i in 1:a20)
    {
        colnames0[i]=paste(param,"=",a2[i],sep="")
    }  
    
    for (i in 1:a10)
    {
        print(i)
        for (j in 1:a20)
        {
            tempname0=paste(testname,"_","test",sep="")
            tempname1=list()
            tempname1[["p"]]=a1[i]
            tempname1[[param]]=a2[j]
            tempname1[["repl"]]=trials0
            tempname1[["distr"]]=distrib
            tempname1[["pval"]]=pvals
            val[i,j,]=do.call(eval(as.symbol(tempname0)),tempname1)
        } 
    }
    val=val/trials0
    for (i in 1:v)
    {
        tempp=round(val[,,i],digits=6)
        rownames(tempp)=rownames0
        colnames(tempp)=colnames0
        # persp3d(b1,b2,tempp,col="green",xlab="Dimensions",ylab=param,zlab=paste(ind0," (",testname," test",")",sep=""))
        # filename1=paste("C:/Users/User/Desktop/NONPARA PLOTS/",distrib,"_",param,"_",testname,"_",i,"_",ind,".stl",sep="")
        # writeSTL(filename1)
        # persp3d(a1,a2,tempp,col="green",xlab="Dimensions",ylab=param,zlab=paste(ind0," (",testname," test",")",sep=""))
        filename0=paste("./",testname,"_",distrib,"_",param,"_mst",i,".csv",sep="")
        write.csv(tempp,filename0)
    }
    val <- val
}
plots=function(dims0,par0,trials0,param,testname,distrib,pvals=F)
{
    a1 <- seq(dims0[1],dims0[2],by=dims0[3])
    a10=length(a1)
    b1 <- (a1-a1[1])/(max(a1)-min(a1))
    a2 <- seq(par0[1],par0[2],by=par0[3])
    a20=length(a2)
    b2 <- (a2-a2[1])/(max(a2)-min(a2))
    ind="rejfrac"
    ind0="Rejection Fraction"
    if (pvals==TRUE){
        ind="pval"
        ind0="Average p-value"
    }
    v=1
    if (testname=="runs") {v=3}
    val=array(0,dim=c(a10,a20,v))
    rownames0=paste("Dim",a1,sep=" ")
    colnames0=paste(param,"=",a2,sep="")
    for (i in 1:a10)
    {
        print(i)
        foreach(j=1:a20,.combine=rbind,.packages=c('igraph','parallel','doParallel'),.export = c('a1','a2','runs_test'))%dopar%{
            tempname0=paste(testname,"_","test",sep="")
            tempname1=list()
            tempname1[["p"]]=a1[i]
            tempname1[[param]]=a2[j]
            tempname1[["repl"]]=trials0
            tempname1[["distr"]]=distrib
            tempname1[["pval"]]=pvals
            set.seed(1729)
            do.call(eval(as.symbol(tempname0)),tempname1)
        }->out 
        val[i,,]=out
    }
    val=val/trials0
    for (i in 1:v)
    {
        tempp=round(val[,,i],digits=6)
        rownames(tempp)=rownames0
        colnames(tempp)=colnames0
        # persp3d(b1,b2,tempp,col="green",xlab="Dimensions",ylab=param,zlab=paste(ind0," (",testname," test",")",sep=""))
        # filename1=paste("C:/Users/User/Desktop/NONPARA PLOTS/",distrib,"_",param,"_",testname,"_",i,"_",ind,".stl",sep="")
        # writeSTL(filename1)
        # persp3d(a1,a2,tempp,col="green",xlab="Dimensions",ylab=param,zlab=paste(ind0," (",testname," test",")",sep=""))
        filename0=paste("./Outputs/",testname,"_",distrib,"_",param,"_mst",i,".csv",sep="")
        write.csv(tempp,filename0)
    }
    val<- val
}

# a=commandArgs(1)
# #Normal theory test
# 
# set.seed(1729)
# if(a==1)system.time({plots(c(1,20,1),c(0,1.3,0.01),1000,"del","normal","normal")})    #This was done
# #plots(c(1,20,1),c(0,1.3,0.05),10,"del","normal","normal")      #Bachha version
# 
# set.seed(1729)
# if(a==2)system.time({plots(c(1,20,1),c(1,1.4,0.005),1000,"sig","normal","normal")})   #This was done
# #plots(c(1,20,1),c(1,1.4,0.02),10,"sig","normal","normal")      #Bachha version
# 
# set.seed(1729)
# if(a==3)system.time({plots(c(1,20,1),c(0,0.5,0.01),1000,"del","normal","lognormal")}) #Needs to be done
# #plots(c(1,10,1),c(0,0.5,0.05),10,"del","normal","lognormal")   #Bachha version
# 
# 
# #Degree test
# 
# set.seed(1729)
# if(a==4)system.time({plots(c(1,20,1),c(0,1.3,0.01),100,"del","deg","normal")})   #Needs to be done
# #plots(c(1,20,1),c(0,1.3,0.05),1,"del","deg","normal")          #Bachha version
# 
# set.seed(1729)
# if(a==5)system.time({plots(c(1,20,1),c(1,1.4,0.1),100,"sig","deg","normal")})      #Needs to be done
# #plots(c(1,20,1),c(1,1.4,0.05),1,"sig","deg","normal")          #Bachha version
# 
# set.seed(1729)
# if(a==6)system.time({plots(c(1,20,1),c(0,0.5,0.01),100,"del","deg","lognormal")})     #Needs to be done
# #plots(c(1,20,1),c(0,0.5,0.05),1,"del","deg","lognormal")       #Bachha version


#Runs Test
    
set.seed(1729)
system.time({plots(c(1,20,1),c(0,1.3,0.01),1000,"del","runs","normal")})        #Needs to be done
# plots(c(1,20,1),c(0,1.3,0.05),1,"del","runs","normal")          #Bachha version

set.seed(1729)
system.time({plots(c(1,20,1),c(1,1.4,0.005),1000,"sig","runs","normal")})       #Needs to be done
#plots(c(1,20,1),c(1,1.4,0.05),1,"sig","runs","normal")          #Bachha version

set.seed(1729)
system.time({plots(c(1,20,1),c(0,0.5,0.01),1000,"del","runs","lognormal")})     #Needs to be done
#plots(c(1,20,1),c(0,0.5,0.05),1,"del","runs","lognormal")       #Bachha version


