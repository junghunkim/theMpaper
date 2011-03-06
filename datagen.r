rm(list=ls())

source("mySEEDer.R")
source("myFILEer.R")

maxTime = 1
nstep = 100
nvertex = 12

myunitRate = 10

myPsi.C1 = rep(c(0.25,-4,0.8),5)
myPsi.C2 = rep(c(0.25,-4,0.3),7)
myPsi = c(myPsi.C1,myPsi.C2)
myPsi = matrix(myPsi,byrow=T,nrow=nvertex)

myinit = as.vector(myPsi[,3])
myIDmat = diag(nvertex)

mysimulator <- function(psi=NULL, init=NULL, LHS=NULL, RHS=NULL){       

  sigma <- function(y) {
    1
  }

  mu <- function(y) {
    cur.vol = psi[1]
    cur.beta = psi[2]
    cur.pi = psi[3]
    
    if( y*cur.vol > pi/2 ) {
      x = 0.999
    } else if (y*cur.vol < -pi/2) {
      x = 0.001
    } else {
      x = (sin(cur.vol*y)+1)/2 
    }
    
    L = (cur.beta/cur.vol)*(x-cur.pi)
    R = (-1/4)*cur.vol*(1-2*x)
    D = sqrt(x-x^2)
    
    (L+R)/D
  }
  
  dT = c(0,rep((RHS-LHS)/nstep,nstep))
  time.grid = LHS + cumsum(dT)
  Y = NULL;Y[1] = asin(2*init-1)/psi[1];
  dY = NULL;dY[1] = 0;

  for (i in 2:(length(time.grid))) {
    drift = mu(Y[i-1])*dT[i]
    diffs = rnorm(1)*sqrt(dT[i])*sigma(Y[i-1])
    
    dY[i] = drift + diffs
    Y[i] = Y[i-1] + dY[i]   
  }
  
  X = (sin(psi[1]*Y)+1)/2
  
  retObj = list(LHS=LHS, RHS=RHS, start=init, end=X[[length(time.grid)]], ts.time=time.grid, ts.state=X)  
}

ts.hyper.latent.ij = vector("list",nvertex)
for(i in 1:nvertex) {
  ts.hyper.latent.ij[[i]] <- mysimulator(myPsi[i,],init=myinit[i],LHS=0,RHS=maxTime)
}

ts.accept.ij = NULL;
ts.itr = 1;
ts.obs.ij = NULL;
ts.latent.ij = NULL;
ts.TKVs = NULL;

for( i in 1:(nvertex-1) ) {
  for( j in (i+1):nvertex) {
#  cat("************ Current (i,j) iteration : (", i,",",j,")\n")   
    myrate = maxTime * myunitRate
    n.ij = rpois(1,lambda= myrate)
    times = sort(runif(n.ij))*maxTime
    marks = vector("numeric",length(times))
    
    if(length(times) > 0 ) {
        for(k in 1:length(times)) {
          t.cur.i = min(which(ts.hyper.latent.ij[[i]]$ts.time > times[k]))
          t.cur.j = min(which(ts.hyper.latent.ij[[j]]$ts.time > times[k]))
    
          p.cur.i = ts.hyper.latent.ij[[i]]$ts.state[t.cur.i]
          pq.cur.i = c(p.cur.i,1 - p.cur.i)
          p.cur.j = ts.hyper.latent.ij[[j]]$ts.state[t.cur.j]
          pq.cur.j = c(p.cur.j,1 - p.cur.j)
          accept.ij = (pq.cur.i*pq.cur.j)
          
          if (runif(1) <= sum(accept.ij)) {
            marks[k] = sample(c(1,2),1,prob=accept.ij)
            ts.obs.ij = rbind(ts.obs.ij,c(times[k],marks[k]))
            ts.TKVs = rbind(ts.TKVs,c(times[k],marks[k],myIDmat[i,]+myIDmat[j,]))
          } else {
            ts.TKVs = rbind(ts.TKVs,c(times[k],0,myIDmat[i,]+myIDmat[j,]))
          }
          
          ts.accept.ij[ts.itr] = sum(accept.ij)
              ts.itr = ts.itr + 1
        }
    } 
    
  }
}

ts.TKVs.all = ts.TKVs[order(ts.TKVs[,1]),]
ts.TKVs.obs = ts.TKVs.all[ts.TKVs.all[,2]>0,]

outOBJ = ts.TKVs.obs[,c(1,(3:(2+nvertex)),2)]
write.table(outOBJ,filename,sep=" ",col.names=F,row.names=F);
