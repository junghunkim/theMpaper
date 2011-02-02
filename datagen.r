rm(list=ls())
source("simulator.r")

maxTime = 1
nstep = 100
nvertex = 4

myunitRate = 10
myPsi.C1 = rep(c(1,-2,0.5),2)
myPsi.C2 = rep(c(1,-2,0.5),2)
myPsi = c(myPsi.C1,myPsi.C2)
myPsi = matrix(myPsi,byrow=T,nrow=nvertex)

myinit = as.vector(myPsi[,3])
myIDmat = diag(nvertex)


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
  cat("************ Current (i,j) iteration : (", i,",",j,")\n")   
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

time.ts = ts.TKVs.obs[,1]
nobs = length(time.ts)
topic.id = c("red","black")
plot(time.ts,rep(0,nobs),type="n",xlim=c(0,1),ylim=c(0.5,nvertex+0.5),xlab="TIME",ylab="ACTOR")
for(itr in 1:nobs) {
	cur.data = ts.TKVs.obs[itr,]
	cur.time = cur.data[1]
	cur.topic = topic.id[cur.data[2]]
        cur.pair = which(cur.data[3:(2+nvertex)] == 1) 
	points(cur.time,cur.pair[1],col=cur.topic,cex=0.5)
	points(cur.time,cur.pair[2],col=cur.topic,cex=0.5)
}
axis(2,at=1:nvertex)


nobs = nrow(ts.TKVs.obs)
ngrids = 10

ts.TKVs.grid = matrix(0,nrow=nobs*ngrids,ncol=(2+nvertex))
itr.row = 1

for(itr.nobs in 1:nobs){
  cat("@ ITR=", itr.row," out of ",ngrids*nobs,"\n",sep="")
  ts.TKVs.grid[itr.row,] = ts.TKVs.obs[itr.nobs,]
  itr.row= itr.row+1
  
  LHS = ts.TKVs.obs[itr.nobs,1]
  if(itr.nobs < nobs) {
    RHS = ts.TKVs.obs[itr.nobs+1,1]
  } else {
    RHS = maxTime
  }
  
  dT = (RHS-LHS)/ngrids

  for(itr.grid in 1:(ngrids-1)) {
    cat("@ ITR=", itr.row," out of ",ngrids*nobs,"\n",sep="")
    lhs = LHS+(itr.grid)*dT
    ts.TKVs.grid[itr.row,] = c(lhs,rep(0,nvertex+1))
    itr.row = itr.row + 1
  }
  
}

save(file="mydata.RData")

