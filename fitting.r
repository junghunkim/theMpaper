maxEMstep = 50
maxMCrep = 50
mylist = vector("list",8)
mylist[[1]] = list(C1=c(1,2,3,4),C2=NULL)
mylist[[2]] = list(C1=c(1,2,3),C2 = c(4))
mylist[[3]] = list(C1=c(1,2,4),C2 = c(3))
mylist[[4]] = list(C1=c(1,3,4),C2 = c(2))
mylist[[5]] = list(C1=c(2,3,4),C2 = c(1))
mylist[[6]] = list(C1=c(1,2),C2=c(3,4))
mylist[[7]] = list(C1=c(1,3),C2=c(2,4))
mylist[[8]] = list(C1=c(1,4),C2=c(2,3))
mylist.est  = vector("list",8)
nstep = 1
for(list.itr in 1:8){
  C1.tmp = mylist[[list.itr]]$C1
  C2.tmp = mylist[[list.itr]]$C2
  toprecord = c(0.5,0.5) ###  first guess
  baserate.est = 10 #### first guess
  for(em.ITR in 1:maxEMstep) {
    cat("\n************ Starting the EM expriement #", em.ITR, "\n", sep="") 
    mystart.time = timestamp()
    myPsi.C1.tmp = rep(c(1,-2,toprecord[1]),length(C1.tmp))
    myPsi.C2.tmp = rep(c(1,-2,toprecord[2]),length(C2.tmp))
    myPsi.tmp = c(myPsi.C1.tmp,myPsi.C2.tmp)
#### E-step
    clock.mc = vector("list",nvertex)
    path.mc = vector("list",nvertex)
#   source("simulator.r")
    for(mcrep in 1:maxMCrep) {
      cond.sim.out = mysimulator.cond(ts.TKVs.obs = ts.TKVs.grid, cur.psi = matrix(myPsi.tmp,nrow=nvertex,byrow=T), maxTime=maxTime, mybaseRate=baserate.est,mcitr=mcrep)
      for(itr.v in 1:nvertex) {
        path.mc[[itr.v]] = cbind(path.mc[[itr.v]],cond.sim.out[[itr.v]]$state)
        clock.mc[[itr.v]] =cbind(clock.mc[[itr.v]],cond.sim.out[[itr.v]]$clock)
      }
    }
###########################
    temprecord = numeric(4)
    for(itr.v in 1:4) {
      temprecord[itr.v] = mean(apply(path.mc[[itr.v]],2,mean))
    }
#
    toprecord[1] = mean(temprecord[C1.tmp])
    toprecord[2] = mean(temprecord[C2.tmp])
    cat(toprecord,"\n")
#
    myend.time = timestamp()
    cat("************** Ending the EM experiment #", em.ITR, "\n", sep="") 
    em.ITR = em.ITR + 1
  }
  temp = matrix(0,nrow=nvertex,ncol=3)
  temp[,1] = 1
  temp[,2] = -2
  for(itr.v in C1.tmp) {
    temp[itr.v,3] = toprecord[1]
  }
  for(itr.v in C2.tmp) {
    temp[itr.v,3] = toprecord[2]
  }
  mylist.est[[list.itr]] =  temp
}

myloglik = vector("numeric",2)
for(vec.itr in 1:2) {
  myPsi.tmp = mylist.est[[vec.itr]]
  maxMCrep = 10
  for(mcrep in 1:maxMCrep) {
    cond.sim.out = mysimulator.cond(ts.TKVs.grid,
      myPsi.tmp, maxTime, baserate.est,mcrep)
    time.for.all = cond.sim.out[[1]]$clock
    temp.path = matrix(0,nrow=length(time.for.all),ncol=nvertex)
    for(itr.v in 1:nvertex) {      
      temp.path[,itr.v] = cond.sim.out[[itr.v]]$state
      myloglik[vec.itr] = myloglik[vec.itr] +
        (mydensity.vertex(mylist.est[[vec.itr]][itr.v,],temp.path[,itr.v],time.for.all))
    }
    for(itr.u in 1:(nvertex-1)) {
      for(itr.v in (itr.u+1):nvertex) {
        itr.u.X = cond.sim.out[[itr.u]]$state
        itr.u.Y = 1 - itr.u.X
        itr.v.Y = cond.sim.out[[itr.v]]$state
        itr.v.Y = 1 - itr.v.Y
        myloglik[vec.itr]=myloglik[vec.itr]-sum(time.for.all*(itr.u.X *itr.v.X + itr.u.Y*itr.v.Y))
      }
    }
  }
  index.itr = seq(1,240,by=10)
  for(obs.itr in index.itr) {
        cur.data = ts.TKVs.grid[obs.itr,]
        cur.sim = temp.path[obs.itr,]
        cur.uv = which(cur.data[3:6]==1) 
        u.tmp = cur.sim[cur.uv[1]]
        v.tmp = cur.sim[cur.uv[2]]
        u.tmp = c(u.tmp,1-u.tmp)
        v.tmp = c(v.tmp,1-v.tmp)
        myloglik[vec.itr] = myloglik[vec.itr]+log(u.tmp[cur.data[2]]*v.tmp[cur.data[2]])
  }
}



########################################################################

    M = matrix(0,nrow=nvertex,ncol=3) 
    for(id in 1:nvertex) {
      for(mc in 1:maxMCrep) {
        xt = path.mc[[id]][,mc]
	dT = c(0,diff(clock.mc[[id]][,mc]))
        xT = xt[length(xt)]
        x0 = xt[1]
        M[id,1] = M[id,1] + log((xT/(1-xT))/(x0/(1-x0)))
        M[id,2] = M[id,2] + sum(1/(xt-xt^2)*dT)
        M[id,3] = M[id,3] + sum((1/2-xt)/(xt-xt^2)*dT)
      }
    }
    
    myfcnC1 <- function(g) {
      b=-2
      s= 1
      if(length(C1.tmp) == 0) 
        return(0)
      COEF = vector("numeric",3)
      for( id in C1.tmp ) {
        for( coef in 1:3) {
          COEF[coef] = COEF[coef] + M[id,coef]
        }
      }
      retVal =0
      retVal = retVal + (-b/s^2) * COEF[1]
      retVal = retVal + -0.5*((b/s)^2*(1-2*g)*(-1)*COEF[2]+(b/s+s/2)*(b/s)*COEF[3])
      retVal = retVal + -0.5*b*COEF[3]
      return(retVal)
    }
    
    myfcnC2 <- function(g) {
      b = -2
      s = 1
      if(length(C2.tmp) == 0) 
        return(0)
      COEF = vector("numeric",3)
      for( id in C2.tmp ) {
        for( coef in 1:3) {
          COEF[coef] = COEF[coef] + M[id,coef]
        }
      }
      retVal =0
      retVal = retVal + (-b/s^2) * COEF[1]
      retVal = retVal + -0.5*((b/s)^2*(1-2*g)*(-1)*COEF[2]+(b/s+s/2)*(b/s)*COEF[3])
      retVal = retVal + -0.5*b*COEF[3]
      return(retVal)
    }
    
    toprecord[1] = uniroot(myfcnC1,c(0,1))$root
    toprecord[2] = uniroot(myfcnC2,c(0,1))$root
    timestamp()
