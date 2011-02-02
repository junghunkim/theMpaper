source("simulator.r")

for(mcrep in 1:maxMCrep) {
    cat("\n************** Starting the ", mcrep, "-th experiment *********** \n", sep="") 
    cat(timestamp(),"\n")
    cond.sim.out = mysimulator.cond(ts.TKVs.obs = ts.TKVs.grid, cur.psi = matrix(myPsi.est,nrow=nvertex,byrow=T), maxTime=maxTime, mybaseRate=baserate.est)
   for(itr.v in 1:nvertex) {
        path.mc[[itr.v]] = cbind(path.mc[[itr.v]],cond.sim.out[[itr.v]]$state)
        clock.mc[[itr.v]] =cbind(clock.mc[[itr.v]],cond.sim.out[[itr.v]]$clock)
    }
    cat("\n")
}

M = matrix(0,nrow=nvertex,ncol=7) 
for(id in 1:nvertex) {
  for(mc in 1:maxMCrep) {
    xt = path.mc[[id]][mc,]
	dT = clock.mc[[id]][mc,]

    xT = xt[length(xt)]
    x0 = xt[1]

    M[id,1] = M[id,1] + (1/2)*log((xT-xT^2)/(x0-x0^2))
    M[id,2] = M[id,2] + log((xT/x0)/((1-xT)/(1-x0)))
    M[id,3] = M[id,3] + sum(1/(1-xt)*dT)
    M[id,4] = M[id,4] + sum(1/xt*dT)
    M[id,5] = M[id,5] + sum(1/sqrt(xt-xt^2)*dT)
    M[id,6] = M[id,6] + sum((1/sqrt(xt-xt^2))*(1/(1-xt))*(1/2-xt)*dT)
    M[id,7] = M[id,7] + sum((1/sqrt(xt-xt^2))^3*(1/2-xt)*dT)
  }
}

myfcnC1 <- function(b,s,g,C1=1:n) {
	  if(length(C1) == 0) 
		return(0);

      COEF = vector("numeric",7)

      for( id in C1 ) {
        for( coef in 1:7) {
          COEF[coef] = COEF[coef] + M[id,coef]
        }
      }
      
      retVal = -(b/s^2 + 1/2) * COEF[1] + 
      b*(1-2*g)/(2*s^2) * COEF[2]
+
      (-0.5)*( (-b/s + s/2)^2*(tRHS-tLHS)*nvetex.C1 +
      ((b/s+s/2) - (b*g/s + s/4))^2 * COEF[3] +
      (b*g/s + s/4)^2 *COEF[4] )
+
      (-0.5)* ( (b/s+s/2) *(COEF[5]- COEF[6]) + 
      (b*g/s+s/4)* COEF[7] )
	  
	  return(retVal)
}

myfcnC2 <- function(b,s,g,C2=NULL) {
	  if(length(C1) == 0) 
		return(0);

      COEF = vector("numeric",7)

      for( id in C2 ) {
        for( coef in 1:7) {
          COEF[coef] = COEF[coef] + M[id,coef]
        }
      }
      
      retVal = -(b/s^2 + 1/2) * COEF[1] + 
      b*(1-2*g)/(2*s^2) * COEF[2]
+
      (-0.5)*( (-b/s + s/2)^2*(tRHS-tLHS)*nvetex.C1 +
      ((b/s+s/2) - (b*g/s + s/4))^2 * COEF[3] +
      (b*g/s + s/4)^2 *COEF[4] )
+
      (-0.5)* ( (b/s+s/2) *(COEF[5]- COEF[6]) + 
      (b*g/s+s/4)* COEF[7] )

	  return(retVal)
}

myobjfcn <-function(myPsi) {
      A = myfcnC1(myPsi[1],myPsi[2],myPsi[3])
      B = myfcnC2(myPsi[4],myPsi[5],myPsi[6])
      A+B
}

output = optim(c(,,,,,),myfcn,method="SANN")
