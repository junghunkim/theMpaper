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

master.clock = cond.sim.out[[1]]$clock

myobjfcn<- function(psi) {
    psi = matrix(psi,nrow=nvertex,byrow=T)
    retVal = 0

    for(mcitr in 1:maxMCrep) {
          for(itr.v in 1:nvertex) {
                retVal = retVal + mydensity.vertex(psi[itr.v,],clock.mc[[itr.v]][,mcitr],path.mc[[itr.v]][,mcitr])
          }
    }
    retVal/maxMCrep
}
