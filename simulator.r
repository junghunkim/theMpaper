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

mysimulator.cond <- function(ts.TKVs.obs=NULL,cur.psi=NULL,maxTime=NULL,mybaseRate=NULL,mcitr = NULL) {
  ts.TKVs.obs = rbind(ts.TKVs.obs,c(maxTime,0,rep(0,nvertex)),deparse.level=0)
    
    color.id.C1 = heat.colors(length(C1.tmp))
    color.id.C2 = grey.colors(length(C2.tmp))
    color.id = c(color.id.C1,color.id.C2)
    retObj = vector("list",nvertex)
    
    nobs = length(ts.TKVs.obs[,1])
    lhs = 0
    rhs = ts.TKVs.obs[1,1]

    plot(0,0,xlim=c(0,maxTime),ylim=c(0,1),main=paste("CUR-MC=",mcitr,sep=""))
#	plot(0,0,xlim=c(0,maxTime),ylim=c(0,1),main="",xlab="TIME", ylab="TOPIC 1");   
    cur.init = numeric(nvertex)
    for(itr.v in 1:nvertex) {
      cur.init[itr.v] = cur.psi[itr.v,3]
      }
    
    for(myitr in 1:nobs) {
      MORE = TRUE
                                        #            cat("----> simulating for the interval (",formatC(lhs,format="f",digit=7),",",formatC(rhs,format="f",digit=7),"]\n",sep="")
            nrej = 0
            while(MORE) {
              
              cur.path = NULL
                for(itr.v in 1:nvertex) {
                  cur.path[[itr.v]] = mysimulator(cur.psi[itr.v,],init=cur.init[itr.v],LHS=lhs,RHS=rhs)
                  }
                
                mycumsum = 0
                
                for(itr.i in 1:(nvertex-1)){
                  for(itr.j in (itr.i+1):(nvertex)){
                      
                      cur.path.i = cur.path[[itr.i]]
                    cur.path.j = cur.path[[itr.j]]
                    temp.path.i = cbind(cur.path.i$ts.state,1-cur.path.i$ts.state,deparse.level=0)
                    temp.path.j = cbind(cur.path.j$ts.state,1-cur.path.j$ts.state,deparse.level=0)
                    
                    temp.lambda.ij.first = as.vector(apply(cbind(temp.path.i[,1],temp.path.j[,1]),1,prod))
                    temp.lambda.ij.second = as.vector(apply(cbind(temp.path.i[,2],temp.path.j[,2]),1,prod))
                    temp.lambda.ij.both = cbind(temp.lambda.ij.first,temp.lambda.ij.second,deparse.level=0)                  
                    temp.lambda.ij = as.vector(apply(temp.lambda.ij.both,1,sum))
                    temp.end.i = c(cur.path.i$end,1-cur.path.i$end)
                    temp.end.j = c(cur.path.j$end,1-cur.path.j$end)
                    
                    temp.nrow = nrow(temp.lambda.ij.both)
                    temp.topic = ts.TKVs.obs[myitr,2]
                    
                    if(ts.TKVs.obs[myitr,2+itr.i] == 1 && ts.TKVs.obs[myitr,2+itr.j]==1) {
                      temp.lambda.ij.end = temp.lambda.ij.both[temp.nrow,temp.topic]
                      } else {
                      temp.lambda.ij.end = 1
                      }
                    
                    
                    mycumsum = mycumsum + log(temp.lambda.ij.end) + log(exp(-mybaseRate*mean(temp.lambda.ij)*(rhs-lhs)))
                  }
                  }
                
                if(log(runif(1)) < mycumsum) { # ACCEPT  = runif(1) < mycumsum
                                        #                cat("*****************************************************ACCEPT\n")
              MORE = FALSE
              } else {
                                        #                cat("***** REJECT #",nrej,"\n",sep="")
              nrej = nrej + 1
                MORE = TRUE
              }
              }
            
            if ( rhs < maxTime) {
        lhs = rhs
          rhs = ts.TKVs.obs[myitr+1,1]
        }
            
            for(itr.v in 1:nvertex) {
              cur.init[itr.v] = cur.path[[itr.v]]$end
              retObj[[itr.v]]$clock = c(retObj[[itr.v]]$clock,(cur.path[[itr.v]]$ts.time)[-1])
              retObj[[itr.v]]$state = c(retObj[[itr.v]]$state,(cur.path[[itr.v]]$ts.state)[-1])
              lines(retObj[[itr.v]]$clock,retObj[[itr.v]]$state,col=color.id[itr.v])
            }
          }
    

#  cat("****************************FINISHED\n")

  retObj
}

mydensity.vertex <-function(psi,grid.time,path) {
    first.index = 1
    last.index = length(path)    

    path = (1/psi[1])*asin(2*path-1)

    A = (psi[2]/psi[1]^2+1/2)*log(1/(cos(psi[1]*path)))
    B = psi[2]/psi[1]^2*(1-2*psi[3])*log(1/cos(psi[1]*path)+tan(psi[1]*path))
    C1x = A + B

    A= (psi[2]/psi[1]+psi[1]/2) * tan(psi[1]*path)
    B= psi[2]/psi[1]*(1-2*psi[3])* 1/cos(psi[1]*path)
    C2x = (A + B)^2
    
    A = (psi[2]+psi[1]^2/2)/(cos(psi[1]*path))^2
    B = psi[2]/psi[1]*(1-2*psi[3])*(1/(cos(psi[1]*path)))*tan(psi[1]*path)
    C3x = A + B
#################

    retValx = C1x[last.index] - C1x[first.index]
    
    dretValx = (-0.5)*sum(C2x*c(0, diff(grid.time)))
    retValx = retValx + dretValx
    dretValx = (-0.5)*sum(C3x*c(0, diff(grid.time)))
    retValx = retValx + dretValx

    as.numeric(retValx)
}
