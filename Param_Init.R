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

filename ="myoutputfile.txt"
