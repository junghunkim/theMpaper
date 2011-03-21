myURL = "http://www.cis.jhu.edu/~parky/Enron/coe_enron_execs_one_to_per_line.RData"
myData = dget(myURL)

myK_Labler = function(x) { 
  ifelse(sum(x == c(0,5,17,20,30)),2,1)
}

myV_Labler = function(x) {
  retVec = numeric(184)
  retVec[x[1]] = 1
  retVec[x[2]] = 1
  retVec
}
	
myT = myData$time
myV = cbind(myData$from,myData$to,deparse.level=0)
myK = myData$topic

nmessages = length(myT)

myNewData = matrix(0,nrow=nmessages,ncol=1+184+1)
myNewData[,1] = myT
myNewData[,2:185] = apply(myV,1,myV_Labler)
myNewData[,186] = sapply(myK,myK_Labler)

set.seed(123)
myRED = sort(sample(c(82, 94, 129, 6, 39, 21, 37, 16, 93, 181),5))
myBLA = sort(sample((1:184)[-myRED],7))
myALL = (c(myRED,myBLA))

myNewData[,1] = (myNewData[,1] - min(myNewData[,1]))/3600/24
myINDEX = apply(myNewData[,myALL+2],1,function(x) { ifelse(sum(x)==2,T,F) })
myNewData = myNewData[myINDEX,c(1,myALL+2,186)]

for(ITV.itr in 0:8) {
  myLHS = 7000 + ITV.itr * 120
  myRHS = myLHS + 120
  myINDEX = (myNewData[,1] > myLHS) & (myNewData[,1] < myRHS)
  myNewNewData = myNewData[myINDEX,]
  myNewNewData[,1] = jitter(myNewNewData[,1],amount=0.01)
  myINDEX = order(myNewNewData[,1])
  myNewNewData = myNewNewData[myINDEX,]
  myNewNewData[,1] = (myNewNewData[,1] - myLHS)/120
  plot(c(0,0),xlim=c(0,1),ylim=c(0,14),type="n",ylab="VERTEX ID")
  nmessages = nrow(myNewNewData)  
  for(itr in 1:nmessages) {
    tempT = myNewNewData[itr,1]
    tempUV = myNewNewData[itr,-c(1,14)];tempUV = which(tempUV>0)
    tempK = myNewNewData[itr,14]
    points(rbind(c(tempT,tempUV[1])),col=tempK)
    points(rbind(c(tempT,tempUV[2])),col=tempK)
    cat(tempT,"--",tempUV,"\n")
  }
  fname = paste("DATA/myData_",ITV.itr,".txt",sep="")
  write.table(myNewNewData,fname,sep=" ",col.names=F,row.names=F);
}


myResult = matrix(0,nrow=3,ncol=13)
myResult[1,] = apply(read.table("DONE/myRESULT_0.txt",header=F)[25:50,],2,mean)
myResult[2,] = apply(read.table("DONE/myRESULT_1.txt",header=F)[25:50,],2,mean)
myResult[3,] = apply(read.table("DONE/myRESULT_2.txt",header=F)[25:50,],2,mean)
myResult[1,13] = myResult[1,13]/4

row.names(myResult) <- c("P1","P2","P3")
colnames(myResult) <- c(sapply(seq(1:12),function(x) paste("V",x,sep="")),"LAMBDA")
myResult = format(myResult,digit=5)
write.table(myResult,"myRESULT_summary.txt",sep=" ", row.names=F,col.names=F)




