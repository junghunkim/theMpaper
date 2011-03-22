myK_Labler = function(x) { 
  ifelse(sum(x == c(0,5,17,20,30)),2,1)
}
  
myUV_Expander = function(x) {
  retVec = numeric(184)
  retVec[x[1]] = 1
  retVec[x[2]] = 1
  retVec
}

set.seed(123)
myRED = sort(sample(c(82, 94, 129, 6, 39, 21, 37, 16, 93, 181),10))
myBLA = sort(sample((1:184)[-myRED],24))
myALL = (c(myRED,myBLA))
  
myURL = "http://cis.jhu.edu/~parky/Enron/execs.email.linesnum.ldctopic"
myDATA.ALL = read.table(myURL)
myDATA.ALL$time = as.POSIXct(myDATA.ALL$time,origin=ISOdatetime(1970,1,1,0,0,0))

w18.start.time <- as.POSIXct("2000-12-18 01:00:00 GMT", tz = "GMT")
w38.start.time <- as.POSIXct("2001-05-07 01:00:00 GMT", tz = "GMT")
w58.start.time <- as.POSIXct("2001-09-24 01:00:00 GMT", tz = "GMT")
w78.start.time <- as.POSIXct("2002-02-11 01:00:00 GMT", tz = "GMT")

nmsg = 125409
myDdata = vector("list",4)
for(itr in 1:nmsg) {
  cat("MSG # == ",itr,"/",nmsg,"\n",sep="")
  myT = myDATA.ALL$time[itr]
  myV = cbind(myDATA.ALL$sender[itr], myDATA.ALL$receiver[itr],deparse.level=0)
  myK = myDATA.ALL$LDCtopic[itr]
  myMSG = cbind(myT,myV,myK,deparse.level=0)
  if(w18.start.time <= myT && myT < w38.start.time) {
    tmp = 1
  }
  else if(w38.start.time <= myT && myT < w58.start.time) {
    tmp = 2
  }
  else if(w58.start.time <= myT && myT < w78.start.time){
    tmp = 3
  }
  else {
    tmp = 4
  }
  myDdata[[tmp]] = rbind(myDdata[[tmp]],myMSG)
}  

w38.myDdata = myDdata[[2]]
w58.myDdata = myDdata[[3]]

write.table(w38.myDdata,"w38myDdata.txt",sep=" ", row.names=F,col.names=F)
write.table(w58.myDdata,"w58myDdata.txt",sep=" ", row.names=F,col.names=F)



#############################################################################
w38.myDdata = read.table("w38myDdata.txt", header=F)
nmessages = nrow(w38.myDdata)
myT = as.POSIXct(w38.myDdata[,1],origin=ISOdatetime(1970,1,1,0,0,0),tz="GMT")
myV = cbind(w38.myDdata[,2],w38.myDdata[,3],deparse.level=0)
myK = w38.myDdata[,4]
base.time <- as.POSIXct(ISOdatetime(2001,5,7,1,0,0),tz = "GMT")

w38newDdata = matrix(0,nrow=nmessages,ncol=1+184+1)
for(itr in 1:nmessages){
  cat("MSG # == ",itr,"/",nmessages,"\n",sep="")
  w38newDdata[itr,1] = as.numeric(myT[itr] - base.time,unit="hours")/(150*24)  
  w38newDdata[itr,2:185] = myUV_Expander(myV[itr,])
  w38newDdata[itr,186] = myK_Labler(myK[itr])
}

nmessages = nrow(w38newDdata)
curmsg = w38newDdata[1,]
w38newMSG = NULL
for(itr in 2:nmessages) {
  cat("MSG # == ",itr,"/",nmessages,sep="")
  pasmsg = curmsg
  curmsg = w38newDdata[itr,]
  if(!identical(curmsg[1],pasmsg[1])){
    if(itr < nmessage){
      nxtmsg = w38newDdata[itr+1,]
      if(!identical(curmsg[1],nxtmsg[1])){
        w38newMSG =c(w38newMSG,itr)
        cat("*\n");
      } else {
        cat("\n")
      }
    } else {
      w38newMSG =c(w38newMSG,itr)
      cat("*\n");
    }
  } else {
    cat("\n")
  }
}

w38myDdata = w38newDdata[w38newMSG,c(1,myALL+2,186)]

nmessages = nrow(w38myDdata)
w38pairMSG = NULL
for(itr in 1:nmessages){
  cat("MSG # == ",itr,"/",nmessages,sep="")
  if(sum(w38myDdata[itr,-c(1,36)]) == 2) {
    w38pairMSG = c(w38pairMSG,itr)
    cat("*\n")
  } else {
    cat("\n")
  }  
}

w38myDdata = w38myDdata[w38pairMSG,]
str(w38myDdata)

dev.new()
myNewNewData = w38myDdata
plot(c(0,0),xlim=c(0,0.8),ylim=c(0,34),type="n",ylab="VERTEX ID")
nmessages = nrow(myNewNewData)  
for(itr in 1:nmessages) {
  tempT = myNewNewData[itr,1]
  tempUV = myNewNewData[itr,-c(1,36)];
  tempUV = which(tempUV>0)
  tempK = myNewNewData[itr,36]
  points(rbind(c(tempT,tempUV[1])),col=tempK)
  points(rbind(c(tempT,tempUV[2])),col=tempK)
  cat(tempT,"--",tempUV,"\n")
}
fname = paste("myDdata_w38.txt",sep="")
write.table(myNewNewData,fname,sep=" ",col.names=F,row.names=F);

########################################################################
w58.myDdata = read.table("w58myDdata.txt", header=F)
nmessages = nrow(w58.myDdata)
myT = as.POSIXct(w58.myDdata[,1],origin=ISOdatetime(1970,1,1,0,0,0),tz="GMT")
myV = cbind(w58.myDdata[,2],w58.myDdata[,3],deparse.level=0)
myK = w58.myDdata[,4]
base.time <- as.POSIXct(ISOdatetime(2001,9,24,1,0,0),tz = "GMT")

w58newDdata = matrix(0,nrow=nmessages,ncol=1+184+1)
for(itr in 1:nmessages){
  cat("MSG # == ",itr,"/",nmessages,"\n",sep="")
  w58newDdata[itr,1] = as.numeric(myT[itr] - base.time,unit="hours")/(150*24)  
  w58newDdata[itr,2:185] = myUV_Expander(myV[itr,])
  w58newDdata[itr,186] = myK_Labler(myK[itr])
}

nmessages = nrow(w58newDdata)
curmsg = w58newDdata[1,]
w58newMSG = NULL
for(itr in 2:nmessages) {
  cat("MSG # == ",itr,"/",nmessages,sep="")
  pasmsg = curmsg
  curmsg = w58newDdata[itr,]
  if(!identical(curmsg[1],pasmsg[1])){
    if(itr < nmessage){
      nxtmsg = w58newDdata[itr+1,]
      if(!identical(curmsg[1],nxtmsg[1])){
        w58newMSG =c(w58newMSG,itr)
        cat("*\n");
      } else {
        cat("\n")
      }
    } else {
      w58newMSG =c(w58newMSG,itr)
      cat("*\n");
    }
  } else {
    cat("\n")
  }
}

w58myDdata = w58newDdata[w58newMSG,c(1,myALL+2,186)]

nmessages = nrow(w58myDdata)
w58pairMSG = NULL
for(itr in 1:nmessages){
  cat("MSG # == ",itr,"/",nmessages,sep="")
  if(sum(w58myDdata[itr,-c(1,36)]) == 2) {
    w58pairMSG = c(w58pairMSG,itr)
    cat("*\n")
  } else {
    cat("\n")
  }  
}

w58myDdata = w58myDdata[w58pairMSG,]
w58myDdata = w58myDdata[which(w58myDdata[,1]<0.8),]
str(w58myDdata)

dev.new()
myNewNewData = w58myDdata
plot(c(0,0),xlim=c(0,0.8),ylim=c(0,34),type="n",ylab="VERTEX ID")
nmessages = nrow(myNewNewData)  
for(itr in 1:nmessages) {
  tempT = myNewNewData[itr,1]
  tempUV = myNewNewData[itr,-c(1,36)];
  tempUV = which(tempUV>0)
  tempK = myNewNewData[itr,36]
  points(rbind(c(tempT,tempUV[1])),col=tempK)
  points(rbind(c(tempT,tempUV[2])),col=tempK)
  cat(tempT,"--",tempUV,"\n")
}
fname = paste("myDdata_w58.txt",sep="")
write.table(myNewNewData,fname,sep=" ",col.names=F,row.names=F);


myResult = matrix(0,nrow=3,ncol=13)
myResult[1,] = apply(read.table("DONE/myRESULT_0.txt",header=F)[25:50,],2,mean)
myResult[2,] = apply(read.table("DONE/myRESULT_1.txt",header=F)[25:50,],2,mean)
myResult[3,] = apply(read.table("DONE/myRESULT_2.txt",header=F)[25:50,],2,mean)
myResult[1,13] = myResult[1,13]/4

row.names(myResult) <- c("P1","P2","P3")
colnames(myResult) <- c(sapply(seq(1:12),function(x) paste("V",x,sep="")),"LAMBDA")
myResult = format(myResult,digit=5)
write.table(myResult,"myRESULT_summary.txt",sep=" ", row.names=F,col.names=F)




