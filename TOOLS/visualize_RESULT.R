
myRESULTs = NULL
for(i in 0:15){
  out=read.table(paste("myRESULT_",i,".txt",sep=""))
  myRESULTs = rbind(myRESULTs,apply(out[25:50,],2,mean))
}

for(i in 1:16){
  pdf(file=paste("EXP#_",i,".pdf",sep=""));
  plot(1:12,myRESULTs[i,1:12],ylim=c(0,1),main=paste("EXP#_",i,sep=""),xlab="VERTEX ID",ylab="VERTEX LOCATION (EST)")
  dev.off()
}

