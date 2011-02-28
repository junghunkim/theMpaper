mylist = vector("list",100)

for(i in 1:100) {
	tmp = paste("{",i-1,"}_param_estimate.mat",sep="")
	mylist[[i]] = read.table(tmp,header=F)	
}

tmp = paste("{",1,"}_data_tuvk.mat",sep="")
mytempdata = as.matrix(read.table(tmp,header=F))

dev.new()
topic.id = c(8,21)
plot(mytempdata[,1],rep(0,nrow(mytempdata)),type="n",xlim=c(0,1),ylim=c(0.5,4.5),xlab="TIME",ylab="ACTOR",main="Example of T UV K pattern")
axis(2,at=1:4)
for(itr in 1:nrow(mytempdata)) {
      cur.data = mytempdata[itr,]
      cur.time = cur.data[1]
      cur.topic = topic.id[cur.data[6]]
      cur.pair = which(cur.data[2:5] == 1)
      points(cur.time,cur.pair[1],pch=cur.topic,cex=1)
      points(cur.time,cur.pair[2],pch=cur.topic,cex=1)
}
legend(0.8,4.6,  c("TOPIC 1","TOPIC 2"), pch = c(8,21), text.col = "green4", bg = 'gray90')


dev.new() 
par(mfrow=c(2,2))
for(i in 1:4) {
	plot.ts(TRUEparam[,i],ylim=c(0,1.5),main=paste("VERTEX",i),ylab="Location Parameter")
	lines(ESTMparam[,i],lty="dotted",col="red")
	legend(0, 1.5, c("TRUE","ESTM"), col = c("black","red"), text.col = "green4", lty = c("solid", "dotted"), merge = TRUE, bg = 'gray90')
}

dev.new()
par(mfrow=c(2,2));
for(i in 1:4) { 
	plot(density(TRUEparam[,i]),xlim=c(0,1),ylim=c(0,3),xlab="Location Parameter",main=paste("VERTEX",i));
	abline(v=mean(TRUEparam[,i]),col="black")
	lines(density(ESTMparam[,i]),col="red",lty="dotted");
	abline(v=mean(ESTMparam[,i]),col="red",lty="dotted");
	if( i <=2) {
		x.loc = 0.625;y.loc=3.0;
	} else {
		x.loc = 0;y.loc=3.0;
	}
	legend(x.loc, y.loc, c("TRUE","ESTM"), col = c("black","red"), text.col = "green4", lty = c("solid", "dotted"), merge = TRUE, bg = 'gray90')
}

     
 
