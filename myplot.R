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
