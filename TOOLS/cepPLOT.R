w38 = read.table("http://cis.jhu.edu/~nhlee/theDexperiment/w38_RESULT.txt",header=F)
w58 = read.table("http://cis.jhu.edu/~nhlee/theDexperiment/w58_RESULT.txt",header=F)

y38 = apply(w38[-c(1:25),],2,mean)
y58 = apply(w58[-c(1:25),],2,mean)

X = cbind(rep(2,34),rep(3,34)) 
Y = cbind(y38[1:34],y58[1:34])
XY = cbind(X,Y)

par(mfrow=c(1,2))
plot(1,1,type="n",xlim=c(2,3),ylim=c(0,1),ylab="Xhat",xlab="Time",main="theDmodel(Enron) - GROUP 1");
apply(XY[1:10,],1,function(x) {points(x[1:2],x[3:4],col="red");lines(x[1:2],x[3:4],col="red",lty="dotted")})

plot(1,1,type="n",xlim=c(2,3),ylim=c(0,1),ylab="Xhat",xlab="Time",main="theDmodel(Enron) - GROUP 2");
apply(XY[11:34,],1,function(x) {points(x[1:2],x[3:4],col="black");lines(x[1:2],x[3:4],col="black",lty="dotted")})

GROUP_1 <- XY[1:10,3:4]
GROUP_2 <- XY[11:34,3:4]
Lambda <- c(y38[35],y58[35])

crudeEst = rbind(GROUP_1,GROUP_2,Lambda)
colnames(crudeEst) =c("PERIOD I", "PERIOD II")
write.table(crudeEst,"CrudeEst.txt",row.names=F,col.names=F)

STAT.mean = rbind(c(mean(GROUP_1[,1]),mean(GROUP_1[,2])),
c(mean(GROUP_2[,1]),mean(GROUP_2[,2])))
rownames(STAT.mean) = c("GROUP 1", "GROUP 2")
colnames(STAT.mean) = c("PERIOD I", "PERIOD II")

STAT.sd = rbind(c(sd(GROUP_1[,1]),sd(GROUP_1[,2])),
c(sd(GROUP_2[,1]),sd(GROUP_2[,2])))
rownames(STAT.sd) = c("GROUP 1", "GROUP 2")
colnames(STAT.sd) = c("PERIOD I", "PERIOD II")

STAT.mean
STAT.sd
