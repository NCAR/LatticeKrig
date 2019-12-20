setwd("~/Dropbox/Home/Projects/BRACEData")
load("slopeBoulder.rda")
#image.plot( Longitude, Latitude,TS[,,1] )
#points( Longitude[jBoulder], Latitude[ kBoulder], pch="+" )
library( fields)

fields.style()
# matplot(year, cbind( tempX[,1:2], tempY[,1:2]), 
#        type=c("l","l","p","p"), pch=16,
#         xlab="Year",
#        col=c("green2", "grey40", "green2", "grey40"),
#         lty=1,
#        lwd=3, ylab="Degrees C" )

pdf("pix/BoulderEx1A.pdf")
fields.style()
matplot(year, cbind( tempX[,1], tempY[,1] ), 
            type=c("l","p"), pch=16,
              xlab="Year",
             col=c("grey40", "green3"),
            lty=1,
             lwd=3, ylab="Degrees C" )
dev.off()

pdf("pix/BoulderEx1B.pdf", width=6, height=6)
fields.style()
plot( tempX[,1], tempY[,1] , 
        pch=16,type="p",
        xlab="Global Average",
        ylab="Boulder GridBox",
       col="orange3"
        )
obj<- lm(tempY[,1] ~ tempX[,1] )
abline( obj, col="grey", lwd=3)
dev.off()

pdf("pix/BoulderEx2.pdf")
fields.style()
matplot(year, cbind( tempX[,2], tempY[,2] ), 
        type=c("l","p"), pch=16,
        xlab="Year",
        col=c("grey40", "green3"),
        lty=1,
        lwd=3, ylab="Degrees C" )
dev.off()

pdf("pix/BoulderEx2B.pdf", width=6, height=6)
fields.style()
matplot( tempX[,1:2], tempY[,1:2] , 
      pch=16,type="p",
      xlab="Global Average",
      ylab="Boulder GridBox",
      col=c("orange3", "green4")
)
obj<- lm(tempY[,1] ~ tempX[,1] )
abline( obj, col="orange1", lwd=3)
obj<- lm(tempY[,2] ~ tempX[,2] )
abline( obj, col="green1", lwd=3)

dev.off()

load("slopeBoulder.rda")
b1<-b2<- rep( NA, 30)
for (k in 1:30){
  obj<- lm(tempY[,k] ~ tempX[,k] )
  b2[k]<- obj$coef[2]
  b1[k]<- obj$coef[1]
}
for( k in 1:30){
  abline( b1[k],b2[k], col="grey")
}




