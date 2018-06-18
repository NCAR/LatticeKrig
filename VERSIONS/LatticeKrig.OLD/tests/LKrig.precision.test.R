
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1

# test building MRF
  mx<- 5
my<- 6
m<- mx*my
temp<- matrix( 1:m, mx,my)
   I.c<- temp
   I.B<-  cbind(rep(NA,mx),temp[,1:(my-1)])
   I.T<-  cbind( temp[,2:my], rep(NA,mx))
   I.L<-  rbind( rep(NA, my),temp[ 1:(mx-1),] )
   I.R<-  rbind( temp[2:mx, ], rep( NA, my))
#   t( I.c)[my:1,];  t( I.T)[my:1,]; t( I.B)[my:1,];  t( I.L)[my:1,]; t( I.R)[my:1,]  
   Bi<- rep( 1:30, 5)
   Bj<- cbind( c(I.c), c(I.T), c(I.B), c( I.L), c(I.R))
   atest<- matrix( (1:m)*5, mx,my)
   values<- cbind( c(atest), rep(-1,m),rep(-1,m),rep(-1,m),rep(-1,m) )
   good<- !is.na(c(Bj))
   Bi<- Bi[good]
   Bj<- Bj[good]
   values<- c(values)[c(good)]
   obj<- LKrig.MRF.precision( mx,my,a.wght=atest,stationary=FALSE, edge=FALSE)

   test.for.zero( cbind( Bi,Bj), obj$ind, tag="MRF index")
   test.for.zero( values, obj$ra, tag="MRF value")

   atest<- matrix( 1:(30*5), 5,6)
   obj<- LKrig.MRF.precision( 5,6,a.wght=atest,stationary=FALSE,  edge=FALSE)
   values<- cbind( c(atest), rep(-1,30), rep( -1,30), rep(-1,30), rep( -1,30))
   values<- values[good]
   test.for.zero( values, obj$ra, tag="MRF a.wght as matrix")

# filling out precision matrix to second nearest neighbors.
# 
   mx<- 4
   my<- 5
   m<- mx*my
   temp<- matrix( 1:m, mx,my)
   obj<- LKrig.MRF.precision( mx,my,a.wght= 1:9,stationary=TRUE,  edge=FALSE)
   pmat1<- spind2full(obj)
# hard coded indices up to 2nd order neighbors for 4X5
 ind0<- matrix(
 c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20, 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 6, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20, 1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 2, 3, 4, 6, 7, 8, 10, 11, 12, 14, 15, 16, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 6, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20), ncol=2)
# visual check
set.panel(5,4)
pmat<- matrix(0,20,20)
pmat[ ind0] <-1
par( mar=c(0,1,1,0))
k<- 0
for ( k in 1:20){
  temp<- matrix( pmat[k,],4,5)
  image(1:4, 1:5, temp, axes=FALSE, xlab="", ylab="", col=c("white",rainbow(9)))
  xline((0:4)+.5)
  yline( (0:5)+.5)
}


   obj<- LKrig.MRF.precision( 5,6,a.wght=matrix((1:30)*5,5,6),stationary=FALSE,
                                   edge=FALSE)
   obj2<-spind2full( obj)
   test2<- matrix( obj2[8,], 5,6)
   test0<- matrix( 0, 5,6)
   test0[3,2]<- 5*8
   test0[2,2]<- test0[3,1]<- test0[4,2]<- test0[3,3] <- -1
   test.for.zero( test2, test0, tag="MRF weight placement")
   



## cylinder tests
mx<- 5
my<- 6
m<- mx*my
temp<- matrix( 1:m, mx,my)
   I.c<- temp
   I.B<-  cbind(rep(NA,mx),temp[,1:(my-1)])
   I.T<-  cbind( temp[,2:my], rep(NA,mx))
   I.L<-  temp[c(mx,1:(mx-1)), ]
   I.R<-  temp[c(2:mx,1), ]
#   t( I.c)[my:1,];  t( I.T)[my:1,]; t( I.B)[my:1,];  t( I.L)[my:1,]; t( I.R)[my:1,]  
   Bi<- rep( 1:30, 5)
   Bj<- cbind( c(I.c), c(I.T), c(I.B), c( I.L), c(I.R))
   atest<- matrix( (1:m)*5, mx,my)
   values<- cbind( c(atest), rep(-1,m),rep(-1,m),rep(-1,m),rep(-1,m) )
   good<- !is.na(Bj)
   Bi<- Bi[good]
   Bj<- Bj[good]
   values<- c(values)[c(good)]
   obj<- LKrig.MRF.precision( mx,my,a.wght=atest, stationary=FALSE, distance.type="cylinder", edge=FALSE)

   test.for.zero( cbind( Bi,Bj), obj$ind)
   test.for.zero( values, obj$ra)

   obj<- LKrig.MRF.precision( 5,6,a.wght=5,stationary=TRUE,distance.type="cylinder", edge=FALSE)
   obj2<-spind2full( obj)
   test2<- matrix( obj2[8,], 5,6)
   test0<- matrix( 0, 5,6)
   test0[3,2]<- 5
   test0[2,2]<- test0[3,1]<- test0[4,2]<- test0[3,3] <- -1
   test.for.zero( test2, test0, tag="MRF cyl weight placement")

   test2<-matrix( obj2[1,], 5,6)
   test0<- matrix( 0, 5,6)
   test0[1,1]<- 5
   test0[5,1]<- test0[1,2]<- test0[2,1] <- -1
   test.for.zero( test2, test0, tag="MRF cyl corner weight placement")

   test2<-matrix( obj2[30,], 5,6)
   test0<- matrix( 0, 5,6)
   test0[5,6]<- 5
   test0[4,6]<- test0[1,6]<- test0[5,5] <- -1
   test.for.zero( test2, test0, tag="MRF cyl corner weight placement")

############ end of simple MRF tests
### a visual test for cylinder
#   mx<- 7;  my<- 9; m<- mx*my
#   obj<- spind2full(LKrig.MRF.precision( mx,my,a.wght=5, distance.type="cylinder", edge=TRUE))
#   par( mar=c(0,0,0,0), mfcol=c(my,mx))
#  for( k in c(t(matrix(1:m,mx,my))[my:1,])){ image( matrix(obj[,k],mx,my), axes=FALSE, zlim=c(-1,5), col=terrain.colors(256)); box(col="grey", lwd=2);title(paste(k), line=-1)}


LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                     a.wght=c(5,6,7), alpha=c( 1,1,1),
                     edge=FALSE, normalize=FALSE)

hold<- LKrig.precision( LKinfo, return.B=TRUE)
hold<- spam2full(hold)

test.for.zero( diag(hold), rep( c(5,6,7), LKinfo$mx*LKinfo$my),
                     tag="diagonal elements of precision 3-levels")
hold2<- LKrig.precision( LKinfo, return.B=TRUE, level.index=2)
hold2<- spam2full( hold2)
number.level<-  LKinfo$mx[2]*LKinfo$my[2]
ind2<-ind1<-  (1:number.level) + LKinfo$offset[2]
test.for.zero( c( hold[ind1, ind2]), c(hold2), tag="just level 2 B matrix")

# now test  t(B)%*%B

hold<- LKrig.precision( LKinfo)
hold<- spam2full(hold)
hold2<- LKrig.precision( LKinfo, level.index=2)
hold2<- spam2full( hold2)
number.level<-  LKinfo$mx[2]*LKinfo$my[2]
ind2<-ind1<-  (1:number.level) + LKinfo$offset[2]
test.for.zero( c( hold[ind1, ind2]), c(hold2), tag="just level 2 Q matrix")

# now everything
  LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4, a.wght=c(5,6,7)
                       , alpha=c(1,1,1))

  hold<- LKrig.precision( LKinfo)
  hold<- spam2full(hold)
  hold2<- LKrig.precision( LKinfo, level.index=2)
  hold2<- spam2full( hold2)
  number.level<-  LKinfo$mx[2]*LKinfo$my[2]
  ind2<-ind1<-  (1:number.level) + LKinfo$offset[2]
  test.for.zero( c( hold[ind1, ind2]), c(hold2), tag= "level 2 Q normalization and edge")

  set.seed(123)
  x1<- cbind( runif( 10), runif(10))
  LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                                                 a.wght=c(5,6,7), alpha=c(6,6,6))
  X<- LKrig.basis(x1, LKinfo)
  X<- spam2full(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- spam2full( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  marginal.var<- sum(unlist(LKinfo$alpha))
  test.for.zero( diag(look), rep( marginal.var,10),
                tag="normalization to unit variance at each level")
  look2<- LKrig.cov( x1, LKinfo= LKinfo, marginal=TRUE)
  test.for.zero( look2, rep(marginal.var,10) ,
                tag="normalization based on logic in LKrig.cov")
# check full covariance matrix
  look3<- LKrig.cov( x1, x1,LKinfo= LKinfo)
  test.for.zero( look3, look, tag="full covariance from matrix expressions")
# Now test w/o normalization
  LKinfo$normalize<- FALSE
  X<- LKrig.basis(x1, LKinfo)
  X<- spam2full(X)
# check on normalization to unit variance at each level
  Q<- LKrig.precision( LKinfo)
  Q<- spam2full( Q)
  look<- (X)%*% solve( Q) %*%t(X)
  look3<- LKrig.cov( x1, x1,LKinfo= LKinfo)
  test.for.zero( look3, look,
             tag="full covariance from matrix expressions w/o normalization")
#
# check of component covariance matrices
   LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=5,
                        a.wght=c(5,6,7), alpha=c(4,2,1), edge=FALSE)
  set.seed(123)
  x1<- cbind( runif( 10), runif(10))
  x2<- cbind(0,0)
  comp<- matrix( NA,10, 3)   
  for ( l in 1:3){
    grid.info<- LKinfo$grid.info
    grid.info$delta<- LKinfo$delta[l]
    LKinfo.temp<- LKrig.setup( grid.info=grid.info,
                         nlevel=1, a.wght=LKinfo$a.wght[l],
                         alpha=1, edge=FALSE) 
    comp[,l]<- LKrig.cov(x1,x2,LKinfo.temp )
  }
  look1<- comp%*%c(unlist( LKinfo$alpha))
  look3<- LKrig.cov( x1, x2,LKinfo= LKinfo)
  test.for.zero( look1, look3, tag="comp normalized cov and LKrig.cov")
#
# check construction with spatial a.wght
  cat("Now check spatial a.wght and alpha", fill=TRUE)

  LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=4,
                        alpha=1, edge=FALSE)
  a.wght<-  list( matrix(4 + (1:LKinfo$m)*.1, LKinfo$mx[1], LKinfo$my[1]))
  LKinfo<- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=1, NC=5,
                        a.wght=a.wght,
                        alpha=1, edge=FALSE)

  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- spam2full( look)
  test.for.zero( diag( look2), a.wght[[1]], tag="spatial a.wght 1 level")
# three levels
  LKinfo0 <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=c(5,5,5),
                        alpha=c(1,1,1), edge=FALSE)
  N<- LKinfo$mx*LKinfo0$my
  a.wght<-  list(
                 matrix(4 +  (1:N[1])*.1, LKinfo0$mx[1], LKinfo0$my[1] ),
                 matrix(4 +  (1:N[2])*.1, LKinfo0$mx[2], LKinfo0$my[2] ),
                 matrix(4 +  (1:N[3])*.1, LKinfo0$mx[3], LKinfo0$my[3] )
                )
   LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,
                        a.wght=a.wght,
                        alpha=c(1,1,1), edge=FALSE)
  look<- LKrig.precision( LKinfo=LKinfo, return.B=TRUE)
  look2<- spam2full( look)
  test.for.zero( diag( look2), unlist(a.wght), tag="spatial a.wght 3 levels")
#

# three levels nonzero buffer
  LKinfo0 <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4,NC.buffer=3,
                        a.wght=c(5,5,5),
                        alpha=c(1,1,1), edge=FALSE)
  N<- LKinfo0$mx*LKinfo0$my
  alpha<-  list(  (1:N[1])*.1, (1:N[2])*.1, (1:N[3])*.1)
  LKinfo <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4, NC.buffer=3,
                        a.wght=5,
                        alpha=alpha, edge=FALSE)
  look<- LKrig.precision( LKinfo, return.B=TRUE)
  look2<- spam2full( look)
 
  LKinfo2 <- LKrig.setup( cbind( c( -1,1), c( -1,1) ), nlevel=3, NC=4, NC.buffer=3,
                        a.wght=5,
                        alpha=c(1,1,1), edge=FALSE)
  look3<- spam2full(LKrig.precision( LKinfo2, return.B=TRUE))
  look3<-  diag( 1/sqrt(unlist(alpha)))%*%look3
  test.for.zero( look3, look2, tag=" 3 levels spatial alpha buffer=0")
  look4<-  spam2full(LKrig.precision( LKinfo))
  test.for.zero( t(look3)%*%look3, look4, tag="3 levels spatial alpha Q buffer=3")


# test of buffer grid points.
 NC.buffer<- 7
  LKinfo0 <- LKrig.setup( cbind( c( -1,1), c( 0,3) ), nlevel=3, NC=12,
                        a.wght=c(5,6,7), alpha=c(4,2,1), NC.buffer=0)
  LKinfo7 <- LKrig.setup( cbind( c( -1,1), c( 0,3) ), nlevel=3, NC=12,
                        a.wght=c(5,6,7), alpha=c(4,2,1), NC.buffer=NC.buffer)
  check.sum <- 0
  for(k in 1:3){
    check.sum <- check.sum + length((LKinfo0$grid[[k]])$x) +2*NC.buffer - length((LKinfo7$grid[[k]])$x)
    check.sum <- check.sum + length((LKinfo0$grid[[k]])$y) +2*NC.buffer - length((LKinfo7$grid[[k]])$y)
    
  }
  test.for.zero( check.sum, 0, relative=FALSE, tag="adding  7 buffer points")
# getting the margins and spatial domain right
# this test works because x y aspects are both divided evenly by delta.
   
 check.sum<-0
 for(k in 1:3){
    check.sum <- check.sum + (LKinfo0$grid[[k]])$x[1] - (LKinfo7$grid[[k]])$x[NC.buffer +1] 
    check.sum <- check.sum +  (LKinfo0$grid[[k]])$y[1] - (LKinfo7$grid[[k]])$y[NC.buffer +1]
    m2<- LKinfo7$mx[k]
    n2<- LKinfo7$my[k]
    check.sum <- check.sum + max((LKinfo0$grid[[k]])$x) - (LKinfo7$grid[[k]])$x[m2 - NC.buffer] 
    check.sum <- check.sum + max((LKinfo0$grid[[k]])$y) - (LKinfo7$grid[[k]])$y[n2 -NC.buffer] 
  }
 test.for.zero( check.sum, 0, relative=FALSE,tol=1e-8, tag="correct nesting of buffer=0  grid")
 cat("Done testing LKrig.precision",fill=TRUE)
  options( echo=FALSE)












