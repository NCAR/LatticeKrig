###############################################
# A Function to test mpi.apply with LatticeKrig
#
# Input:
#   numSlaves: number of mpi slaves to spawn
#   timeSteps: amount of ozone2 data to run
#
# Output:
#   a list of two times to determine speed of
#   mpi apply functions
###############################################

rmpiTest <- function(numSlaves, timeSteps){
	require(LatticeKrig)
	require(Rmpi)
	data(ozone2)

	LKtest <- function(index){
		# Get data and remove missing values
		x <- ozone2$lon.lat
		y <- ozone2$y[index,]
		good <- !is.na(y)
		x <- x[good,]
		y <- y[good]

		# Make the LaatticeKrig object and predict a surface (time intensive)
		obj <- LatticeKrig(x,y)
		out <- predictSurface(obj)

		return(out)
	}

	
	# Spawn n slaves
	ptm1 <- proc.time()
	mpi.spawn.Rslaves(nslaves=numSlaves)
	ptm2 <- proc.time()

	# send the data and LKtest to the slaves
	mpi.bcast.cmd(library(LatticeKrig))
	mpi.bcast.Robj2slave(ozone2)
	mpi.bcast.Robj2slave(LKtest)

	# run with applyLB or iapplyLB
	out <- mpi.applyLB(1:timeSteps,LKtest)
	outf <- lapply(out,as.array)
	tm2 <- proc.time() - ptm2
	mpi.close.Rslaves()
	tm1 <- proc.time() - ptm1

	return(list(tm1,tm2))
}