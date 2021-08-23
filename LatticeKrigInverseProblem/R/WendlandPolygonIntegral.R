WendlandPolygonIntegral <- function(vertices, wendlandCenter = c(0,0), wendlandRadius = 1) {
  #translate and scale the points so the wendland function is centered at the origin with radius 1
  vertices <- sweep(vertices, 1, wendlandCenter, '-') / wendlandRadius
  N = dim(vertices)[2]
  xMin <- min(vertices[1,])
  xMax <- max(vertices[1,])
  yMin <- min(vertices[2,])
  yMax <- max(vertices[2,])
  #find all the y coordinates of vertices, which will be the y-coordinates of the lines we evaluate line integrals on
  yVert <- unique(sort(vertices[2,]))
  M = length(yVert)
  yGap = diff(yVert)
  segmentsOnLine <- rep(0L, M)
  intersections <- NULL
  #TODO and maybe also better line-polygon intersection, but maybe not because this is pretty trivial; could port it to FORTRAN though
  
  #this block finds all the intersection points of each horizontal line with the polygon
  for (i in 1:M) {
    yIntersections <- NULL
    for (j in 1:N-1) {
      if (vertices[2,j] == yVert[i]) {
        #first check makes sure that the line goes inside the polygon, not just tangent to the vertex;
        #second check makes sure that if the line overlaps with a side made of 3 or more points, only the outermost get added
        if (xor(vertices[2,(j-1+N)%%N] < vertices[2,j], vertices[2,(j+1)%%N] < vertices[2,j]) || xor(vertices[2,(j-1+N)%%N] == vertices[2,j], vertices[2,(j+1)%%N] == vertices[2,j])) {
          yIntersections <- cbind(yIntersections, vertices[,j])
        }
      } else if (vertices[2,j] < yVert[i] && vertices[2,(j+1)%%N] > yVert[i] || vertices[2,j] > yVert[i] && vertices[2,(j+1)%%N] < yVert[i]) {
        #the line segment from vertex j to j+1 crosses the y coordinate we're checking against
        t = (yVert[i] - vertices[2,j])/(vertices[2, (j+1)%%N] - vertices[2,j])
        x = vertices[1,j] + t*(vertices[1,(j+1)%%N] - vertices[1,j])
        yIntersections <- cbind(yIntersections, rbind(x, yVert[i]))
      }
    }
    #sort the yIntersections before adding them, so the line segments connect correctly
    ord <- order(yIntersections[1,])
    if (length(ord) %% 2 == 1) {
      #odd number of intersections; should be impossible for finite polygons
      stop(sprintf("Odd number of intersections with polygon at y=%.2f", yVert[i]))
    }
    intersections <- cbind(intersections, yIntersections[,ord])
    segmentsOnLine[i] = length(ord)/2
  }
  
  intersections <- matrix(intersections, nrow=4) #rearrange the list of intersections so that the two endpoints of each segment inside the polygon form one column
  segmentIntegrals <- lineIntegral(intersections)
  lineIntegrals = rep(0, M)
  segmentIdx = 0
  #in case there are lines with more than one line segment inside the polygon, combine those segments' integrals to get the line integral
  for (i in 1:M) {
    lineIntegrals[i] = sum(segmentIntegrals[(segmentIdx+1):(segmentIdx+segmentsOnLine[i])])
    segmentIdx = segmentIdx + segmentsOnLine[i]
  }
  #do 1D quadrature with the line integral values; trapezoid rule is used here to start
  totalArea = 0.5*lineIntegrals[1]*yGap[1]
  for (i in 2:(M-1)) {
    totalArea = totalArea + 0.5*lineIntegrals[i] * (yGap[i-1] + yGap[i])
  }
  totalArea = totalArea + 0.5*lineIntegrals[M] * yGap[M-1]
  return(totalArea * wendlandRadius^2) #scale back up by wendlandRadius, since we scaled down by it at the start
}

lineIntegral <- function(lines) {
  #placeholder line integral function; returns the lengths of the line (i.e. integrating a constant over the line) so we get an area approximation
  return(sqrt((lines[3,]-lines[1,])^2 + (lines[4,]-lines[2,])^2))
}