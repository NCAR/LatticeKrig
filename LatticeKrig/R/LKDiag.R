LKDiag <- function(entries, nrow, diags = NULL, ncol = nrow) {
  entries = as.double(entries)
  nEntries = as.integer(length(entries))
  nrow = as.integer(nrow)
  ncol = as.integer(ncol)
  mat <- as.double(matrix(0, nrow=nrow, ncol=ncol))
  
  if (is.null(diags)) {
    if (nEntries %% 2 == 0) {
      diags = c((-nEntries/2) : -1, 1 : (nEntries/2))
    } else {
      diags = (-(nEntries-1)/2) : ((nEntries-1)/2)
    }
  }
  
  diags = as.integer(diags)
  if (length(entries == 1)) {
    entries = as.double(rep(entries[1], length(diags)))
    nEntries = as.integer(length(entries))
  }
  
  if (length(entries) != length(diags)) {
    stop("The length of entries and length of diagonals don't match")
  }
  out <- .Fortran("LKDiag", entries = entries, nEntries = nEntries, diags = diags,
                  nRow = nrow, nCol = ncol, matrix = mat, package="LatticeKrig")
  return(matrix(out$matrix, nrow = nrow, ncol = ncol))
}
