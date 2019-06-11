LKDiag <- function(entries, nrow, diags = NA, ncol = nrow) {
  nEntries = length(entries)
  if (is.na(diags)) {
    if (nEntries %% 2 == 0) {
      diags = c((-nEntries/2) : -1, 1 : (nEntries/2))
    } else {
      diags = (-(nEntries-1)/2) : ((nEntries-1)/2)
    }
  }
  mat <- matrix(0, nrow=nrow, ncol=ncol)
  out <- .Fortran("LKDiag", entries = entries, nEntries = nEntries, diags = diags, nRow = nrow, nCol = ncol, matrix = mat)
  return(out$matrix)
}