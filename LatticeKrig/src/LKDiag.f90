subroutine LKDiag(entries, nEntries, diags, nrow, ncol, matrix)
  integer nEntries, nRow, nCol
  double precision entries(nEntries), diags(nEntries), matrix(nRow, nCol)
  integer idx, currentEntry, currentDiag, diagLength
  
  do idx = 1, nEntries
    currentEntry = entries(idx)
    currentDiag = diags(idx)
    if (currentDiag < 0) then
      diagLength = min(nRow + currentDiag, nCol)
      matrix(1-currentDiag : diagLength-currentDiag, 1:diagLength) = currentEntry
    else
      diagLength = min(nCol - currentDiag, nRow)
      matrix(1:diagLength, 1+currentDiag : diagLength+currentDiag) = currentEntry
    endif
  enddo
end subroutine LKDiag
