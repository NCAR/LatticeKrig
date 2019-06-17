subroutine LKDiag(entries, nEntries, diags, nRow, nCol, matrix)
  implicit none
  integer nEntries, nRow, nCol, diags(nEntries)
  double precision entries(nEntries), matrix(nRow, nCol)
  double precision currentEntry
  integer idx, currentDiag, diagLength, diagElement

  do idx = 1, nEntries
     currentEntry = entries(idx)
     currentDiag = diags(idx)
     if (currentDiag < 0) then
       diagLength = min(nRow + currentDiag, nCol)
       do diagElement = 1,diagLength
         matrix(diagElement-currentDiag, diagElement) = currentEntry
       enddo
     else
       diagLength = min(nCol - currentDiag, nRow)
       do diagElement = 1,diagLength
         matrix(diagElement, diagElement+currentDiag) = currentEntry
       enddo
     endif
   enddo
end 
