EnvelopePlot <- function(x1, y1, x2 = x1, y2, col = NA, lineCol = NA, ...) {
  if(is.na(col)) {
    col = "thistle1"
  }
  if(is.na(lineCol)) {
    lineCol = "thistle3"
  }
  
  polygon(c(x1, rev(x2)), c(y1, rev(y2)), col = col, border = NA, ...)
  lines(x1, y1, lwd = 3, col = lineCol)
  lines(x2, y2, lwd = 3, col = lineCol)
}
