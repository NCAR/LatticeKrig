EnvelopePlot <- function(x1, y1, x2, y2, col = NA, ...) {
  if(is.na(col)) {
    col = "thistle1"
  }
  lines(x1, y1, lwd = 5, col = "thistle3")
  lines(x2, y2, lwd = 5, col = "thistle3")
  polygon(c(x1, rev(x2)), c(y1, rev(y2)), col = col, border = NA, ...)
}
