LKDefaultVarNames<- function(A, tag=NULL){
  if( is.null(A)){
    return(NULL)
  }
  if( ! is.matrix(A)){
    stop("Not a matrix! ")
  }
  colA<- colnames(A)
  if( is.null(colA)){
    if( is.null(tag)){ 
      stop("need a tag for column names")
    }
    colA<- paste0(tag,1:ncol(A))
  }
  return(colA)
}
