#' Determines whether the input vector is a valid mapping result.
#' @param a_Vec a integral/logical vector, represents a mapping table of multi-valued function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system.
#' @return IntegerVector
#' @export
r_CheckValidMulVFun<-function(a_Vec, k, L){
  if(!is.numeric(a_Vec)&&!is.integer(a_Vec)&&!is.logical(a_Vec)){
    stop("Invalid mapping table input. The table should be numeric, integer or logical.\n");}
  else {
    a_vec=as.integer(a_Vec);}
  tmp=log(length(a_vec),base=L);
  if(abs(k-tmp)>1e-7){
    stop("Invalid length of a_Vec. The length should be L^k.\n");}
  return(a_vec);
}