#' @title Multi-valued function checker
#' @description Check whether the input vector is a valid multi-valued function
#' @param a_Vec An IntegerVector representing a multi-valued function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system.
#' @return IntegerVector or prints an error message and execution halts.
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