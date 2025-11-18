#' @title Calculate the sensitivity of a multi-valued function
#' @description Correspodning concepts of multi-valued systems are generalized 
#' from the Boolean system. 
#' @param Map_Vec An IntegerVector that represents a mapping table of 
#' multi-valued function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system.
#' @return A numeric value representing the sensitivity of the given function.
#' @export
#' @examples
#' # Show the canalized structure of a multi-valued nested canalized function.
#' # MulVFun_Sensitivity(c(0,1,2, 1,0,1, 2,2,0), k=2, L=3);
#' # Return 1.777778
#' 
#' MulVFun_Sensitivity(c(0,1,2, 1,0,1, 2,2,0), k=2, L=3)
#' 
MulVFun_Sensitivity<-function(Map_Vec, k, L){
  tmp=r_CheckValidMulVFun(Map_Vec,k,L);
  xx=c_MulVF_Sensitivity(tmp,as.integer(k),as.integer(L),as.integer(L^k));
  return(xx);
}