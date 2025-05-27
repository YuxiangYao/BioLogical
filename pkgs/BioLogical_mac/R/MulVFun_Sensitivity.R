#' @title Calculate the sensitivity of a multi-valued function
#' @description Correspodning concepts of multi-valued systems are generalized 
#' from the Boolean system. 
#' @param Map_Vec a multi-valued vector, represents a mapping table of multi-
#' valued functions.
#' @param k An integer for the number of input-variable.
#' @param L An integer for the level of discrete system.
#' @return double, sensitivity of a multi-valued function.
#' @export
#' @examples
#' # Show the canalized structure of a multi-valued nested canalized function.
#' # MulVFun_Sensitivity(c(0,1,2, 1,0,1, 2,2,0), k=2, L=3)
#' # Return 1.777778
#' MulVFun_Sensitivity(c(0,1,2, 1,0,1, 2,2,0), k=2, L=3)
#' 
MulVFun_Sensitivity<-function(Map_Vec, k, L){
  tmp=r_CheckValidMulVFun(Map_Vec,k,L);
  xx=c_MulVF_Sensitivity(tmp,as.integer(k),as.integer(L),as.integer(L^k));
  return(xx);
}