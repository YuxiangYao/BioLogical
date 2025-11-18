#' @title Calculate the complexity of a multi-valued function
#' @description The analysis depends on average of prime implicants in multi-
#' valued function and its complement(s), based on Quine-McCluskey method for 
#' multi-valued systems. Although this function is downward compatible with 
#' Boolean scenarios, it is recommended to use Boolean-specific counterpart  
#' (\link{BoolFun_Complexity}) directly due to its superior performance.
#' @param aVec An IntegerVector that represents a mapping table of multi-valued function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system (Default: 2).
#' @return A numeric value representing the complexity of a multi-valued function.
#' @export
#' @examples
#' # Calculate the complexity of a Booelan function (Show the compatibility).
#' # MulVFun_Complexity(c(0,0,1,0, 0,1,1,0), k=3)
#' # {(*10)+(101) ==> 1 (00*)+(*00)+(*11) ==> 0 } -> (2+3)=5
#' # Return 5
#' MulVFun_Complexity(c(0,0,1,0, 0,1,1,0), 3) # Compatible with Boolean system
#' 
#' # Calculate the complexity of a Ternary function.
#' # MulVFun_Complexity(c(1,1,1, 0,2,1, 1,0,2), k=2, L=3)
#' # { [(10)|(21)] ==> 0; [(0*)|(12)|(20)] ==> 1; [(11)|(22)] ==> 2} -> 
#' # (2+3+2)=7
#' # Return 7
#' MulVFun_Complexity(c(1,1,1, 0,2,1, 1,0,2), k=2, L=3) # Ternary system
#' 
MulVFun_Complexity<-function(aVec, k, L=2){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  xx=c_MulF_Complexity(a_vec,as.integer(k),as.integer(L));
  return(xx);
}