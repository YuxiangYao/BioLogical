#' @title Calculate the complexity of a multi-valued function
#' @description The analysis depends on average of prime implicants in multi-
#' valued function and its complementary one(s), based on multi-valued system's 
#' Quine-McCluskey method. Although this function is downward compatible with 
#' Boolean scenarios, here recommends to directly use Boolean corresponding 
#' one (\link{BoolFun_Complexity}) because of latter's performance. 
#' @param aVec a vector, represents a mapping table of multi-valued function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system (Default: 2, Boolean).
#' @return double, complexity of a multi-valued function.
#' @export
#' @examples
#' # Calculate the complexity of a Booelan function.
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
MulVFun_Complexity<-function(aVec, k, L=2){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  # if(all(a_vec == a_vec[1])){
  #   xx=0;
  # } else {
  xx=c_MulF_Complexity(a_vec,as.integer(k),as.integer(L));
  #}
  return(xx);
}