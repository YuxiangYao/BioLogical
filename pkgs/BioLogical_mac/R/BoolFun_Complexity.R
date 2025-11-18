#' @title Calculate the complexity of a Boolean function
#' @description The analysis depends on average of prime implicants in a Boolean 
#' function, based on Boolean Quine-McCluskey method.
#' @param Bit_Vec A Boolean vector (or 0/1-valued \code{IntegerVector}) represents
#' the truth table of a Boolean function.
#' @return A numeric value that denotes the complexity of a Boolean function.
#' @export
#' @examples
#' # Calculate the complexity of a Boolean function \code{c(0,0,1,0,0,1,1,0)}.
#' # Based on Quine-McCluskey method:
#' # (*10)+(101) ==> 1 
#' # (00*)+(*00)+(*11) ==> 0 
#' # Thus, 2+3 =5 items as the complexity.
#' # BoolFun_Complexity(c(0,0,1,0,0,1,1,0))
#' # Return 5
#' BoolFun_Complexity(c(0,0,1,0,0,1,1,0))
#' 
BoolFun_Complexity<-function(Bit_Vec){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  if(all(tmp$bits == tmp$bits[1])){# All elements are same.
    xx=0;
  } else {
    xx=c_BF_Complexity(tmp$bits,tmp$kin);
  }
  return(xx);
}