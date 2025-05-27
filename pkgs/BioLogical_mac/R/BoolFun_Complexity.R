#' @title Calculate the complexity of a Boolean function
#' @description The analysis depends on average of prime implicants in Boolean 
#' function and its complementary one, based on Boolean Quine-McCluskey method.
#' @param Bit_Vec a bool (or 0/1-integral) vector, represents a truth table of Boolean function.
#' @return double, complexity of a Boolean function
#' @export
#' @examples
#' # Calculate the complexity of a Boolean function \code{c(0,0,1,0,0,1,1,0)}.
#' # BoolFun_Complexity(c(0,0,1,0,0,1,1,0))
#' # {(*10)+(101) ==> 1 (00*)+(*00)+(*11) ==> 0 } -> 2.5
#' # Return 5
#' BoolFun_Complexity(c(0,0,1,0,0,1,1,0))
BoolFun_Complexity<-function(Bit_Vec){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  if(all(tmp$bits == tmp$bits[1])){# All elements are same.
    xx=0;
  } else {
    xx=c_BF_Complexity(tmp$bits,tmp$kin);
  }
  return(xx);
}