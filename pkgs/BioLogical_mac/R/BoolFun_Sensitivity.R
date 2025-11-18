#' @title Calculate the sensitivity of a Boolean function
#' @description The analysis depends on enumerating input vectors and examining 
#' the mapping result of one-bit perturbations, such as \eqn{f(..., x_i=0, ...)} is 
#' identical to \eqn{f(..., x_i=1, ...)} or not. Detail concepts can see the paper 
#' [\href{https://doi.org/10.1103/PhysRevLett.93.048701}{Shmulevich2004}].
#' @param Bit_Vec A Boolean vector (or 0/1-valued \code{IntegerVector}) represents the 
#' truth table of a Boolean function.
#' @return A numeric value that denotes the sensitivity of a Boolean function.
#' @export
#' @examples
#' # Calculate the sensitivity of a Boolean function c(1,1,0,1,1,1,0,0).
#' # BoolFun_Sensitivity(c(1,1,0,1,1,1,0,0));
#' # Return 1.25000
#' BoolFun_Sensitivity(c(1,1,0,1,1,1,0,0))
#' 
BoolFun_Sensitivity<-function(Bit_Vec){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  xx=c_BF_Sensitivity(tmp$bits,tmp$kin);
  return(xx);
}