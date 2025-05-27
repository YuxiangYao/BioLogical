#' @title Calculate the sensitivity of a Boolean function
#' @description The analysis depends on the enumeration of input vectors and 
#' observe its sensitivity of one-bit perturbations, such as f(..., x_i=0, ...)
#' is identical to f(..., x_i=1, ...) or not in all possible inputs. Detail 
#' concepts can see the paper [\href{https://doi.org/10.1103/PhysRevLett.93.048701}{Shmulevich,2004}].
#' @param Bit_Vec a bool (or 0/1-integral) vector, represents a truth table of Boolean function.
#' @return double, sensitivity of a Boolean function
#' @export
#' @examples
#' # Calculate the sensitivity of a Boolean function \code{c(1,1,0,1,1,1,0,0)}.
#' # BoolFun_Sensitivity(c(1,1,0,1,1,1,0,0))
#' # Return 1.25000
#' BoolFun_Sensitivity(c(1,1,0,1,1,1,0,0))
#' 
BoolFun_Sensitivity<-function(Bit_Vec){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  xx=c_BF_Sensitivity(tmp$bits,tmp$kin);
  return(xx);
}