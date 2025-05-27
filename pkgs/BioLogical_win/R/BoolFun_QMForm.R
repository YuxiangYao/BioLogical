#' @title Show a Boolean function's Quine-McCluskey form
#' @description The function employs our native C++ standard Quine-McCluskey algorithm 
#' to obtain prime implicants. Because of algorithms and system properties, the
#' result is not unique.
#' @param Bit_Vec a bool (or 0/1-integral) vector, represents a truth table of Boolean function.
#' @param VarsName a character vector, represents each variable's name
#' @details Function prints out the Quine-McCluskey form in the terminal, where 
#' upercase/lowercase denotes the corresponding variable is one or zero. From low 
#' to high bits are shown as \code{A/a}, \code{B/b}, \code{C/c}, ... (\code{is.null(VarsName)} 
#' is \code{TRUE}). If provided (enough length or throw an error), the character 
#' vector represents the low to high bits. Symbols in each clause denote zero values.
#' @return Show Quine-McCluskey forms.
#' @export
#' @examples
#' # Calculate the Quine-McCluskey form of a Boolean function \code{c(0,0,0,1,1,0,1,1)}.
#' # (a=0 and c=1) or (a=1 and b=1) ==> 1 ; (a->c, from low to high digit)
#' # BoolFun_QMForm(c(0,0,0,1,1,0,1,1))
#' # Return { 1 = aC + AB }
#' BoolFun_QMForm(c(0,0,0,1,1,0,1,1))
#' 
#' # Calculate the Quine-McCluskey form of a Boolean function \code{c(1,1,0,1,1,1,0,0)}.
#' # (a=1 and c=0) or (b=0) ==> 1 ; (a->c, from low to high digit)
#' # BoolFun_QMForm(c(1,1,0,1,1,1,0,0),c("X1", "X2", "X3"))
#' # Return { 1= X1*~X3 + ~X2} (identical to { 1 = Ac + b })
#' BoolFun_QMForm(c(1,1,0,1,1,1,0,0),c("X1", "X2", "X3"))
#' 
BoolFun_QMForm<-function(Bit_Vec, VarsName=NULL){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  bit_vec=tmp$bits;
  k_in=tmp$kin;
  if(all(bit_vec==bit_vec[1])){
    if(1==bit_vec[1]){
      cat("It is a 1-constant Boolean function.\n");
    } else {
      cat("It is a 0-constant Boolean function.\n");}}
  else {
    var.names=c(NA,NA,"a");
    if(is.null(VarsName)){;}
    else {
      if(length(VarsName)<as.integer(k_in)){
        stop("Not enough variable's names.\n");
      }else {
        var.names=VarsName;}}
    c_BF_QuineMcCluskey(bit_vec,k_in,var.names);}
}