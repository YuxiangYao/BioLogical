#' @title Display the Quine-McCluskey form of a Boolean function
#' @description The function employs our native C++ standard Quine-McCluskey algorithm
#' to obtain prime implicants.
#' The function utilizes our native C++ implementation of the Quine-McCluskey
#' algorithm to compute the prime implicants.
#' @param Bit_Vec A Boolean vector (or 0/1-valued \code{IntegerVector}) represents the
#' truth table of a Boolean function.
#' @param VarsName A \code{CharacterVector} that represents each name of variable.
#' @details The function prints the Quine-McCluskey form in the terminal, using
#' uppercase and lowercase letters to indicate whether a variable is 1 or 0,
#' respectively. Form low to high bits, variables are represented as \code{A/a},
#' \code{B/b}, \code{C/c}, ... If \code{VarsName} is provided, it must be of
#' sufficient length to represent the bits from low to high; otherwise, it
#' returns an error message and execution halts. Where the symbol "\code{~}"
#' denotes the logical NOT.
#' @return Show the Quine-McCluskey form of the given Boolean function.
#' @export
#' @examples
#' # Calculate the Quine-McCluskey form of c(0,0,0,1,1,0,1,1).
#' # (a=0 and c=1) or (a=1 and b=1) ==> 1; (a->c, from low to high bit)
#' # BoolFun_QMForm(c(0,0,0,1,1,0,1,1));
#' # Return { 1 = aC + AB }
#' BoolFun_QMForm(c(0,0,0,1,1,0,1,1))
#'
#' # Calculate the Quine-McCluskey form of c(1,1,0,1,1,1,0,0).
#' # (a=1 and c=0) or (b=0) ==> 1; (a->c, from low to high bit)
#' # BoolFun_QMForm(c(1,1,0,1,1,1,0,0), c("X1", "X2", "X3"));
#' # Return { 1= X1*~X3 + ~X2} (identical to { 1 = Ac + b })
#' BoolFun_QMForm(c(1,1,0,1,1,1,0,0), c("X1", "X2", "X3"))
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
