#' @title Determines whether the input vector is a valid mapping result.
#' @param Bit_Vec a bool vector, represents a truth table of Boolean function.
#' @return a List, [[1]] Boolean vector, [[2]] in-degree
#' @export
r_CheckValidBoolFun<-function(Bit_Vec){
  if(!is.numeric(Bit_Vec)&&!is.integer(Bit_Vec)&&!is.logical(Bit_Vec)){
    stop("Invalid mapping table input. The table should be numeric, integer or logical.\n");}
  else {
    bit_vec=as.logical(Bit_Vec);}
  if(length(bit_vec)<=0||log2(length(bit_vec))>16){
    stop("Invalid length of bit_vec. The number of variables is less than 16.\n");}
  tmp=log2(length(bit_vec));
  if(abs(tmp-round(tmp))>1e-7){
    stop("Invalid length of bit_vec. The length should be 2^k.\n");}
  res=list(bits=bit_vec,kin=round(tmp));
  return(res);
}