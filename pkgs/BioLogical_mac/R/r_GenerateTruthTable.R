#' @title Generate the header of truth table
#' @description Generate state enumeration table for the given number of input variables
#' @param inDegree An integer that denotes the number of input variables.
#' @param SpinLike A logical value. Should output the spin-like form?
#' @return A data.frame of input case enumeration 
#' @export
r_GenerateTruthTable<-function(inDegree=2L, SpinLike=FALSE){
  if(as.integer(inDegree)<=0||as.integer(inDegree)>16){
    stop("Invalid in-degree. Please select an integrate between [1,16].\n");}
  lens=bitwShiftL(1,inDegree)-1;
  tabs=NULL;
  for(ii in c(0:lens)){
    tmp=as.integer(as.vector(intToBits(ii)));
    tabs=rbind(tabs,tmp[inDegree:1]);}
  colnames(tabs)=paste0("v",(c(inDegree:1)-1));
  if(SpinLike){# Converted to spin-like
    tabs=2*tabs-1;}
  return(tabs);
}