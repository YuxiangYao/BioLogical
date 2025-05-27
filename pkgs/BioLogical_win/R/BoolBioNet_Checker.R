#' @title Check the Boolean network is appropriate
#' @description It compares and analyzes the lengths of node's names, input/output links, and Boolean function configuration. 
#' Check the Boolan network is structurally appropriate. For example, one node has three predecessors with 4-bit truth table is invalid.
#' @param OneBioNet a 4-element formatted list, record Boolean network information
#' @return a logical value, \code{TRUE}, for other downstream analysis. 
#' Any abnormal results and invalid inputs return error messages.
#' @export 
#' @examples
#' # Analysis of one real Boolean genetic network.
#' # BoolBioNet_Checker(BoolGRN_CellCollective$c_1557);
#' # Return TRUE (Denote all lists are appropriate)
#' BoolBioNet_Checker(BoolGRN_CellCollective$c_1557)
BoolBioNet_Checker<-function(OneBioNet){
  # Boolean networks not distinguish biologically genetic or manually generated systems.
  # RealBioNet should be standard four-list format data that records information of bool-nets.
  if(4>length(OneBioNet)){
    stop("Not enough liste provided.\n");}
  genname=OneBioNet[["AllMember"]]; # Gene's names
  inDegs=OneBioNet[["InEdge"]]; # In-degree frames.
  otDegs=OneBioNet[["OutEdge"]]; # Out-degree frames.
  bnlist=OneBioNet[["BoolFun"]]; # Boolean function list.
  # Check point-1 (Same lengths and names of in-degree, out-degree, & boolfunc list-slots).
  # Should be FALSE
  LogFlag1=is.null(genname)||is.null(inDegs)||is.null(otDegs)||is.null(bnlist);
  # Should be TRUE.
  LogFlag2=(length(genname)==length(inDegs))&&(length(genname)==length(otDegs))&&(length(genname)==length(bnlist));
  # Should be TRUE. (All names should be same and mathced!)
  LogFlag3=identical(genname,names(inDegs))&&identical(genname,names(otDegs))&&identical(genname,names(bnlist));
  if(LogFlag1||(!LogFlag2)||(!LogFlag3)){
    stop("Missing key information list; or diffferent lengths of each list. 
      Please check your \"OneBioNet\" object.\n");}
  n.gene=length(OneBioNet[[1]]);
  # Check point-2 (appropriate input and Boolean functions.)
  for(ii in c(1:n.gene)){
    if(is.na(inDegs[[ii]][1])){# The genes are external factor (their in-degrees are zero).
      if(!is.na(bnlist[[ii]][1])){# Corresponding lists should be NA.
        stop("External factors should be NA.\n");}}
    else {# Have input.
      # The length should be in-degree[ii] equals to log2(boolfun[ii]);
      if(length(inDegs[[ii]])!=as.integer(log2(length(bnlist[[ii]])))){
        stop("Mismatching lengths of node: ",genname[ii],"\n");}
      # Check rationality of Boolean functions.
      if(length(inDegs[[ii]])>0){
        LogFlag1=!(sum(sapply(bnlist[[ii]],is.integer))==length(bnlist[[ii]]));# All integrate?
        LogFlag2=!(sum(sapply(bnlist[[ii]],is.logical))==length(bnlist[[ii]]));# All logical?
        if(!(LogFlag1||LogFlag2)){
          stop("Some invalid elements exists in Boolean function of ",genname[ii],"\n");}}}}
  xx=TRUE;# Return logical label.  
  return (xx);
}
