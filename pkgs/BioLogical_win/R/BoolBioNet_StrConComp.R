#' @title Caluclate the network's strongly connected components
#' @description This function provides a brief analysis of the strongly connected 
#' components to facilitate comparison in core dynamic analysis.
#' @param OneNet A four-element formatted list that records information of a Boolean network.
#' @return A list in which each element represents a strongly connected component.
#' Where those labeled as \code{_TRUE} in the name are non-trivial SCCs, while 
#' \code{_FALSE} indicates terminal or relay nodes.
#' @export
#' @examples
#' # Analyze the strongly connected components of a network.
#' # Here is an example of the network [c_2202] from the CellCollective Set
#' # BoolBioNet_StrConComp(BoolGRN_CellCollective$c_2202)
#' # Return a List
#' # [[1]]scc_1_FALSE   <-- a terminal node
#' # [1] "Exocytosis"
#' # [[2]]scc_2_FALSE   <-- a relay node 
#' # [1] "Packaging_Proteins"
#' # [[3]]scc_3_TRUE <-- a "real" SCC that meets the requirement of signal feedback.
#' # [1] "Protein_Phosphatase_1" "DARPP32"               "Calcineurin"          
#' # [4] "Calcium"               "Glutamate_Receptor"
#' # ... ... 
#' 
#' BoolBioNet_StrConComp(BoolGRN_CellCollective$c_2202)
#' 
BoolBioNet_StrConComp<-function(OneNet){
  Size=length(OneNet[[1]]);
  ind=otd=rep(0L,Size);
  con.vals=rep(666L,Size);# Here is real controlled values.
  NonIputer=NULL;# Label recording.
  for(ii in c(1:Size)){
    if(is.na(OneNet[[2]][[ii]][1])||is.null(OneNet[[2]][[ii]][1])){
      ind[ii]=0L;
      NonIputer=c(NonIputer,ii);}# Insert Exponents
    else {
      ind[ii]=length(OneNet[[2]][[ii]]);}
    if(is.na(OneNet[[3]][[ii]][1])||is.null(OneNet[[3]][[ii]][1])){
      otd[ii]=0L;}
    else {
      otd[ii]=length(OneNet[[3]][[ii]]);}}
  # Conver to node name.
  tmp=c_StrongConnectComponent(OneNet,ind,otd);
  tmp1=tmp[[1]]; tmp2=tmp[[2]];
  for(ii in c(1:length(tmp1))){
    tmp1[[ii]]=OneNet[[1]][tmp1[[ii]]+1];
  }
  tmp2=as.character(factor(tmp2, levels=0:1, labels=c("FALSE","TRUE")));
  names(tmp1)=paste0("scc_",c(1:length(tmp1)),"_",tmp2);
  return (tmp1);
}