#' @title Simplify networks
#' @description All isolated/useless/ineffective nodes are removed from the systems 
#' (Especially after analyzing). Return a simplified network. 
#' @param OneNet a 4-element formatted list, record biological network information
#' @return a 4-element formatted list, record the simplified network
#' @export
#' @examples
#' # Analysis of one real Boolean genetic network. The following 
#' # "Analyzed_Net" and "Simplified_Net" respectively represents the simplified 
#' # or not after analyzing process. Compare them and understand the usage.
#' set.seed(1234)
#' Analyzed_Net=BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_2084, AnalysisMod="scaling",
#'   ExternalSelfLoop=FALSE, ResidualNet = TRUE)[[4]];
#' Simplified_Net=BoolBioNet_SimplifiedNet(Analyzed_Net)
#' # Original network.
#' BoolBioNet_Visualization(BoolGRN_CellCollective$c_2084, "print")
#' # Anlayzed network.
#' BoolBioNet_Visualization(Analyzed_Net, "print")
#' # Simplified network.
#' BoolBioNet_Visualization(Simplified_Net, "print")
#' 
BoolBioNet_SimplifiedNet<-function(OneNet){
  # Split object as four lists.
  member=OneNet[[1]];
  indegs=OneNet[[2]];
  otdegs=OneNet[[3]];
  boolfn=OneNet[[4]];
  # Inner function of remove pointed node.
  NodeRemove<-function(vx,id){
    if(is.na(vx[1])){
      return(NA);}
    else {
      vv=vx; index=(vx>=id);
      vv[index]=vv[index]-1;return(vv);}}
  # Execute repeatedly >>>
  ii=1;
  while(ii<=length(member)){
    if(is.na(indegs[[ii]][1])&&is.na(otdegs[[ii]][1])){# Isolated node
      indegs=sapply(indegs,NodeRemove,id=ii-1);
      indegs[[ii]]=NULL;
      otdegs=sapply(otdegs,NodeRemove,id=ii-1);
      otdegs[[ii]]=NULL;
      boolfn[[ii]]=NULL;
      member=member[-ii];
    } else {# Check the next 
      ii=ii+1;}}
  SimNets=list(
    AllMember=member,
    InEdge=indegs,
    OutEdge=otdegs,
    BoolFun=boolfn);
  return(SimNets);
}