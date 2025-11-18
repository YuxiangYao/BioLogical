#' @title Simplify networks
#' @description All isolated, redundant, or useless nodes are removed from the
#' system, particularly following the analysis of core dynamic components, 
#' resulting in a simplified network. 
#' @param OneNet A four-element formatted list that records information of a Boolean network.
#' @return A four-element formatted list that records information of a Boolean network.
#' @export
#' @examples 
#' # Analysis of a real Boolean genetic network. The terms "Analyzed_Net" and 
#' # "Simplified_Net" denote the network before and after simplification, respectively. 
#' # Compare the two to understand their differences and the purpose of the 
#' # simplification process.
#' set.seed(1234)
#' Analyzed_Net=BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_2084, AnalysisMod="relevant",
#'   ExtSelfLoop=FALSE, ResidualNet=TRUE)[[4]]
#' Simplified_Net=BoolBioNet_SimpNet(Analyzed_Net)
#' # Original network.
#' BoolBioNet_Visualize(BoolGRN_CellCollective$c_2084, "print")
#' # Anlayzed network.
#' BoolBioNet_Visualize(Analyzed_Net, "print")
#' # Simplified network.
#' BoolBioNet_Visualize(Simplified_Net, "print")
#' 
BoolBioNet_SimpNet<-function(OneNet){
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