#' @title Generate the header of mapping table
#' @description Generate state enumeration table for the given number of input 
#' variables (multi-valued version)
#' @param inDegree An integer that denotes the number of input variables (Default: 2).
#' @param Lsystem An integer that represents the level of multi-valued system (Default: 2).
#' @return a data.frame of input case enumeration 
#' @export
r_GenerateTruthTable_M<-function(inDegree=2L, Lsystem=2L){
  #if(as.integer(inDegree)<=0||as.integer(Lsystem^inDegree)>10000){
  #  stop("Invalid in-degre or too long table! Please select appropriate numbers of in-degree and L-level.\n");}
  tabs=matrix(NA,as.integer(Lsystem^inDegree),inDegree);
  colnames(tabs)=paste0("v",(c(inDegree:1)-1));
  l_sys=c(1:Lsystem)-1;
  for(ii in c(1:inDegree)){
    tabs[,ii]=rep(
      rep(l_sys,each=as.integer(Lsystem^(inDegree-ii))),
      time=as.integer(Lsystem^(ii-1)));}
  tabs=as.matrix(tabs); #tabs=as.data.frame(tabs);
  return(tabs);
}