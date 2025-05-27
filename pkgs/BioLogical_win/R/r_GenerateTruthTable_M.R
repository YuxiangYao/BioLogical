#' Generate a truth table with pointed in-degree function.
#' @param inDegree an integrate value, represents a truth table of Boolean function.
#' @param Lsystem an integer, L-level multi-valued system.
#' @return a data.frame of input case enumeration 
#' @export
r_GenerateTruthTable_M<-function(inDegree=2L,Lsystem=2L){
  #if(as.integer(inDegree)<=0||as.integer(Lsystem^inDegree)>10000){
  #  stop("Invalid in-degre or too long table! Please select appropriate numbers of in-degree and L-level.\n");}
  tabs=matrix(NA,as.integer(Lsystem^inDegree),inDegree);
  colnames(tabs)=paste0("v",(c(inDegree:1)-1));
  l_sys=c(1:Lsystem)-1;
  for(ii in c(1:inDegree)){
    tabs[,ii]=rep(
      rep(l_sys,each=as.integer(Lsystem^(inDegree-ii))),
      time=as.integer(Lsystem^(ii-1)));}
  #tabs=as.data.frame(tabs);
  tabs=as.matrix(tabs);
  return(tabs);
}