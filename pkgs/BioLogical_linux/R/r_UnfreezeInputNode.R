#' Treat all input/external nodes as self-loop node.
#' @param RealBioNet List[[4]], a formatted data, record biological network information
#' @return List[[4]], a formatted data, record biological network information
#' @export
r_UnfreezeInputNode<-function(RealBioNet){
  nodes=RealBioNet[[1]];
  indgs=RealBioNet[[2]];
  otdgs=RealBioNet[[3]];
  bbnns=RealBioNet[[4]];
  logis=is.na(indgs);
  for(ii in c(1:length(nodes))){
    if(logis[ii]){
      indgs[[ii]]=c(as.integer(ii-1));# 
      otdgs[[ii]]=c(otdgs[[ii]],as.integer(ii-1));
      bbnns[[ii]]=c(0L,1L);
    }
  }
  res=list(nodes,indgs,otdgs,bbnns);
  names(res)=names(RealBioNet);
  return(res);
}