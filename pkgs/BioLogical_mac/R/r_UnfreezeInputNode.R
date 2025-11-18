#' @title Add self-loop circuits to all external nodes
#' @description Convert all input and external nodes into self-loop nodes
#' @param RealBioNet A four-element formatted list that records information of a
#' Boolean network (See \link{BoolGRN_CellCollective})
#' @return A four-element formatted list that records information of the modified
#' Boolean network (See \link{BoolGRN_CellCollective})
#' @export
r_UnfreezeInputNode<-function(RealBioNet){
  nodes=RealBioNet[[1]];
  indgs=RealBioNet[[2]];
  otdgs=RealBioNet[[3]];
  bbnns=RealBioNet[[4]];
  logis=is.na(indgs);
  for(ii in c(1:length(nodes))){
    if(logis[ii]){
      indgs[[ii]]=c(as.integer(ii-1));
      otdgs[[ii]]=c(otdgs[[ii]],as.integer(ii-1));
      bbnns[[ii]]=c(0L,1L);# Default: positive regualtions.
    }
  }
  res=list(nodes,indgs,otdgs,bbnns);
  names(res)=names(RealBioNet);
  return(res);
}