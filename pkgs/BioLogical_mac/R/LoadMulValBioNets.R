#' @title Load a multi-valued genetic network
#' @description Load genetic network information from external files.
#' @param NetName A string that denotes a valid file path or name for the target 
#' network.
#' @return A four-element formatted list that records information of a multi-valued
#' network. Note that the 4-th elements save weights rather than the Boolean funciton.
#' @export
#'
LoadMulValBioNets<-function(NetName){
  net1=BioLogical::LoadBoolBioNetwork(NetName,"Threshold",ThresMode = "one");
  net2=read.table(NetName, sep="\t", header=TRUE);
  for(genes in net1$AllMember){
    net2[net2[,2]==genes,3]->xxx;
    n.p=sum(xxx==1);# How many positive
    n.n=sum(xxx==2);# How many negative 
    xxx[xxx==2]=-1;
    #is.all.n=(n.n==length(xxx))# Is all negative
    # if(n.p>0)xxx[xxx>0]=xxx[xxx>0]/n.p;
    # if(n.n>0)xxx[xxx<0]=xxx[xxx<0]/n.n;
    if(length(xxx)>0){
      net1$BoolFun[[genes]]=xxx;
    } else {
      net1$BoolFun[[genes]]=0;
    }
  }
  return(net1);
}
