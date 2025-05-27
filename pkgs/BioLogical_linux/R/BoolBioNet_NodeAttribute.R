#' @title Show a node's attribute
#' @description This function returns the predecessors, successors, and 
#' regulatory patterns of given nodes. 
#' @param RealBioNet a 4-element formatted list, record biological network information.
#' @param NodeName Node/Gene's name.
#' @return Attribute/information/feature of the node.
#' @export
#' @examples
#' # Print the information of Node "Akt" in gene network [c_11863] of 
#' # CellCollective Set.
#' # BoolBioNet_NodeAttribute(BoolGRN_CellCollective$c_11863,"Akt")
#' # Return information:
#' # $TruthTable
#' #      PI3K [Akt]
#' # [1,]    0     0
#' # [2,]    1     1
#' # $Successor
#' # [1] "IKK"  "mTOR"
#' BoolBioNet_NodeAttribute(BoolGRN_CellCollective$c_11863,"Akt")
BoolBioNet_NodeAttribute<-function(RealBioNet, NodeName){
  if(BoolBioNet_Checker(RealBioNet)){;}
  if(!(NodeName%in%RealBioNet[[1]])){
    stop("The input node does not exist!\n");
  } else {
    # Input:
    if(is.na(RealBioNet$InEdge[[NodeName]][1])||is.null(RealBioNet$InEdge[[NodeName]][1])){
      ttt="This is an external signal node.";
    } else {
      in_nod=RealBioNet$InEdge[[NodeName]];
      varnum=length(in_nod);
      n.leng=bitwShiftL(1,varnum);
      ttt=matrix(NA,n.leng,varnum);
      colnames(ttt)=RealBioNet[[1]][in_nod+1];
      for(ii in c(0:(n.leng-1))){
        ttt[ii+1,]=as.numeric(intToBits(ii)[varnum:1]);
      }
      ttt=cbind(ttt,RealBioNet$BoolFun[[NodeName]]);
      colnames(ttt)[varnum+1]=paste0("[",NodeName,"]");
    }
    if(is.na(RealBioNet$OutEdge[[NodeName]][1])||is.null(RealBioNet$OutEdge[[NodeName]][1])){
      # No input 
      Successor=NA;
    } else {
      Successor=RealBioNet[[1]][RealBioNet$OutEdge[[NodeName]]+1];
    }
  }
  Res=list(TruthTable=ttt,Successor=Successor);
  return(Res);
}

