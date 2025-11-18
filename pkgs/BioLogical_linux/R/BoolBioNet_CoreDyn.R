#' @title Analyze dynamic core components of Boolean networks
#' @description This analysis only emphasizes the static feature of the network.
#' The function recursively analyzes dynamic core components (\code{dyncore}) or 
#' relevant nodes (\code{relevant}) that contribute to systematic dynamic feature. 
#' These two indicators consider coupled feedforward loop or not, respectively.
#' @param RealBioNet A four-element formatted list that records information of a Boolean network.
#' @param AnalysisMod A string specifying either \code{dyncore} for dynamic core
#' analysis or \code{relevant} for identifying all relevant nodes.
#' @param ExtSelfLoop A logical value. Should inputs be converted into
#' self-loops? (Default: \code{TRUE})
#' @param Controller An integer or character vector indicating the nodes to be
#' controlled. (Default: \code{NULL})
#' @param ConVals A Boolean vector specifying the values of controlled of
#' controlled nodes. (Default: \code{NULL})
#' @param ExternalNode A Boolean vector specifying the states of non-input nodes.
#' (Default: \code{NULL}, indicating random assignment; see details for further information)
#' @param NodeAttri A logical value. Should return node's attribute? (Default: \code{FALSE},
#' return \code{NA})
#' @param ResidualNet A logical value. Should return the residual network? (Default: 
#' \code{FALSE}, return \code{NA})
#' @param Times A positive integer representing the recursive number of analyses; 
#' users are advised not to modify this value.
#' @details Some nodes in genetic networks lack inputs and are referred to as 
#' free nodes. If \code{ExtSelfLoop} is \code{TRUE}, self-loops are added 
#' to all free-nodes. When it is \code{FALSE}, the configuration of free nodes 
#' depends on \code{ExternalNode}. If \code{ExternalNode} is \code{NULL}, the 
#' states of the free nodes are fixed at randomly generated Boolean values (0 or
#' 1). If provided (a Boolean vector), only the first \code{free node} elements 
#' are used; if the length of \code{ExternalNode} is less than \code{free node},
#' the function returns an error message and execution halts.
#' 
#' The \code{Controller} parameter specifies the manually controlled nodes, 
#' where either node indices or names are acceptable. Note that the function does
#' not validate the correctness of the indices; users must ensure their accuracy.
#' \code{ConVals} provides the corresponding values for the controlled 
#' genes. If \code{ConVals} is not specified (\code{NULL}), values are 
#' generated randomly. If provided (a Boolean vector), only the first 
#' \code{|{Controller}|} values are used; if the length of \code{ConVals} 
#' is less than that of \code{Controller}, the function returns an error message 
#' and execution halts. If a gene/node appears in both \code{Controller} and 
#' \code{ExternalNode}, the value specified in \code{Controller} takes precedence
#' over that in \code{ExternalNode}.
#' 
#' \code{ResidualNet} retain all nodes of the original systems for 
#' comparative purposes (if \code{ResidualNet} is \code{TRUE}). Stable nodes, 
#' useless nodes, and their corresponding edges are removed in subsequent analyses.
#' Only the remaining nodes, edges, and Boolean functions contribute to the 
#' terminal dynamic behaviors.
#' 
#' @return A list containing four elements: 
#' \itemize{
#' \item \code{[[1]] Overview}: Sizes of [1] stable nodes, [2] invalid nodes, [3] 
#' valid nodes, [4] external components, and [5] terminal component.
#' \item \code{[[2]] Node_S01U}: Node attribute 1 (\code{IntegerVector} or 
#' \code{NA}). Where "1", "2", "3" represent the node being stable at state 0, 
#' stable at state 1, and unstable, respectively.
#' \item \code{[[3]] Node_SUE}: Node attribute 2 (\code{IntegerVector} or 
#' \code{NA}). Where "1", "0", "-1" represent that the node belongs to stable, 
#' useless, and engage types, respectively.
#' \item \code{[[4]] ResidualNetwork}: The residual network is the analyzed 
#' network containing all nodes but with all invalid edges removed for 
#' comparative analysis.}
#' @export
#' @examples
#' # Analysis of dynamic core components of a Boolean genetic network.
#' # Treat external nodes as self-loop!
#' # BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, ExtSelfLoop=TRUE);
#' # Return [[1]] [1] 0, 28, 23, 0, 0
#' #        [[2]] ~ [[4]] NA
#' BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, ExtSelfLoop=TRUE)
#' 
#' # set.seed(1234);
#' # Treat free-node as fixed nodes!
#' # BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, ExtSelfLoop=FALSE);
#' # Return [[1]] [1] 51, 0, 0, 2, 0
#' #        [[2]] ~ [[4]] NA
#' 
#' set.seed(1234);
#' BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, ExtSelfLoop=FALSE)
#' 

BoolBioNet_CoreDyn<-function(RealBioNet, AnalysisMod=c("dyncore","relevant"), ExtSelfLoop=TRUE, Controller=NULL, 
  ConVals=NULL, ExternalNode=NULL, NodeAttri=FALSE, ResidualNet=FALSE, Times=100L){
  # Prepare for transmitting parameters to C-prototype function.
  Size=length(RealBioNet[[1]]);
  ind=otd=rep(0L,Size);
  con.vals=rep(666L,Size);# Here is real controlled values.
  # Check no input nodes.
  NonIputer=NULL;# Label recording.
  if(ExtSelfLoop){
    arealbionet=r_UnfreezeInputNode(RealBioNet);
  } else {
    arealbionet=RealBioNet;
  }
  for(ii in c(1:Size)){
    if(is.na(arealbionet[[2]][[ii]][1])||is.null(arealbionet[[2]][[ii]][1])){
      ind[ii]=0L;
      NonIputer=c(NonIputer,ii);}# Insert Exponents
    else {
      ind[ii]=length(arealbionet[[2]][[ii]]);}
    if(is.na(arealbionet[[3]][[ii]][1])||is.null(arealbionet[[3]][[ii]][1])){
      otd[ii]=0L;}
    else {
      otd[ii]=length(arealbionet[[3]][[ii]]);}}
  
  # Check external nodes and Controller settings.
  # First, external node
  if(!is.null(NonIputer)){# Has external nodes. Note that 
    if(is.null(ExternalNode)){# Not provided, random setting.
      con.vals[NonIputer]=runif(length(NonIputer))>0.5;}
    else {# Provided.
      if(length(NonIputer)<=length(ExternalNode)){# Enough long
        con.vals[NonIputer]=as.logical(ExternalNode)[1:length(NonIputer)];}
      else {
        stop("Length of provided 'ExternalNode' is insufficient.");}}}
  
  # Second, Controlled genes.
  con.id=-666;
  if(!is.null(Controller)){# Has controllers.
    con.id=Controller-1;# Note the differences in the index between C++ and R.
    if(is.null(ConVals)){# Not provided, random setting.
      con.vals[Controller]=runif(length(Controller))>0.5;}
    else {# Provided.
      if(length(Controller)<=length(ConVals)){# Enough long
        con.vals[Controller]=as.logical(ConVals)[1:length(Controller)];}
      else {
        stop("Length of provided 'ConVals' is insufficient.");}}}
  
  # Execute analysis >>>
  if('dyncore'==AnalysisMod[1]){
    xx=c_CoreDynamicNode(arealbionet, con.id, con.vals, ind, otd, 
      as.integer(NodeAttri), as.integer(ResidualNet), as.integer(Times));
  } else if('relevant'==AnalysisMod[1]){
    xx=c_ScalingLaw_RealNet(arealbionet, con.id, con.vals, ind, otd, 
      as.integer(NodeAttri), as.integer(ResidualNet));
  } else {
    stop("Input of 'AnalysisMod' is invalid. Only 'coredyn' and 'scaling'.");
  }
  names(xx)=c("Overview","Node_S01U","Node_SUE","ResidualNetwork");
  if('dyncore'==AnalysisMod[1]){
    names(xx[[1]])=c("Stable", "Relay", "DynCore", "ExternalNode","Terminal");
  } else if('relevant'==AnalysisMod[1]){
    names(xx[[1]])=c("Stable", "Relay", "Relevant", "ExternalNode","Terminal");
  } else {
    ;
  }
  if(ResidualNet&&(xx[[1]][3]>0)){# Should return ResNet (Exist engaged nodes).
    names(xx[[4]][[2]])=xx[[4]][[1]];
    names(xx[[4]][[3]])=xx[[4]][[1]];
    names(xx[[4]][[4]])=xx[[4]][[1]];}
  return (xx);
}