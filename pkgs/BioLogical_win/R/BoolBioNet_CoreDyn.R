#' @title Analyze core dynamic of Boolean networks
#' @description This analysis only emphasizes the static feature of the network.
#' The function recursively analyzes core dynamic components (\code{coredyn}) or 
#' engaged nodes (\code{scaling}) that contribute to systematic dynamic feature. 
#' These two indicators consider coupled feedforward loop or not, respectively.
#' @param RealBioNet a 4-element formatted list, record biological network information.
#' @param AnalysisMod string, \code{coredyn} for core dynamic analysis, \code{scaling} for searching all engaged nodes.
#' @param ExternalSelfLoop logical, should convert all external inputs as self-loop nodes? (Default: \code{TRUE})
#' @param Controller integer or character vector, denote which nodes should be controlled. (Default: \code{NULL})
#' @param ControllerValue Boolean vector, controlled nodes' values. (Default: \code{NULL})
#' @param ExternalNode Boolean vector, set non input nodes' state (Default: \code{NULL}, random setting, see details)
#' @param NodeAttri logical, should return node's attribute? (Default: \code{FALSE}, ruturn value: \code{NA})
#' @param ResidualNet logical, should return residual network? (Default: \code{FALSE}, ruturn value: \code{NA})
#' @param Times an positive integer, recursive number of analysis; not recommended for users to change it.
#' @details Some nodes within bio-networks have no input, called as free-nodes.
#' If \code{ExternalSelfLoop} is \code{TRUE}, all free-nodes are added a self-loop
#' circuits. When it is \code{FALSE}, free-nodes depend on following parameters.
#' \code{ExternalNode} configure all free-nodes. If \code{ExternalNode} is \code{NULL},
#' they can be set via series random Boolean values (0/1). If provided, function only use the [1,N_FreeNode] 
#' values to set; when its length shorter than N_FreeNode, function return an error message and stop. 
#' 
#' \code{Controller} denotes manual setting of controlling nodes. Node's indexes or names are both acceptable. Note
#' that function does not check index's correctness. Please ensure the correct indexes. \code{ControllerValue} denotes 
#' corresponding values of controlled genes. If not provided (\code{NULL}), they are generated randomly. If provided,
#' only the [1,N_Control] are utilized; shorter than \code{N_Control} return an error message and stop. One gene is
#' both included in \code{Controller} and \code{ExternalNode}, the configurated value of the former covers the latter one. 
#' 
#' \code{Residual networks} still contain all nodes within original systems for comparison. Stable nodes, useless 
#' nodes, and corresponding egdes would be removed and cut-off. Remaing nodes, edges, and Boolan functions
#' are engaged in terminal dynamic behaviors.
#' 
#' @return List[4], [[1]]Overview: [1] size stable node, [2] size of invalid node, [3] size of valid nodes, 
#' [4] size of external component, [5] size of terminal component
#' [[2]]Node_S01U: Node's attribute 1(\code{IntegerVector} or \code{NA}), 0-stable (as 1), 1-stable (as 2), unstable (as 3).
#' [[3]]Node_SUE: Node's attribute 2(\code{IntegerVector} or \code{NA}), stable (as 1), useless (as 0), engaged (as -1).
#' [[4]]ResidualNetwork: residual network (\code{List[4]} or \code{NA}), still contains all nodes.
#' @export
#' @examples
#' # Analysis of core dynamic components of a Boolean genetic network.
#' # BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, 
#' #   ExternalSelfLoop=TRUE);# Treat external nodes as self-loop!
#' # Return [[1]] [1] 0, 28, 23, 0, 0
#' #        [[2]] ~ [[4]] NA
#' BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, 
#'   ExternalSelfLoop=TRUE);
#' 
#' # set.seed(1234);
#' # BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, 
#' #   ExternalSelfLoop=FALSE);# Treat free-node as fixed nodes!
#' # Return [[1]] [1] 51, 0, 0, 2, 0
#' #        [[2]] ~ [[4]] NA
#' set.seed(1234)
#' BoolBioNet_CoreDyn(BoolGRN_CellCollective$c_11863, 
#'   ExternalSelfLoop=FALSE)
#' 
BoolBioNet_CoreDyn<-function(RealBioNet, AnalysisMod=c("coredyn","scaling"), ExternalSelfLoop=TRUE, Controller=NULL, 
  ControllerValue=NULL, ExternalNode=NULL, NodeAttri=FALSE, ResidualNet=FALSE, Times=100L){
  # Prepare for transmitting parameters to C-prototype function.
  Size=length(RealBioNet[[1]]);
  ind=otd=rep(0L,Size);
  con.vals=rep(666L,Size);# Here is real controlled values.
  # Check no input nodes.
  NonIputer=NULL;# Label recording.
  if(ExternalSelfLoop){
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
    if(is.null(ControllerValue)){# Not provided, random setting.
      con.vals[Controller]=runif(length(Controller))>0.5;}
    else {# Provided.
      if(length(Controller)<=length(ControllerValue)){# Enough long
        con.vals[Controller]=as.logical(ControllerValue)[1:length(Controller)];}
      else {
        stop("Length of provided 'ControllerValue' is insufficient.");}}}
  
  # Execute analysis >>>
  if('coredyn'==AnalysisMod[1]){
    xx=c_CoreDynamicNode(arealbionet, con.id, con.vals, ind, otd, 
      as.integer(NodeAttri), as.integer(ResidualNet), as.integer(Times));
  } else if('scaling'==AnalysisMod[1]){
    xx=c_ScalingLaw_RealNet(arealbionet, con.id, con.vals, ind, otd, 
      as.integer(NodeAttri), as.integer(ResidualNet));
  } else {
    stop("Input of 'AnalysisMod' is invalid. Only 'coredyn' and 'scaling'.");
  }
  names(xx)=c("Overview","Node_S01U","Node_SUE","ResidualNetwork");
  names(xx[[1]])=c("Stable","Useless","Engaged","ExternalNode","Terminal");
  if(ResidualNet&&(xx[[1]][3]>0)){# Should return ResNet (Exist engaged nodes).
    names(xx[[4]][[2]])=xx[[4]][[1]];
    names(xx[[4]][[3]])=xx[[4]][[1]];
    names(xx[[4]][[4]])=xx[[4]][[1]];}
  return (xx);
}