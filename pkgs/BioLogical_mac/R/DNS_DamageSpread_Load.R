#' @title Analyze damage spread within system (Loaded systems)
#' @description Same as \link{DNS_DamageSpread}, but the network is loaded. 
#' @param RealBioNet A four-element formatted list that records information of a Boolean network.
#' @param ExtSelfLoop A logical value. Should inputs be converted into
#' self-loops? (Default: \code{TRUE})
#' @param Controller An integer or character vector indicating the nodes to be
#' controlled. (Default: \code{NULL})
#' @param ConVals A Boolean vector specifying the values of controlled of
#' controlled nodes. (Default: \code{NULL})
#' @param ExternalNode A Boolean vector specifying the states of non-input nodes.
#' (Default: \code{NULL}, indicating random assignment)
#' @param Times an integer, steps of simulating system.
#' @param Init_Dist a numeric value, initial normalized Hamming distance of two vector.
#' @param Init_1_Ratio a NumvericVector, proportion of each value (its length is 
#' \code{NumSys}).
#' @param NumSys an integer, number of discrete value.
#' @param UpdateRule an integer, update rule: 1=synchronous, 2=asynchronous, 
#' 3=quick-asynchronous.
#' @details See \link{DNS_DamageSpread}.
#' @return A numeric value that representing the normalized Hamming distance.
#' @export
#' @examples
#' # set.seed(20250101L);
#' # DNS_DamageSpread(Size=1000L, SimStep=1000L, OBF_Type='C',
#' #   OBF_iPara1=1, OBF_iPara2=-1, OBF_Ratio=0.5);
#' # Return 0.335
#' # set.seed(20250101L)
#' # DNS_DamageSpread(Size=1000L, SimStep=1000L)
#' 
DNS_DamageSpread_Load<-function(RealBioNet, ExtSelfLoop=TRUE, Controller=NULL, 
  ConVals=NULL, ExternalNode=NULL, Times=1000L, Init_Dist=0.10,
  Init_1_Ratio=c(0.5,0.5), NumSys=2L, UpdateRule=1L){
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
  xx=c_Derrida_Load(arealbionet, con.id, con.vals, ind, otd, Init_Dist,
    as.numeric(Init_1_Ratio), as.integer(Times), as.integer(UpdateRule) );
  return (xx);
}
