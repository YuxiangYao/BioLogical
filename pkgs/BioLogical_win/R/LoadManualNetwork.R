#' @title Load an user-defined network
#' @description Instead of importing from files, this function allows users to 
#' define regulatory relationships and network structures to explore various 
#' system features.
#' @param Object a formatted data, record regulatory patterns (Support \code{matrix}, \code{data.frame}, \code{igraph-object})
#' @param WhichCols integer Vector, Which columns should be source and target columns (Default: \code{c(1,2)})
#' @param BoolFunList list, Boolean function, its length equals to size of system, and each element denotes a Boolean function. (Default: \code{NULL})
#' @param OBF_ratio a numeric in \code{[0,1]}, proportion of ordered  Boolean function; it will be invalid if \code{is.null(BoolFunList)} is \code{FALSE}.
#' @param OBF_type a character, type of ordered Boolean function, acceptable parameters can be found document of \link{BoolFun_Generator}.
#' @return a list contains member names of system, topological structure, and Boolean function.
#' @export
#' # Return a List of Boolean network with functions.
#'
LoadManualNetwork<-function(Object, WhichCols=c(1,2), BoolFunList=NULL,
                  OBF_ratio=0, OBF_type='R'){
  # Filter valid edges.
  if(class(Object)[1]=="matrix"||class(Object)[1]=="data.frame"){
    in.object=Object[,WhichCols];# Only two columns
    in.object=in.object[!duplicated(in.object),];# Remove repeative.
    Sources=as.character(in.object[,1]);
    Targets=as.character(in.object[,2]);
  }else if (class(Object)[1]=="igraph"){# Is an igraph object.
    in.object=igraph::as_data_frame(Object);
    in.object=in.object[!duplicated(in.object),];# Remove repeative.
    Sources=as.character(in.object[,1]);
    Targets=as.character(in.object[,2]);
  } else {
    stop("Unknown object type. Please see LoadManualNetwork() help documents.");}
  # Now, configure the list frame of network.
  AllMember=union(Sources,Targets);
  InEdge=OutEdge=BoolFun=as.list(rep(NA,length(AllMember)));
  names(InEdge)=names(OutEdge)=names(BoolFun)=AllMember;
  # Loop doing, record the networks.
  for(ii in c(1:length(Sources))){
    if(is.na(OutEdge[[Sources[ii]]][1])){# Set out-edge
      OutEdge[[Sources[ii]]]=which(Targets[ii]==AllMember)-1;}
    else {
      OutEdge[[Sources[ii]]]=c(OutEdge[[Sources[ii]]],which(Targets[ii]==AllMember)-1);}
    if(is.na(InEdge[[Targets[ii]]][1])){# Set in-edge
      InEdge[[Targets[ii]]]=which(Sources[ii]==AllMember)-1;}
    else {
      InEdge[[Targets[ii]]]=c(InEdge[[Targets[ii]]],which(Sources[ii]==AllMember)-1);}}
  # Set Boolean functions
  n.fun=length(AllMember);
  # if(is.null(RandSeed)){
  #   ran.seed=as.integer(Sys.time());}
  # else {
  #   ran.seed=RandSeed;}
  if(is.null(BoolFunList)){# Null slot for Boolean functions
    if(OBF_ratio>1||OBF_ratio<0){
      stop("Invalid input of [OBF_ratio].");}
    obf.config=sample.int(n.fun,round(n.fun*OBF_ratio));
    for(ii in c(1:n.fun)){
      if(!is.na(InEdge[[ii]])){# Has input edges
        if(ii %in% obf.config){# Should configure OBFs
          BoolFun[[ii]]=BoolFun_Generator(bf_type=OBF_type, VarNum=length(InEdge[[ii]]));}
        else {# Should configure RBFs
          BoolFun[[ii]]=BoolFun_Generator(bf_type="R", VarNum=length(InEdge[[ii]]));}}}}
  else {
    if(length(BoolFun)!=length(AllMember)){
      stop("Different lengths of size and function list.\nPlease check it.");}
    # Note: users provided data can partial match demanded lengths, our codes automatically complete or truncate.
    std.leng=c(2L,4L,8L,16L, 32L,64L,128L,256L, 512L,1024L,2048L,4096L, 8192L,16384L,32768L,65536L);
    for(ii in c(1:n.fun)){
      if(!is.na(InEdge[[ii]])){# Has inputs
        real.leng=length(BoolFunList[[ii]]);
        need.leng=std.leng[InEdge[[ii]]];
        if(real.leng>=need.leng){# Enough long for configuration
          BoolFun[[ii]]=sign(abs(BoolFunList[[ii]][1:need.leng]));}
        else {# Shorter than needed length
          tmp.bf=runif(need.leng-real.leng)>0.5;
          BoolFun[[ii]]=c(BoolFunList[[ii]],tmp.bf);}}}}
  return(list(AllMember=AllMember, InEdge=InEdge, OutEdge=OutEdge, BoolFun=BoolFun));
}
