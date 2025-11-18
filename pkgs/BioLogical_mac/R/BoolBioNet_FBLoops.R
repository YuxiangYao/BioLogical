#' @title Calculate the network's feedback loop
#' @description This function provides a brief analysis of feedback loops within 
#' the systems.
#' @param OneNet A four-element formatted list that records information of a Boolean network.
#' @param AnalyType A string indicating the feature of the feedback loop being analyzed.
#' Only those "real" strongly connected components were used for analysis. 
#' (See \link{BoolBioNet_StrConComp})
#' \itemize{
#'   \item \code{loop}: Return the list of feedback loops. 
#'   \item \code{length}: Count the distribution of feedback loop lengths.
#'   \item \code{edgeattr}: Analyze the fundamental properties of each edge in the
#' loop, including the number of embedded loops, sign, edge activity, effectiveness.
#'   \item \code{loopsign}: Return the sign property of each cycle.
#' }
#' The sign attributes are categorized into three types: positive (+), negative 
#' (-), and hybrid (?). For the definitions of sign, edge activity, and effectiveness,
#' refer to the following papers: 
#' [\href{https://doi.org/10.1007/s11538-008-9304-7}{Aracena2008}],
#' [\href{https://doi.org/10.1073/pnas.0407783101}{Kauffman2004}],
#' [\href{https://doi.org/10.1073/pnas.2022598118}{Gates2021}].
#' @param Isolated A logical value. Should keep self-loops? (Default: \code{TRUE})
#' @return Depending on the value of \code{AnalyType}, the following types of 
#' return values may result:
#' \itemize{
#'   \item \code{loop}: a list, with the names of the elements indicating the 
#' SCC composition and the feedback loop.
#'   \item \code{length}: a matrix storing the distribution of feedback loop lengths.
#'   \item \code{edgeattr}: a datamframe characterizing the properties of each 
#' edge in the loop.
#'   \item \code{loopsign}: a matrix recording the sign attributes of each edge 
#' within each cycle.
#' }
#' @export
#' @examples
#' # Analyze various properties of feedback loops in gene networks.
#' # Here is an example of the network [k_1753781] from Kadelka's paper.
#' # BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "loop")
#' # [[1]]scc1_loop1
#' # [1] "A"  "T2" "T1" "B"  "A" 
#' # [[2]]scc1_loop2
#' # [1] "A"  "T2" "T1" "H"  "B"  "A" 
#' # ... ... 
#' 
#' # BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "length")
#' #      Loop_Length Count
#' # [1,]           1     6
#' # [2,]           2     8
#' # ... ... 
#' 
#' # BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "edgeattr")
#' #    from to support SignType activity effectiveness
#' # 1     A T1       4        - 0.015625    0.14955357
#' # 2     A T2       8        + 0.031250    0.09505208
#' # ... ... 
#' 
#' # BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "loopsign")
#' #      Positive Negative Hybrid
#' # [1,]        3        1      0
#' # [2,]        4        1      0
#' # ... ... 
#'
#' BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "loop")
#' BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "length")
#' BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "edgeattr")
#' BoolBioNet_FBLoops(BoolGRN_KadelkaSet$k_1753781, "loopsign")
#' 
BoolBioNet_FBLoops<-function(OneNet, 
    AnalyType=c("loop", "length", "edgeattr","loopsign"),
    Isolated=TRUE){
  if(!(AnalyType[1]%in%c("loop", "length", "edgeattr","loopsign"))){
    stop("Input of 'AnalyType' is invalid. Please check the argument.");}
  Size=length(OneNet[[1]]);
  ind=otd=rep(0L,Size);
  con.vals=rep(666L,Size);# Here is real controlled values.
  NonIputer=NULL;# Label recording.
  for(ii in c(1:Size)){
    if(is.na(OneNet[[2]][[ii]][1])||is.null(OneNet[[2]][[ii]][1])){
      ind[ii]=0L;
      NonIputer=c(NonIputer,ii);}# Insert Exponents
    else {
      ind[ii]=length(OneNet[[2]][[ii]]);}
    if(is.na(OneNet[[3]][[ii]][1])||is.null(OneNet[[3]][[ii]][1])){
      otd[ii]=0L;}
    else {
      otd[ii]=length(OneNet[[3]][[ii]]);}}
  # Main code >>>
  if("loop"==AnalyType[1]){# All FBLs 
    tmp0=c_ObtainFeedBackLoop(OneNet,ind,otd,Isolated,1);
    #remove_na_from_second_layer <- function(x){ lapply(x, function(sublist)
    #   {Filter(Negate(is.null), sublist)[!sapply(sublist, function(y) is.atomic(y) && is.na(y))]})}
    #tmp0=remove_na_from_second_layer(tmp1);
    tmp=list();
    for(ii in c(1:length(tmp0))){
      tmpx=tmp0[[ii]];# A SCC's FBLs.
      for(jj in c(1:length(tmpx))){
        tmpx[[jj]]=OneNet[[1]][tmp0[[ii]][[jj]]+1];}
      names(tmpx)=paste0("scc",ii,"_loop",c(1:length(tmpx)));
      tmp=c(tmp,tmpx);
    }
  } else if("length"==AnalyType[1]){# Loop's length
    tmp=c_ObtainFeedBackLoop(OneNet,ind,otd,Isolated,2)[[1]];
    colnames(tmp)=c("Loop_Length","Count");
  } else if("edgeattr"==AnalyType[1]){# Edge's attributes.
    tmp=c_ObtainFeedBackLoop(OneNet,ind,otd,Isolated,3)[[1]];
    tmp=cbind.data.frame(tmp[[1]],tmp[[2]]);
    tmp$from=OneNet[[1]][tmp$from+1];
    tmp$to=OneNet[[1]][tmp$to+1];
    tmp$SignType=as.character(factor(tmp$SignType, 
      levels=0:2, labels=c("+","-","?")));# +,
  } else if("loopsign"==AnalyType[1]){# Loop+Edge's sign
    tmp=c_ObtainFeedBackLoop(OneNet,ind,otd,Isolated,4)[[1]];
    rownames(tmp)=paste0("loop_",c(1:nrow(tmp)));
  } else {
    ;}
  return (tmp);
}