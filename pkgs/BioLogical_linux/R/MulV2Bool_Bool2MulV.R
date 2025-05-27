#' @title Boolean to multi-valued and multi-valued to Boolean transformation.
#' @param aVec IntegerVector, represents a mapping table of multi-valued function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system (Default: 3).
#' @param Thres IntegerVector, Thresholds of multi-valued system. Its length is k+1;
#' each element is [1,L-1].
#' @param Bool2MulV a logical values, Original map is Boolean (Default: TRUE)
#' @param MappingTable logical value. Is the mapping table also returned? (Default: \code{FALSE}) 
#' @return an IntegerVector or IntegerMatrix, mapping table of target system.
#' @export
#' @examples
#' # Transform a Boolean to a multi-valued one.
#' # the Boolean is f(x,y)=x & y, the binary threshold is [0,1), [1,3), see results
#' # set.seed(2025)
#' # MulV2Bool_Bool2MulV(c(0,0,0,1), 2, 3, c(1,1,1), MappingTable = TRUE)
#' # Return
#' #       v1 v0 xx
#' #  [1,]  0  0  0
#' #  [2,]  0  1  0
#' #  [3,]  0  2  0
#' #  [4,]  1  0  0
#' #  [5,]  1  1  2
#' #  [6,]  1  2  2
#' #  [7,]  2  0  0
#' #  [8,]  2  1  1
#' #  [9,]  2  2  2
#' set.seed(2025)
#' MulV2Bool_Bool2MulV(c(0,0,0,1), 2, 3, c(1,1,1), MappingTable = TRUE)
#' 
MulV2Bool_Bool2MulV<-function(aVec, k=0, L=3, Thres=NA, Bool2MulV=TRUE, 
  MappingTable=FALSE){
  if(is.na(Thres[1])){
    thres=sample.int((L-1), as.integer(k+1), replace=TRUE);
  } else {
    if((k+1)!=length(Thres)){
      stop("Should give k+1 thresholds of multi-valued system\n");
    }
    if(any(Thres<=0)||any(Thres>=L)){
      stop("Each threshold should be [1,L-1].\n");
    }
    thres=as.integer(Thres);
  }
  if(Bool2MulV){# a Boolen 2 Multi-valued
    a_fun=r_CheckValidBoolFun(aVec)[[1]];
  } else {
    a_fun=r_CheckValidMulVFun(aVec, as.integer(k), as.integer(L));
  }
  xx=c_MulV2Bool_Bool2MulV(a_fun, as.integer(k), as.integer(L), 
    thres, Bool2MulV);
  if(MappingTable){
    if(Bool2MulV){
      xx=cbind(r_GenerateTruthTable_M(k,L), xx);
    } else {
      xx=cbind(r_GenerateTruthTable(k, FALSE), xx);
    }
  }
  return(xx);
}